#!/usr/bin/env python3
"""Audit a completed bpRNA chunk run; retain small, reproducible audit artifacts.

Uses only the Python standard library. Does not modify production outputs.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import subprocess

LENS = (1, 3, 7, 10)
EXPECTED = {
    'engine': 'lincapr', 'effective_energy_model': 'turner2004',
    'requested_energy_model': 'turner2004', 'beam': '100',
    'c_multi': '30', 'c_hairpin': '30', 'access_lens': '1,3,7,10',
    'output_mode': 'probability', 'normalize_profiles': 'true',
}


def fasta(path):
    result = []
    name, seq = None, []
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    result.append((name, ''.join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        result.append((name, ''.join(seq)))
    return result


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def distribution(values):
    if not values:
        return None
    values = sorted(values)
    def quantile(p):
        x = (len(values) - 1) * p
        lo, hi = math.floor(x), math.ceil(x)
        return values[lo] + (values[hi] - values[lo]) * (x - lo)
    return dict(n=len(values), min=values[0], median=statistics.median(values),
                mean=statistics.mean(values), p90=quantile(.9),
                p95=quantile(.95), p99=quantile(.99), max=values[-1],
                sum=sum(values))


def audit_output(path, seq_id, length):
    errors = Counter()
    metadata, headers = {}, []
    pos_count = 0
    stats = {k: dict(count=0, zero=0, one=0, negative=0, above_one=0,
                     outside_tolerance_1e_12=0, outside_tolerance_1e_6=0,
                     min=None, max=None) for k in LENS}
    if not path.is_file():
        return Counter(missing_output=1), metadata, stats
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if '=' in line:
                    key, value = line[1:].strip().split('=', 1)
                    metadata[key] = value
                continue
            if line.startswith('>'):
                headers.append(line[1:].split()[0])
                continue
            try:
                pos, payload = line.split('\t')
                pos = int(pos)
                pairs = [item.split(',') for item in payload.rstrip(';').split(';')]
                keys = [int(k) for k, _ in pairs]
                values = [float(v) for _, v in pairs]
            except (ValueError, TypeError):
                errors['malformed_row'] += 1
                continue
            if pos != pos_count:
                errors['position_order'] += 1
            if keys != [k for k in LENS if pos + k <= length]:
                errors['window_coverage'] += 1
            pos_count += 1
            for k, v in zip(keys, values):
                if not math.isfinite(v):
                    errors['nonfinite_probability'] += 1
                    continue
                if not 0 <= v <= 1:
                    errors['probability_out_of_range'] += 1
                if k in stats:
                    s = stats[k]
                    s['count'] += 1
                    s['zero'] += v == 0
                    s['one'] += v == 1
                    s['negative'] += v < 0
                    s['above_one'] += v > 1
                    s['outside_tolerance_1e_12'] += v < -1e-12 or v > 1 + 1e-12
                    s['outside_tolerance_1e_6'] += v < -1e-6 or v > 1 + 1e-6
                    s['min'] = v if s['min'] is None else min(s['min'], v)
                    s['max'] = v if s['max'] is None else max(s['max'], v)
    if pos_count != length:
        errors['position_count'] += 1
    if headers != [seq_id]:
        errors['sequence_header'] += 1
    for k, v in EXPECTED.items():
        if metadata.get(k) != v:
            errors['metadata_' + k] += 1
    return errors, metadata, stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--job-id', default='8661385')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--accounting-file', type=Path,
                        help='Use an existing qacct -j JOB_ID capture instead of querying UGE')
    args = parser.parse_args()
    root, out = args.repo.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    inputs = root / 'benchmark/data'
    run = root / 'benchmark/results/bprna_final'
    clean_path = inputs / 'bpRNA_1m_90_filtered_clean.fa'
    records = fasta(clean_path)
    source_path = inputs / 'bpRNA_1m_90.fasta'
    source = fasta(source_path)
    normalized = [(n, s.upper().replace('T', 'U')) for n, s in source]
    filtered = [(n, s) for n, s in normalized if len(s) >= 50 and re.fullmatch('[ACGU]+', s)]
    issues, metadata_counts, families = Counter(), {}, Counter()
    if filtered != records:
        issues['filter_reproduction_mismatch'] += 1
    if len(set(n for n, _ in records)) != len(records):
        issues['duplicate_input_ids'] += 1
    chunks = sorted((inputs / 'bprna_chunks').glob('chunk_*.fa'))
    combined = [record for p in chunks for record in fasta(p)]
    if combined != records:
        issues['chunk_content_or_order_mismatch'] += 1
    prob_stats = {k: dict(count=0, zero=0, one=0, negative=0, above_one=0,
                          outside_tolerance_1e_12=0, outside_tolerance_1e_6=0,
                          min=None, max=None) for k in LENS}
    all_rows, issue_rows, lengths, times = [], [], [], []
    output_hash = hashlib.sha256()
    log_hash = hashlib.sha256()
    warn_count = warn_sequences = other_log_lines = 0
    warn_min, warn_max = math.inf, -math.inf
    artifacts = []
    def artifact(p):
        return dict(path=str(p.relative_to(root)), bytes=p.stat().st_size, sha256=sha256(p))
    for chunk in chunks:
        task = int(chunk.stem.split('_')[1])
        batch = run / f'batch_{task}'
        summary = batch / 'summary.tsv'
        artifacts.extend([artifact(chunk), artifact(summary)])
        with summary.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        counts = Counter(r['seq_id'] for r in rows)
        chunk_records = fasta(chunk)
        if counts != Counter(n for n, _ in chunk_records):
            issues['summary_id_coverage'] += 1
        by_id = {r['seq_id']: r for r in rows}
        for number, (name, seq) in enumerate(chunk_records, 1):
            err = Counter()
            row = by_id.get(name)
            elapsed = None
            if row:
                if row['status'] != 'SUCCESS':
                    err['unsuccessful_status'] += 1
                if (int(row['length']), int(row['chunk_id']), int(row['seq_num'])) != (len(seq), task, number):
                    err['summary_input_mismatch'] += 1
                elapsed = float(row['elapsed_sec'])
                if not math.isfinite(elapsed) or elapsed < 0:
                    err['invalid_elapsed'] += 1
                    elapsed = None
            else:
                err['missing_summary_row'] += 1
            safe = name.replace('/', '_')
            p = batch / (safe + '.txt')
            output_errors, meta, stats = audit_output(p, name, len(seq))
            err.update(output_errors)
            if p.is_file():
                output_hash.update(f'{p.relative_to(root)}\t{sha256(p)}\n'.encode())
            for k, v in meta.items():
                metadata_counts.setdefault(k, Counter())[v] += 1
            for k, s in stats.items():
                agg = prob_stats[k]
                for field in ('count', 'zero', 'one', 'negative', 'above_one',
                              'outside_tolerance_1e_12', 'outside_tolerance_1e_6'):
                    agg[field] += s[field]
                for field, fn in [('min', min), ('max', max)]:
                    if s[field] is not None:
                        agg[field] = s[field] if agg[field] is None else fn(agg[field], s[field])
            warnings = 0
            log = run / 'logs' / (safe + '.log')
            if log.is_file():
                log_hash.update(f'{log.relative_to(root)}\t{sha256(log)}\n'.encode())
                for line in log.read_text().splitlines():
                    match = re.fullmatch(r'warn: prob_sum\[\d+\]=(.+)', line)
                    if match:
                        warnings += 1
                        v = float(match[1])
                        if math.isfinite(v):
                            warn_min, warn_max = min(warn_min, v), max(warn_max, v)
                        else:
                            err['nonfinite_prob_sum_warning'] += 1
                    elif line.strip():
                        other_log_lines += 1
            else:
                err['missing_sequence_log'] += 1
            warn_count += warnings
            warn_sequences += warnings > 0
            family = name.split('_')[1]
            families[family] += 1
            lengths.append(len(seq))
            if elapsed is not None:
                times.append(elapsed)
            all_rows.append([task, number, name, family, len(seq), row['status'] if row else 'MISSING', elapsed, warnings, sum(err.values())])
            for kind, count in err.items():
                issue_rows.append([name, kind, count])
            issues.update(err)
        print(f'Audited chunk {task}/{len(chunks)}', flush=True)
    if args.accounting_file:
        accounting_text = args.accounting_file.read_text()
        accounting_returncode = 0
    else:
        accounting = subprocess.run(['qacct', '-j', args.job_id], capture_output=True, text=True)
        accounting_text = accounting.stdout + accounting.stderr
        accounting_returncode = accounting.returncode
    (out / 'qacct.txt').write_text(accounting_text)
    accounting_records = []
    for block in re.split(r'^=+\s*$', accounting_text, flags=re.MULTILINE):
        fields = dict(line.split(None, 1) for line in block.splitlines() if len(line.split(None, 1)) == 2)
        fields = {key: value.strip() for key, value in fields.items()}
        if 'taskid' in fields:
            accounting_records.append(fields)
    if accounting_returncode or Counter(r['taskid'] for r in accounting_records) != Counter(str(int(p.stem.split('_')[1])) for p in chunks):
        issues['accounting_coverage'] += 1
    for r in accounting_records:
        if r.get('jobnumber') != args.job_id:
            issues['accounting_job_id_mismatch'] += 1
        if r.get('failed') != '0' or r.get('exit_status') != '0':
            issues['accounting_failure'] += 1
    for filename, header, rows in [
        ('per_sequence.tsv', ['chunk_id', 'seq_num', 'seq_id', 'family', 'length_nt', 'status', 'elapsed_sec', 'prob_sum_warnings', 'audit_errors'], all_rows),
        ('issues.tsv', ['seq_id', 'issue', 'count'], issue_rows),
    ]:
        with (out / filename).open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(header)
            writer.writerows(rows)
    summary = dict(
        schema_version=1, job_id=args.job_id, sequences=len(records),
        validation_passed=not issues, issues=dict(issues),
        lengths_nt=distribution(lengths), elapsed_sec=distribution(times),
        elapsed_zero_count=times.count(0), families=dict(families), probabilities=prob_stats,
        output_metadata=metadata_counts,
        log_warnings=dict(prob_sum_count=warn_count, sequences=warn_sequences,
                          min=None if warn_min == math.inf else warn_min,
                          max=None if warn_max == -math.inf else warn_max,
                          other_nonempty_lines=other_log_lines),
        accounting=dict(tasks=len(accounting_records),
                        failed=dict(Counter(r.get('failed') for r in accounting_records)),
                        exit_status=dict(Counter(r.get('exit_status') for r in accounting_records)),
                        maxvmem_raw=dict(Counter(r.get('maxvmem', 'unavailable') for r in accounting_records))),
        limitations=[
            'No reference comparison: range and coverage checks do not establish numerical accuracy.',
            'elapsed_sec uses integer wall-clock timestamps during 10-way concurrency; zero means below timestamp resolution, not zero cost.',
            'Sum of per-sequence elapsed times is neither job wall time nor CPU time.',
            'No per-sequence peak RSS instrumentation. qacct maxvmem is task-level virtual memory, not per-sequence RSS.',
            'prob_sum warnings concern structural-profile normalization in the current source; they need separate numerical investigation.',
            'Output source_revision is unknown; current executable/source hashes do not prove the historical executable identity.',
        ])
    dump = lambda p, data: p.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')
    dump(out / 'summary.json', summary)
    binary = root / 'vendor/LinearRaccess/build-raccess/LinRacc'
    provenance = dict(
        schema_version=1, experiment_id='bprna-production-8661385-retrospective-v1',
        purpose='Completion, probability validity, length and runtime audit; no comparator',
        observed_at_utc=datetime.now(timezone.utc).isoformat(), job_id=args.job_id,
        inputs=[artifact(source_path), artifact(clean_path)], chunk_and_summary_files=artifacts,
        filtering=dict(min_length=50, uppercase=True, T_to_U=True, alphabet='ACGU',
                       source_sequences=len(source), retained=len(filtered),
                       too_short=sum(len(s) < 50 for _, s in normalized),
                       invalid_alphabet_at_least_50=sum(len(s) >= 50 and not re.fullmatch('[ACGU]+', s) for _, s in normalized),
                       reproduced_exact_records_and_order=filtered == records),
        parameters=EXPECTED, repetitions=1, warmup='not implemented in job script',
        resources=dict(cpu_16=1, h_rt='04:00:00', array='1-27', task_concurrency=10, sequence_concurrency=10),
        current_files=[artifact(root / p) for p in ['benchmark/scripts/audit_bprna.py', 'benchmark/scripts/uge/bprna_chunk_job.sh', 'benchmark/scripts/filter_bprna_clean.py', 'benchmark/scripts/split_bprna.py']] + [artifact(binary)],
        current_source_revision=subprocess.check_output(['git', '-C', str(root / 'vendor/LinearRaccess'), 'rev-parse', 'HEAD'], text=True).strip(),
        current_source_status=subprocess.check_output(['git', '-C', str(root / 'vendor/LinearRaccess'), 'status', '--porcelain'], text=True),
        historical_source_revision='unknown (output metadata)',
        provenance_scope='Retrospective observations, not an execution-time frozen manifest. Current files may differ from those at submission.',
        output_digest=dict(sha256=output_hash.hexdigest(), algorithm='SHA256 of path<TAB>file_sha256<LF>, in chunk and input record order'),
        sequence_log_digest=dict(sha256=log_hash.hexdigest(), algorithm='same as output_digest'),
        audit_artifacts={p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()},
        comparator=None, archival_identifier=None,
    )
    dump(out / 'manifest.json', provenance)
    print(json.dumps(dict(validation_passed=not issues, issues=dict(issues), sequences=len(records))))
    return 0 if not issues else 1


if __name__ == '__main__':
    raise SystemExit(main())
