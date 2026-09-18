#!/usr/bin/env python3
"""Compare fast/exact logsumexp on deterministic bpRNA audit examples."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import time

from audit_bprna import fasta, sha256


def read_values(path):
    values, metadata = {}, {}
    for line in path.read_text().splitlines():
        if line.startswith('#') and '=' in line:
            k, v = line[1:].strip().split('=', 1)
            metadata[k] = v
        elif line and not line.startswith(('#', '>')):
            pos, payload = line.split('\t')
            for item in payload.rstrip(';').split(';'):
                k, v = item.split(',')
                values[int(pos), int(k)] = float(v)
    return values, metadata


def describe(values):
    data = list(values.values())
    return dict(count=len(data), min=min(data), max=max(data),
                negative=sum(v < 0 for v in data), above_one=sum(v > 1 for v in data),
                outside_1e_12=sum(v < -1e-12 or v > 1 + 1e-12 for v in data),
                outside_1e_6=sum(v < -1e-6 or v > 1 + 1e-6 for v in data),
                nonfinite=sum(not math.isfinite(v) for v in data))


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path,
                        default=root / 'benchmark/results/bprna_logsumexp_diagnostic')
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(exist_ok=False)
    raw = out / 'raw'
    raw.mkdir()
    audit = root / 'benchmark/results/bprna_audit_8661385'
    with (audit / 'per_sequence.tsv').open() as h:
        rows = list(csv.DictReader(h, delimiter='\t'))
    worst_max, worst_min, shortest = (-math.inf, None), (math.inf, None), (math.inf, None)
    selected = {}
    for row in rows:
        if not int(row['audit_errors']):
            continue
        name = row['seq_id']
        path = root / f"benchmark/results/bprna_final/batch_{row['chunk_id']}/{name}.txt"
        values, _ = read_values(path)
        lo, hi = min(values.values()), max(values.values())
        if hi > worst_max[0]:
            worst_max = hi, name
        if lo < worst_min[0]:
            worst_min = lo, name
        if hi > 1 + 1e-6 and int(row['length_nt']) < shortest[0]:
            shortest = int(row['length_nt']), name
    for reason, (_, name) in [('largest_probability', worst_max), ('most_negative_probability', worst_min), ('shortest_overshoot_above_1e_6', shortest)]:
        selected.setdefault(name, []).append(reason)
    control = next(r['seq_id'] for r in rows if not int(r['audit_errors']) and int(r['length_nt']) <= 100)
    selected.setdefault(control, []).append('first_clean_control_at_most_100nt')
    seqs = dict(fasta(root / 'benchmark/data/bpRNA_1m_90_filtered_clean.fa'))
    by_id = {r['seq_id']: r for r in rows}
    binary = root / 'vendor/LinearRaccess/build-raccess/LinRacc'
    runs = []
    for name, reasons in selected.items():
        inp = out / (name + '.fa')
        inp.write_text(f'>{name}\n{seqs[name]}\n')
        baseline_path = root / f"benchmark/results/bprna_final/batch_{by_id[name]['chunk_id']}/{name}.txt"
        baseline, _ = read_values(baseline_path)
        variants = [(100, False), (100, True)]
        if len(seqs[name]) <= 200:
            variants += [(0, False), (0, True)]
        results = {}
        for beam, exact in variants:
            label = f"beam{beam}_{'exact' if exact else 'fast'}"
            path = raw / f'{name}_{label}.txt'
            cmd = [str(binary), f'-seqfile={inp}', f'-outfile={path}',
                   '-access_len=1,3,7,10', '-engine=lincapr', '-energy=turner2004',
                   f'-beam={beam}', '-c_multi=30', '-c_hairpin=30', '-probabilities', '-metadata']
            if exact:
                cmd.append('-no-fast-logsumexp')
            start = time.perf_counter()
            with (raw / f'{name}_{label}.log').open('w') as log:
                proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=180)
            elapsed = time.perf_counter() - start
            if proc.returncode:
                raise RuntimeError(f'{label} {name} failed: {proc.returncode}')
            vals, meta = read_values(path)
            if vals.keys() != baseline.keys():
                raise RuntimeError('Window coverage mismatch')
            if meta['logsumexp'] != ('exact' if exact else 'fast'):
                raise RuntimeError('Logsumexp mode mismatch')
            results[label] = vals
            runs.append(dict(seq_id=name, length=len(seqs[name]), selection=reasons,
                             variant=label, stats=describe(vals), elapsed_sec=elapsed,
                             prob_sum_warning_lines=sum(line.startswith('warn: prob_sum[') for line in (raw / f'{name}_{label}.log').read_text().splitlines()),
                             max_abs_difference_from_production=max(abs(vals[k] - baseline[k]) for k in vals),
                             metadata=meta, output_sha256=sha256(path),
                             command=[str(Path(arg).relative_to(root)) if arg.startswith(str(root)) else arg.replace(str(root) + '/', '') for arg in cmd]))
            print(name, label, json.dumps(describe(vals)), flush=True)
        for run in runs:
            if run['seq_id'] == name and run['variant'].endswith('exact'):
                vals = results[run['variant']]
                fast = results[run['variant'].replace('exact', 'fast')]
                run['max_abs_difference_from_same_beam_fast'] = max(abs(vals[k] - fast[k]) for k in vals)
    report = dict(schema_version=1, selection=selected, runs=runs,
                  binary_sha256=sha256(binary), script_sha256=sha256(Path(__file__)),
                  inputs={p.name: sha256(p) for p in out.glob('*.fa')},
                  source_revision=subprocess.check_output(['git', '-C', str(root / 'vendor/LinearRaccess'), 'rev-parse', 'HEAD'], text=True).strip(),
                  limitations=['Selected examples, not a full bpRNA rerun.',
                               'beam=0 disables beam pruning, but loop caps remain 30.',
                               'Switching logsumexp at beam=100 may change retained beam states.',
                               'Single-run elapsed times are diagnostics, not performance benchmarks.'])
    (out / 'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
