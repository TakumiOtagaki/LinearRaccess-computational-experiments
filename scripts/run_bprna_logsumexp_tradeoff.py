#!/usr/bin/env python3
"""Prepare/run a stratified, paired fast versus standard logsumexp experiment."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import time

from audit_bprna import fasta, distribution, sha256

ROOT = Path(__file__).resolve().parents[2]
LENS = (1, 3, 7, 10)
BINS = [(50, 100), (100, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 5000)]


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')


def load(path):
    return json.loads(path.read_text())


def parse_output(path, names, length, mode):
    """Validate every repeated record and return one record's probability map."""
    records, meta, current = {}, {}, None
    for line in path.read_text().splitlines():
        if not line:
            continue
        if line.startswith('#'):
            if '=' in line:
                key, value = line[1:].strip().split('=', 1)
                meta[key] = value
        elif line.startswith('>'):
            current = line[1:]
            if current in records:
                raise ValueError('Duplicate output record')
            records[current] = {}
        else:
            pos, payload = line.split('\t')
            for item in payload.rstrip(';').split(';'):
                width, value = item.split(',')
                key = int(pos), int(width)
                value = float(value)
                if current is None or key in records[current] or not math.isfinite(value):
                    raise ValueError('Invalid or nonfinite probability output')
                records[current][key] = value
    expected = {(i, k) for k in LENS for i in range(length - k + 1)}
    if list(records) != names or meta.get('logsumexp') != mode:
        raise ValueError('Record identity or mode mismatch')
    for key, value in {'beam': '100', 'c_multi': '30', 'c_hairpin': '30',
                       'engine': 'lincapr', 'effective_energy_model': 'turner2004',
                       'output_mode': 'probability', 'access_lens': '1,3,7,10'}.items():
        if meta.get(key) != value:
            raise ValueError(f'Metadata mismatch: {key}')
    first = records[names[0]]
    for record in records.values():
        if record.keys() != expected or record != first:
            raise ValueError('Window coverage or repeated-record reproducibility mismatch')
    return first


def bounds(values):
    values = list(values)
    return dict(count=len(values), min=min(values), max=max(values),
                negative=sum(v < 0 for v in values), above_one=sum(v > 1 for v in values),
                outside_1e_12=sum(v < -1e-12 or v > 1 + 1e-12 for v in values),
                outside_1e_6=sum(v < -1e-6 or v > 1 + 1e-6 for v in values))


def physical_cpus():
    selected, seen = [], set()
    for cpu in sorted(os.sched_getaffinity(0)):
        path = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list')
        siblings = path.read_text().strip() if path.exists() else str(cpu)
        if siblings not in seen:
            seen.add(siblings)
            selected.append(cpu)
    return selected


def prepare(out):
    out.mkdir(parents=True, exist_ok=False)
    (out / 'raw').mkdir()
    dataset = ROOT / 'benchmark/data/bpRNA_1m_90_filtered_clean.fa'
    records = fasta(dataset)
    seqs = dict(records)
    rng = random.Random(20260914)
    selected = []
    for lo, hi in BINS:
        candidates = sorted(n for n, seq in records if lo <= len(seq) < hi)
        for name in rng.sample(candidates, 4):
            selected.append(dict(seq_id=name, length=len(seqs[name]),
                                 stratum=f'{lo}-{hi-1}', cohort='random'))
    for name in ['bpRNA_RFAM_3471', 'bpRNA_RFAM_42909', 'bpRNA_CRW_38238']:
        if name not in {r['seq_id'] for r in selected}:
            selected.append(dict(seq_id=name, length=len(seqs[name]), stratum='sentinel', cohort='sentinel'))
    for row in selected:
        row['batch_size'] = 16 if row['length'] < 200 else 1
        row['sequence_sha256'] = hashlib.sha256(seqs[row['seq_id']].encode()).hexdigest()
    (out / 'selected.fa').write_text(''.join(f">{r['seq_id']}\n{seqs[r['seq_id']]}\n" for r in selected))
    binary = ROOT / 'vendor/LinearRaccess/build-raccess/LinRacc'
    shutil.copy2(binary, out / 'raw/LinRacc')
    from benchmark_common import verify_binary_build, load_json_object
    build_info = verify_binary_build(out / 'raw/LinRacc', load_json_object(ROOT / 'software.lock.json')['software']['commit'])
    m = dict(schema_version=1, experiment_id=out.name, created_utc=datetime.now(timezone.utc).isoformat(),
             purpose='Paired speed and full-output difference of fast/standard logsumexp at beam=100',
             dataset=dict(path=str(dataset.relative_to(ROOT)), sha256=sha256(dataset)),
             selected_fasta_sha256=sha256(out / 'selected.fa'), selection=selected,
             sampling=dict(seed=20260914, per_stratum=4, population=len(records),
                           rule='Uniform sample without replacement from sorted IDs in each length stratum; sentinels excluded from primary aggregates'),
             executable=dict(original_path=str(binary.relative_to(ROOT)), sha256=sha256(binary),
                             observed_source_revision=subprocess.check_output(['git', '-C', str(binary.parents[1]), 'rev-parse', 'HEAD'], text=True).strip(),
                             observed_source_status=subprocess.check_output(['git', '-C', str(binary.parents[1]), 'status', '--porcelain'], text=True),
                             build_source_revision=build_info['source_revision'], build_info=build_info),
             script_sha256=sha256(Path(__file__)), helper_sha256=sha256(ROOT / 'benchmark/scripts/audit_bprna.py'),
             job_script_sha256=sha256(ROOT / 'benchmark/scripts/uge/bprna_logsumexp_tradeoff.sh'),
             settings=dict(engine='lincapr', energy='turner2004', beam=100, c_multi=30, c_hairpin=30, access_lens=list(LENS)),
             measurement=dict(warmups_per_mode=1, measured_repetitions=4,
                              order='Initial mode chosen by sequence SHA256 parity; reverse each measured repetition',
                              short_sequence_batch_size=16, short_threshold_nt=200,
                              timer='time.perf_counter around GNU time subprocess; includes process startup and output IO',
                              memory='GNU time maximum resident set size (KiB) for a process, not per-record memory',
                              workers=4, affinity='One distinct physical CPU per worker; shared node memory and filesystem',
                              timeout_per_process_sec=180, threads_per_process=1),
             interpretation=dict(comparator='Same beam/caps with -no-fast-logsumexp; not an exact folding reference',
                                 numerical_tolerance=1e-12, probabilities_not_clipped=True,
                                 relative_errors='Not used as a primary metric near zero probability',
                                 potential_confound='At beam=100 a logsumexp change can also change retained states'))
    dump(out / 'manifest.json', m)
    print(f'Prepared {len(selected)} sequences in {out}', flush=True)


def run_lane(out_str, cpu, rows, m):
    out = Path(out_str)
    os.sched_setaffinity(0, {cpu})
    seqs = dict(fasta(out / 'selected.fa'))
    completed = []
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    for row in rows:
        name, size = row['seq_id'], row['batch_size']
        folder = out / 'raw' / name
        folder.mkdir()
        names = [f'{name}__copy{i:02d}' for i in range(size)]
        inp = folder / 'input.fa'
        inp.write_text(''.join(f'>{n}\n{seqs[name]}\n' for n in names))
        first_mode = int(row['sequence_sha256'][-1], 16) % 2
        modes = ['fast', 'exact'] if first_mode == 0 else ['exact', 'fast']
        observations, probability_maps, hashes = [], {}, {}
        for repeat in range(-1, m['measurement']['measured_repetitions']):
            order = modes if repeat < 0 or repeat % 2 == 0 else modes[::-1]
            for position, mode in enumerate(order):
                label = f'{repeat}_{mode}'
                output = folder / (label + '.txt')
                timing = folder / (label + '.time')
                logfile = folder / (label + '.log')
                cmd = ['/usr/bin/timeout', str(m['measurement']['timeout_per_process_sec']),
                       '/usr/bin/time', '-f', '%U\t%S\t%M', '-o', str(timing),
                       str(out / 'raw/LinRacc'), f'-seqfile={inp}', f'-outfile={output}',
                       '-engine=lincapr', '-energy=turner2004', '-beam=100', '-c_multi=30',
                       '-c_hairpin=30', '-access_len=1,3,7,10', '-probabilities', '-metadata']
                if mode == 'exact':
                    cmd.append('-no-fast-logsumexp')
                start = time.perf_counter()
                with logfile.open('w') as log:
                    proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
                elapsed = time.perf_counter() - start
                if proc.returncode:
                    raise RuntimeError(f'{name} {label}: exit {proc.returncode}')
                vals = parse_output(output, names, row['length'], mode)
                if mode in probability_maps and probability_maps[mode] != vals:
                    raise ValueError(f'Nonreproducible output: {name} {mode}')
                probability_maps[mode] = vals
                hashes[label] = dict(output=sha256(output), log=sha256(logfile))
                user, system, rss = timing.read_text().strip().split('\t')
                observations.append(dict(repetition=repeat, warmup=repeat == -1, mode=mode,
                                         position=position, wall_sec=elapsed, per_sequence_wall_sec=elapsed / size,
                                         user_sec=float(user), system_sec=float(system), maxrss_kib=int(rss)))
        fast, exact = probability_maps['fast'], probability_maps['exact']
        stats_by_lens = {}
        for k in LENS:
            keys = [key for key in fast if key[1] == k]
            delta = [abs(fast[key] - exact[key]) for key in keys]
            stats_by_lens[k] = dict(abs_error=distribution(delta),
                                    fast=bounds(fast[key] for key in keys), exact=bounds(exact[key] for key in keys))
        measured = [o for o in observations if not o['warmup']]
        fast_times = [o['per_sequence_wall_sec'] for o in measured if o['mode'] == 'fast']
        exact_times = [o['per_sequence_wall_sec'] for o in measured if o['mode'] == 'exact']
        ratios = [next(o['wall_sec'] for o in measured if o['repetition'] == i and o['mode'] == 'exact') /
                  next(o['wall_sec'] for o in measured if o['repetition'] == i and o['mode'] == 'fast') for i in range(4)]
        result = dict(**row, cpu=cpu, observations=observations, stats_by_lens=stats_by_lens,
                      fast_median_sec=statistics.median(fast_times), exact_median_sec=statistics.median(exact_times),
                      speedup_exact_over_fast=statistics.median(exact_times) / statistics.median(fast_times),
                      paired_speedup=distribution(ratios), output_and_log_hashes=hashes)
        dump(out / f'{name}.json', result)
        completed.append(result)
        print(f"Completed {name} ({row['length']} nt), speedup {result['speedup_exact_over_fast']:.3f}x", flush=True)
    return completed


def summarize(out, results):
    groups = {}
    for label in sorted({r['stratum'] for r in results}) + ['all_random']:
        rows = [r for r in results if (r['cohort'] == 'random' if label == 'all_random' else r['stratum'] == label)]
        groups[label] = dict(n=len(rows), speedup=distribution([r['speedup_exact_over_fast'] for r in rows]),
                             max_abs_probability_difference=max(s['abs_error']['max'] for r in rows for s in r['stats_by_lens'].values()),
                             fast_outside_1e_12=sum(s['fast']['outside_1e_12'] for r in rows for s in r['stats_by_lens'].values()),
                             exact_outside_1e_12=sum(s['exact']['outside_1e_12'] for r in rows for s in r['stats_by_lens'].values()),
                             windows=sum(s['fast']['count'] for r in rows for s in r['stats_by_lens'].values()))
    with (out / 'per_sequence.tsv').open('w') as h:
        columns = ['seq_id', 'cohort', 'stratum', 'length', 'batch_size', 'fast_median_sec', 'exact_median_sec', 'speedup_exact_over_fast', 'max_abs_error', 'fast_outside_1e_12', 'exact_outside_1e_12']
        writer = csv.DictWriter(h, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for r in results:
            row = {k: r[k] for k in columns if k in r}
            row.update(max_abs_error=max(s['abs_error']['max'] for s in r['stats_by_lens'].values()),
                       fast_outside_1e_12=sum(s['fast']['outside_1e_12'] for s in r['stats_by_lens'].values()),
                       exact_outside_1e_12=sum(s['exact']['outside_1e_12'] for s in r['stats_by_lens'].values()))
            writer.writerow(row)
    dump(out / 'summary.json', dict(completed=True, sequences=len(results), groups=groups,
                                   result_hashes={p.name: sha256(p) for p in sorted(out.glob('bpRNA_*.json'))},
                                   limitations=['24 stratified random sequences, not a population-weighted estimate; 3 sentinels reported separately.',
                                                'Single node; four pinned workers share memory and filesystem.',
                                                'Elapsed time includes IO and process startup, amortized over 16 copies for sequences below 200 nt.',
                                                'Differences are relative to standard logsumexp at beam=100, not exact folding probabilities.']))


def run(out):
    m = load(out / 'manifest.json')
    for path, expected in [(out / 'raw/LinRacc', m['executable']['sha256']),
                           (out / 'selected.fa', m['selected_fasta_sha256']),
                           (Path(__file__), m['script_sha256']),
                           (ROOT / 'benchmark/scripts/audit_bprna.py', m['helper_sha256']),
                           (ROOT / 'benchmark/scripts/uge/bprna_logsumexp_tradeoff.sh', m['job_script_sha256'])]:
        if sha256(path) != expected:
            raise ValueError(f'Frozen input changed: {path}')
    if (out / 'run_info.json').exists():
        raise ValueError('Run already started; use a new bundle')
    cpus = physical_cpus()[:m['measurement']['workers']]
    if len(cpus) != m['measurement']['workers']:
        raise ValueError('Insufficient distinct physical CPUs')
    dump(out / 'run_info.json', dict(job_id=os.environ.get('JOB_ID'), hostname=platform.node(),
                                   started_utc=datetime.now(timezone.utc).isoformat(), cpus=cpus,
                                   python=platform.python_version(), uname=list(platform.uname()),
                                   cpuinfo=Path('/proc/cpuinfo').read_text().split('\n\n')[0]))
    lanes = [[] for _ in cpus]
    for i, row in enumerate(sorted(m['selection'], key=lambda r: -r['length'])):
        lanes[i % len(cpus)].append(row)
    results = []
    with ProcessPoolExecutor(max_workers=len(cpus)) as pool:
        tasks = [pool.submit(run_lane, str(out), cpu, rows, m) for cpu, rows in zip(cpus, lanes)]
        for task in as_completed(tasks):
            results.extend(task.result())
    summarize(out, sorted(results, key=lambda r: (r['cohort'], r['length'], r['seq_id'])))
    print('Experiment complete', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'run'])
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    {'prepare': prepare, 'run': run}[args.command](args.out.resolve())


if __name__ == '__main__':
    main()
