#!/usr/bin/env python3
"""Run one bpRNA chunk with RSS (memory) measurement."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import psutil


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def parse_fasta(text, minimum_length):
    records, current, parts = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current is not None:
                records.append((current, ''.join(parts)))
            if not line[1:].split():
                raise ValueError('Empty sequence ID')
            current, parts = line[1:].split()[0], []
        elif current is None:
            raise ValueError('Sequence before FASTA header')
        else:
            parts.append(line.upper().replace('T', 'U'))
    if current is not None:
        records.append((current, ''.join(parts)))
    if not records:
        raise ValueError('Empty input chunk')
    paths, ids = set(), set()
    for name, seq in records:
        safe = re.sub(r'[^A-Za-z0-9_.-]', '_', name)
        if name in ids or safe in paths:
            raise ValueError('Duplicate sequence ID or output filename collision')
        if len(seq) < minimum_length or not re.fullmatch('[ACGU]+', seq):
            raise ValueError(f'Invalid or too-short sequence: {name}')
        ids.add(name)
        paths.add(safe)
    return records


def validate_output(path, name, length, settings, build_info):
    metadata, headers, position = {}, [], 0
    outside = 0
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
            pos, payload = line.split('\t')
            if int(pos) != position:
                raise ValueError('Missing, duplicate or reordered output position')
            pairs = [item.split(',') for item in payload.rstrip(';').split(';')]
            if [int(k) for k, _ in pairs] != [k for k in settings['access_lens'] if position + k <= length]:
                raise ValueError('Missing or unexpected accessibility window')
            for _, raw in pairs:
                probability = float(raw)
                if not math.isfinite(probability):
                    raise ValueError('Nonfinite output probability')
                outside += probability < -1e-12 or probability > 1 + 1e-12
            position += 1
    if position != length or headers != [name]:
        raise ValueError('Incomplete output or wrong sequence identity')
    expected = dict(engine=settings['engine'], effective_energy_model=settings['energy'],
                    beam=str(settings['beam']), c_multi=str(settings['c_multi']),
                    c_hairpin=str(settings['c_hairpin']),
                    access_lens=','.join(map(str, settings['access_lens'])),
                    logsumexp='fast' if settings['fast_logsumexp'] else 'exact',
                    output_mode='probability', normalize_profiles='true',
                    software_version=build_info['software_version'],
                    source_revision=build_info['source_revision'])
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f'Output metadata mismatch: {key}')
    return outside


def run_sequence(number, name, sequence, args, batch, binary, settings, info):
    safe = re.sub(r'[^A-Za-z0-9_.-]', '_', name)
    output = batch / f'{safe}.txt'
    partial = batch / f'{safe}.txt.partial'
    log = batch / 'logs' / f'{safe}.log'
    row = dict(chunk_id=args.task_id, seq_num=number, seq_id=name, length=len(sequence),
               status='FAILED', elapsed_sec=0.0, peak_rss_mb=None, exit_code=None, error='',
               outside_tolerance_1e_12=None, output_sha256=None)
    
    start = time.perf_counter()
    peak_rss = 0
    
    try:
        with tempfile.TemporaryDirectory(dir=batch / '.work') as temporary:
            inp = Path(temporary) / 'sequence.fa'
            inp.write_text(f'>{name}\n{sequence}\n')
            command = [str(binary), f'-seqfile={inp}', f'-outfile={partial}',
                       '-access_len=' + ','.join(map(str, settings['access_lens'])),
                       f"-engine={settings['engine']}", f"-energy={settings['energy']}",
                       f"-beam={settings['beam']}", f"-c_multi={settings['c_multi']}",
                       f"-c_hairpin={settings['c_hairpin']}", '-probabilities', '-metadata']
            if not settings['fast_logsumexp']:
                command.append('-no-fast-logsumexp')
            
            # Run with RSS monitoring
            with log.open('w') as handle:
                proc = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT)
                
                # Monitor memory usage
                try:
                    p = psutil.Process(proc.pid)
                    while proc.poll() is None:
                        try:
                            rss = p.memory_info().rss
                            peak_rss = max(peak_rss, rss)
                            time.sleep(0.01)  # 10ms sampling
                        except psutil.NoSuchProcess:
                            break
                except psutil.NoSuchProcess:
                    pass
                
                proc.wait(timeout=args.timeout)
                completed = proc
            
            row['exit_code'] = completed.returncode
            row['peak_rss_mb'] = round(peak_rss / (1024 * 1024), 2)  # Convert to MB
            
            if completed.returncode:
                raise ValueError(f'LinRacc exited with status {completed.returncode}')
            row['outside_tolerance_1e_12'] = validate_output(
                partial, name, len(sequence), settings, info)
            partial.replace(output)
            row['output_sha256'] = digest(output)
            row['status'] = 'SUCCESS'
    except subprocess.TimeoutExpired:
        row['status'] = 'TIMEOUT'
        row['error'] = f'Exceeded {args.timeout} seconds'
    except Exception as error:
        row['error'] = f'{type(error).__name__}: {error}'
    
    row['elapsed_sec'] = time.perf_counter() - start
    return row


def run(args):
    if args.parallel < 1 or args.timeout <= 0 or not math.isfinite(args.timeout) or args.task_id < 1:
        raise ValueError('Task ID, parallelism and timeout must be positive')
    settings = dict(beam=100, c_multi=30, c_hairpin=30, access_lens=[1, 3, 7, 10],
                    engine='lincapr', energy='turner2004', fast_logsumexp=not args.no_fast_logsumexp)
    source_bytes = args.chunk.read_bytes()
    records = parse_fasta(source_bytes.decode(), max(settings['access_lens']))
    batch = args.outdir / f'batch_{args.task_id}'
    batch.mkdir(parents=True, exist_ok=False)
    for directory in ['.inputs', '.work', 'logs']:
        (batch / directory).mkdir()
    frozen_input = batch / '.inputs/chunk.fa'
    frozen_input.write_bytes(source_bytes)
    binary = batch / '.inputs/LinRacc'
    initial_hash = digest(args.linracc)
    shutil.copy2(args.linracc, binary)
    if digest(binary) != initial_hash or digest(args.linracc) != initial_hash:
        raise ValueError('Executable changed while being copied')
    result = subprocess.run([str(binary), '--build-info'], capture_output=True, text=True,
                            check=True, timeout=30)
    info = json.loads(result.stdout)
    if not re.fullmatch('[0-9a-f]{40}', info['source_revision']):
        raise ValueError('Production jobs require a clean Git revision embedded at build time')
    write_json(batch / 'manifest.json', dict(
        schema_version=1, job_id=os.environ.get('JOB_ID'), task_id=args.task_id,
        started_utc=datetime.now(timezone.utc).isoformat(), hostname=platform.node(),
        input=dict(original_path=str(args.chunk), snapshot='.inputs/chunk.fa',
                   sha256=digest(frozen_input), sequences=len(records),
                   normalization='uppercase, T to U'),
        binary=dict(original_path=str(args.linracc), snapshot='.inputs/LinRacc',
                    sha256=initial_hash, build_info=info), settings=settings,
        runner_sha256=digest(Path(__file__)), python_version=platform.python_version(),
        concurrency=args.parallel, timeout_sec=args.timeout,
        completion_rule='Every input record must exit zero and pass output identity, coverage, finite-value and metadata checks',
        probability_policy='Record values outside [-1e-12, 1+1e-12]; do not fail or clip finite approximation errors'))
    rows = []
    with (batch / 'summary.tsv').open('w') as handle, ThreadPoolExecutor(max_workers=args.parallel) as pool:
        writer = csv.DictWriter(handle, fieldnames=['chunk_id', 'seq_num', 'seq_id', 'length', 'status', 'elapsed_sec', 'peak_rss_mb'],
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        handle.flush()
        futures = {pool.submit(run_sequence, i, name, seq, args, batch, binary, settings, info): (i, name, seq)
                   for i, (name, seq) in enumerate(records, 1)}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            writer.writerow({key: row[key] for key in writer.fieldnames})
            handle.flush()
    success = sum(row['status'] == 'SUCCESS' for row in rows)
    complete = len(rows) == len(records) and success == len(records)
    write_json(batch / 'run_summary.json', dict(
        schema_version=1, completed=complete, expected=len(records), recorded=len(rows),
        success=success, failed=len(rows) - success, finished_utc=datetime.now(timezone.utc).isoformat(),
        records=sorted(rows, key=lambda row: row['seq_num'])))
    print(f'Chunk {args.task_id}: {success}/{len(records)} successful; completed={complete}', flush=True)
    return 0 if complete else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chunk', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--linracc', type=Path, required=True)
    parser.add_argument('--task-id', type=int, required=True)
    parser.add_argument('--parallel', type=int, default=10)
    parser.add_argument('--timeout', type=float, default=3600)
    parser.add_argument('--no-fast-logsumexp', action='store_true')
    args = parser.parse_args()
    for key in ['chunk', 'outdir', 'linracc']:
        setattr(args, key, getattr(args, key).resolve())
    try:
        return run(args)
    except Exception as error:
        print(f'ERROR: {type(error).__name__}: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
