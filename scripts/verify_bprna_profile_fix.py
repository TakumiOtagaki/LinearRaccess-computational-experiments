#!/usr/bin/env python3
"""Compare profile warnings and accessibility against frozen pre-fix outputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

from audit_bprna import fasta, sha256
from diagnose_bprna_logsumexp import read_values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    baseline, binary, out = args.baseline.resolve(), args.binary.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / 'raw').mkdir()
    manifest = json.loads((baseline / 'manifest.json').read_text())
    seqs = dict(fasta(baseline / 'selected.fa'))
    build_info = json.loads(subprocess.check_output([str(binary), '--build-info'], text=True))

    def compare(row):
        name = row['seq_id']
        inp = out / 'raw' / (name + '.fa')
        inp.write_text(f'>{name}\n{seqs[name]}\n')
        results = []
        for mode in ('fast', 'exact'):
            output = out / 'raw' / f'{name}_{mode}.txt'
            log = out / 'raw' / f'{name}_{mode}.log'
            command = [str(binary), f'-seqfile={inp}', f'-outfile={output}',
                       '-engine=lincapr', '-energy=turner2004', '-beam=100',
                       '-c_multi=30', '-c_hairpin=30', '-access_len=1,3,7,10',
                       '-probabilities', '-metadata']
            if mode == 'exact':
                command.append('-no-fast-logsumexp')
            with log.open('w') as handle:
                subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True, timeout=180)
            before_file = baseline / 'raw' / name / f'0_{mode}.txt'
            old, _ = read_values(before_file)
            new, metadata = read_values(output)
            if old.keys() != new.keys() or metadata['source_revision'] != build_info['source_revision']:
                raise ValueError('Coverage or embedded revision mismatch')
            before_log = baseline / 'raw' / name / f'0_{mode}.log'
            before_count = sum(l.startswith('warn: prob_sum[') for l in before_log.read_text().splitlines())
            if before_count % row['batch_size']:
                raise ValueError('Baseline warning count not divisible by batch size')
            after_count = sum(l.startswith('warn: prob_sum[') for l in log.read_text().splitlines())
            results.append(dict(seq_id=name, length=row['length'], mode=mode,
                                warnings_before=before_count // row['batch_size'], warnings_after=after_count,
                                max_abs_accessibility_change=max(abs(new[k] - old[k]) for k in old),
                                baseline_output_sha256=sha256(before_file), output_sha256=sha256(output),
                                log_sha256=sha256(log)))
        print(f'Compared {name}', flush=True)
        return results

    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = [r for group in pool.map(compare, manifest['selection']) for r in group]
    summary = dict(schema_version=1, build_info=build_info, binary_sha256=sha256(binary),
                   script_sha256=sha256(Path(__file__)), baseline_manifest_sha256=sha256(baseline / 'manifest.json'),
                   comparisons=len(rows), warnings_before=sum(r['warnings_before'] for r in rows),
                   warnings_after=sum(r['warnings_after'] for r in rows),
                   max_abs_accessibility_change=max(r['max_abs_accessibility_change'] for r in rows),
                   rows=rows,
                   limitations='27 selected sequences in two modes; not a full bpRNA rerun or timing benchmark.')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k:summary[k] for k in ('comparisons','warnings_before','warnings_after','max_abs_accessibility_change')}))
    return 0 if summary['max_abs_accessibility_change'] == 0 and summary['warnings_after'] == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
