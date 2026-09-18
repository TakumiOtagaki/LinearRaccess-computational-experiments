#!/usr/bin/env python3
"""Plot paired speed and probability differences from a completed experiment."""
import argparse
import json
from pathlib import Path
import statistics

from audit_bprna import distribution
from run_bprna_logsumexp_tradeoff import parse_output, bounds, dump

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    args = parser.parse_args()
    summary = json.loads((args.bundle / 'summary.json').read_text())
    if not summary['completed']:
        raise ValueError('Incomplete experiment')
    rows = [json.loads((args.bundle / name).read_text()) for name in summary['result_hashes']]
    pooled = {cohort: {k: {'delta': [], 'fast': [], 'exact': []} for k in (1, 3, 7, 10)}
              for cohort in ('random', 'sentinel')}
    for row in rows:
        names = [f"{row['seq_id']}__copy{i:02d}" for i in range(row['batch_size'])]
        folder = args.bundle / 'raw' / row['seq_id']
        fast = parse_output(folder / '0_fast.txt', names, row['length'], 'fast')
        exact = parse_output(folder / '0_exact.txt', names, row['length'], 'exact')
        for key in fast:
            bucket = pooled[row['cohort']][key[1]]
            bucket['delta'].append(abs(fast[key] - exact[key]))
            bucket['fast'].append(fast[key])
            bucket['exact'].append(exact[key])
    metrics = {}
    for cohort, widths in pooled.items():
        widths['all'] = {key: [v for bucket in widths.values() for v in bucket[key]]
                         for key in ('delta', 'fast', 'exact')}
        metrics[cohort] = {}
        for width, bucket in widths.items():
            metrics[cohort][str(width)] = dict(abs_difference=distribution(bucket['delta']),
                                         fast=bounds(bucket['fast']), exact=bounds(bucket['exact']),
                                         abs_difference_over_1e_4=sum(v > 1e-4 for v in bucket['delta']),
                                         abs_difference_over_1e_3=sum(v > 1e-3 for v in bucket['delta']))
    memory = {}
    for mode in ('fast', 'exact'):
        peaks = [statistics.median(o['maxrss_kib'] / 1024 for o in r['observations']
                                   if not o['warmup'] and o['mode'] == mode)
                 for r in rows if r['cohort'] == 'random']
        memory[mode] = distribution(peaks)
    dump(args.bundle / 'pooled_metrics.json', dict(
        metrics=metrics, process_peak_rss_mib=memory,
        interpretation='Window-weighted within the selected cohort; repeated copies counted once; not a population-weighted estimate.'))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for cohort, label, color, marker in [('random', 'Stratified random (n=24)', '#0072B2', 'o'),
                                          ('sentinel', 'Diagnostic examples (n=3)', '#D55E00', 'D')]:
        selected = [r for r in rows if r['cohort'] == cohort]
        lengths = [r['length'] for r in selected]
        speedups = [r['speedup_exact_over_fast'] for r in selected]
        lo = [r['speedup_exact_over_fast'] - r['paired_speedup']['min'] for r in selected]
        hi = [r['paired_speedup']['max'] - r['speedup_exact_over_fast'] for r in selected]
        axes[0].errorbar(lengths, speedups, yerr=[lo, hi], fmt=marker, color=color,
                         markersize=4, alpha=.85, capsize=2, linewidth=.7, label=label)
        errors = [max(s['abs_error']['max'] for s in r['stats_by_lens'].values()) for r in selected]
        axes[1].scatter(lengths, errors, marker=marker, color=color, s=23, alpha=.85, label=label)
    axes[0].axhline(1, color='0.4', linestyle='--', linewidth=1)
    axes[0].set_ylabel('Speedup (standard / fast elapsed time)')
    axes[0].set_title('A  Speed of the fast approximation', loc='left', fontsize=10)
    axes[1].set_ylabel('Maximum absolute probability difference')
    axes[1].set_yscale('log')
    axes[1].set_title('B  Difference from standard log-sum-exp', loc='left', fontsize=10)
    for ax in axes:
        ax.set_xscale('log')
        ax.set_xlabel('Sequence length (nt)')
        ax.grid(True, which='major', alpha=.2)
    axes[0].legend(fontsize=8, loc='best')
    fig.savefig(args.bundle / 'tradeoff.png', dpi=180)
    fig.savefig(args.bundle / 'tradeoff.pdf')
    plt.close(fig)
    print('Wrote tradeoff.png and tradeoff.pdf; bars show min/max of four paired ratios, not confidence intervals.')


if __name__ == '__main__':
    main()
