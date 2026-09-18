#!/usr/bin/env python3
"""
Compute bpRNA statistics by length stratum for Table 2
"""

import json
import argparse
from pathlib import Path
import numpy as np


def load_run_data(results_dir: Path) -> list[dict]:
    """Load all sequence records from batch run_summary.json files."""
    records = []
    
    batch_dirs = sorted(results_dir.glob("batch_*"))
    
    for batch_dir in batch_dirs:
        summary_file = batch_dir / "run_summary.json"
        if not summary_file.exists():
            continue
            
        with open(summary_file) as f:
            data = json.load(f)
        
        for record in data.get("records", []):
            if record.get("status") == "SUCCESS":
                records.append({
                    'seq_id': record['seq_id'],
                    'length': record['length'],
                    'time': record['elapsed_sec']
                })
    
    return records


def compute_strata_stats(records: list[dict]) -> dict:
    """Compute statistics by length stratum."""
    
    # Define strata
    strata = [
        (50, 99, "50-99"),
        (100, 199, "100-199"),
        (200, 399, "200-399"),
        (400, 9999, "400+")
    ]
    
    results = {}
    
    for min_len, max_len, label in strata:
        subset = [r for r in records if min_len <= r['length'] <= max_len]
        
        if not subset:
            continue
        
        times = [r['time'] for r in subset]
        lengths = [r['length'] for r in subset]
        
        results[label] = {
            'count': len(subset),
            'length_median': np.median(lengths),
            'time_median': np.median(times),
            'time_p95': np.percentile(times, 95),
            'time_max': np.max(times),
        }
    
    return results


def print_table(records: list[dict]):
    """Print statistics table."""
    
    # Overall stats
    all_times = [r['time'] for r in records]
    all_lengths = [r['length'] for r in records]
    
    print("=" * 70)
    print("bpRNA Dataset Statistics")
    print("=" * 70)
    print(f"\nTable 1: Dataset Overview")
    print(f"  Total sequences: {len(records):,}")
    print(f"  Length range: {min(all_lengths)} - {max(all_lengths)} nt")
    print(f"  Median length: {np.median(all_lengths):.0f} nt")
    print(f"  Mean length: {np.mean(all_lengths):.2f} nt")
    
    # Count families
    families = {}
    for r in records:
        family = r['seq_id'].split('_')[1]
        families[family] = families.get(family, 0) + 1
    print(f"  Families: {len(families)} ({', '.join(sorted(families.keys()))})")
    
    # Strata stats
    print(f"\nTable 2: Runtime Statistics by Length Stratum")
    print("-" * 70)
    print(f"{'Length (nt)':<12} {'Count':>8} {'Median (s)':>12} {'P95 (s)':>12} {'Max (s)':>12}")
    print("-" * 70)
    
    strata_results = compute_strata_stats(records)
    
    for label in ["50-99", "100-199", "200-399", "400+"]:
        if label in strata_results:
            s = strata_results[label]
            print(f"{label:<12} {s['count']:>8,} {s['time_median']:>12.3f} {s['time_p95']:>12.3f} {s['time_max']:>12.3f}")
    
    print("-" * 70)
    
    # Additional stats
    print(f"\nAdditional Statistics:")
    print(f"  Overall median time: {np.median(all_times):.3f} sec")
    print(f"  Overall P95 time: {np.percentile(all_times, 95):.3f} sec")
    print(f"  Overall max time: {np.max(all_times):.3f} sec")
    print(f"  Sequences < 1 sec: {sum(1 for t in all_times if t < 1):,} ({100*sum(1 for t in all_times if t < 1)/len(all_times):.1f}%)")
    print(f"  Sequences < 10 sec: {sum(1 for t in all_times if t < 10):,} ({100*sum(1 for t in all_times if t < 10)/len(all_times):.1f}%)")
    
    print(f"\nFamily breakdown:")
    for family in sorted(families.keys()):
        count = families[family]
        print(f"  {family}: {count:,} ({100*count/len(records):.1f}%)")


def main():
    parser = argparse.ArgumentParser(description='Compute bpRNA statistics')
    parser.add_argument('--results-dir', type=Path, required=True,
                        help='Directory containing batch results')
    
    args = parser.parse_args()
    
    print(f"Loading data from: {args.results_dir}")
    records = load_run_data(args.results_dir)
    print(f"Loaded {len(records):,} sequences\n")
    
    print_table(records)


if __name__ == '__main__':
    main()
