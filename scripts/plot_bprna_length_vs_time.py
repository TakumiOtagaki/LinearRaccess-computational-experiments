#!/usr/bin/env python3
"""
Plot Figure 2: bpRNA Length vs Runtime
X-axis: Sequence length (nt)
Y-axis: Runtime (seconds, log scale)
Scatter plot with linear fit
"""

import json
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def load_run_data(results_dir: Path) -> tuple[list[int], list[float]]:
    """Load sequence lengths and runtimes from all batch run_summary.json files."""
    lengths = []
    times = []
    
    batch_dirs = sorted(results_dir.glob("batch_*"))
    
    for batch_dir in batch_dirs:
        summary_file = batch_dir / "run_summary.json"
        if not summary_file.exists():
            continue
            
        with open(summary_file) as f:
            data = json.load(f)
        
        for record in data.get("records", []):
            if record.get("status") == "SUCCESS":
                lengths.append(record["length"])
                times.append(record["elapsed_sec"])
    
    return lengths, times


def plot_length_vs_time(lengths: list[int], times: list[float], output_path: Path):
    """Create scatter plot of length vs time with linear fit."""
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Convert to numpy arrays
    lengths_arr = np.array(lengths)
    times_arr = np.array(times)
    
    # Scatter plot
    ax.scatter(lengths_arr, times_arr, alpha=0.3, s=10, c='steelblue', edgecolors='none')
    
    # Linear fit (as per requirements: "LinearRaccessはlinearだから")
    slope, intercept, r_value, p_value, std_err = stats.linregress(lengths_arr, times_arr)
    
    # Fit line
    x_fit = np.linspace(lengths_arr.min(), lengths_arr.max(), 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, 'r-', linewidth=2, label=f'Linear fit (R² = {r_value**2:.4f})')
    
    # Linear scale for Y axis (to show linearity clearly)
    # ax.set_yscale('log')  # Commented out for linear scale
    
    # Labels
    ax.set_xlabel('Sequence length (nt)', fontsize=12)
    ax.set_ylabel('Runtime (seconds)', fontsize=12)
    ax.set_title(f'LinearRaccess Runtime on bpRNA (n = {len(lengths):,})', fontsize=14)
    
    # Grid
    ax.grid(True, alpha=0.3, which='both')
    
    # Legend
    ax.legend(loc='upper left')
    
    # Add fit equation text (simplified, fewer significant figures)
    equation_text = f'y = {slope:.2e}x + {intercept:.2e}'
    ax.text(0.95, 0.05, equation_text, transform=ax.transAxes, 
            fontsize=10, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    
    print(f"Figure saved to: {output_path}")
    print(f"Figure saved to: {output_path.with_suffix('.pdf')}")
    print(f"\nStatistics:")
    print(f"  Total sequences: {len(lengths):,}")
    print(f"  Length range: {min(lengths)} - {max(lengths)} nt")
    print(f"  Time range: {min(times):.3f} - {max(times):.3f} sec")
    print(f"  Linear fit R²: {r_value**2:.4f}")
    print(f"  Slope: {slope:.4e} sec/nt")
    
    return slope, intercept, r_value**2


def main():
    parser = argparse.ArgumentParser(description='Plot bpRNA Length vs Runtime')
    parser.add_argument('--results-dir', type=Path, required=True,
                        help='Directory containing batch results (e.g., benchmark/results/bprna_runs/8666393)')
    parser.add_argument('--output', type=Path, default=Path('figures/figure2_length_vs_time.png'),
                        help='Output figure path')
    
    args = parser.parse_args()
    
    print(f"Loading data from: {args.results_dir}")
    lengths, times = load_run_data(args.results_dir)
    
    print(f"Loaded {len(lengths):,} sequences")
    
    plot_length_vs_time(lengths, times, args.output)


if __name__ == '__main__':
    main()
