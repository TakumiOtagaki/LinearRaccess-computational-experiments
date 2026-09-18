#!/usr/bin/env python3
"""Plot Figure 3: Sequence Length vs Memory Usage (RSS) for LinearRaccess."""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import json


def load_memory_data(results_dir: Path) -> tuple[list[int], list[float]]:
    """Load memory usage data from batch run_summary.json files."""
    lengths = []
    memory = []
    
    batch_dirs = sorted(results_dir.glob("batch_*"))
    
    for batch_dir in batch_dirs:
        summary_file = batch_dir / "run_summary.json"
        if not summary_file.exists():
            continue
            
        with open(summary_file) as f:
            data = json.load(f)
        
        for record in data.get("records", []):
            if record.get("status") == "SUCCESS" and record.get("peak_rss_mb") is not None:
                lengths.append(record["length"])
                memory.append(record["peak_rss_mb"])
    
    return lengths, memory


def plot_memory(lengths: list, memory: list, output_path: Path, params: dict = None):
    """Create scatter plot of length vs memory with linear fit."""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    lengths_arr = np.array(lengths)
    memory_arr = np.array(memory)
    
    # Scatter plot
    ax.scatter(lengths_arr, memory_arr, alpha=0.5, s=20, 
               c='steelblue', edgecolors='none',
               label=f'Data points (n={len(lengths):,})')
    
    # Linear fit
    if len(lengths_arr) > 1:
        slope, intercept, r_value, p_value, std_err = stats.linregress(lengths_arr, memory_arr)
        x_fit = np.linspace(lengths_arr.min(), lengths_arr.max(), 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, color='crimson', linewidth=2, alpha=0.8,
               label=f'Linear fit (R²={r_value**2:.3f})')
    
    # Labels
    ax.set_xlabel('Sequence length (nt)', fontsize=12)
    ax.set_ylabel('Peak RSS memory (MB)', fontsize=12)
    
    # Title with parameters
    if params:
        beam = params.get('beam', 'N/A')
        ax.set_title(f'LinearRaccess (b={beam}) Memory Usage on bpRNA Dataset', fontsize=14)
    else:
        ax.set_title('LinearRaccess Memory Usage on bpRNA Dataset', fontsize=14)
    
    # Grid
    ax.grid(True, alpha=0.3, which='both')
    
    # Legend
    ax.legend(loc='upper left', fontsize=10)
    
    plt.tight_layout()
    
    # Save as PNG and PDF
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    
    return len(lengths), (min(lengths), max(lengths)), (min(memory), max(memory)), slope if len(lengths) > 1 else None


def load_params(results_dir: Path) -> dict:
    """Load parameters from manifest.json."""
    manifest_file = results_dir / "manifest.json"
    if manifest_file.exists():
        with open(manifest_file) as f:
            manifest = json.load(f)
        return manifest.get("settings", {})
    return {}


def main():
    parser = argparse.ArgumentParser(description='Plot Figure 3: Sequence Length vs Memory Usage')
    parser.add_argument('--results-dir', type=Path, 
                        default=Path('benchmark/results/bprna_memory_runs/8676148'),
                        help='Directory containing memory measurement batch results')
    parser.add_argument('--output', type=Path, 
                        default=Path('manuscript/figures/figure3_memory.png'),
                        help='Output figure path')
    
    args = parser.parse_args()
    
    # Check if results directory exists
    if not args.results_dir.exists():
        print(f"ERROR: Results directory not found: {args.results_dir}")
        print("\nNote: Job 8676148 is still running. Results will be available when completed.")
        print("Check job status with: qstat | grep 8676148")
        return 1
    
    # Load parameters
    params = load_params(args.results_dir)
    
    # Load data
    print(f"Loading memory data from: {args.results_dir}")
    lengths, memory = load_memory_data(args.results_dir)
    
    if not lengths:
        print("ERROR: No data loaded. Job may not be complete yet.")
        print(f"\nCheck if any batch directories exist:")
        batch_dirs = list(args.results_dir.glob("batch_*"))
        print(f"  Found {len(batch_dirs)} batch directories")
        if batch_dirs:
            for bd in batch_dirs[:3]:
                print(f"    {bd.name}")
        return 1
    
    # Plot
    count, length_range, memory_range, slope = plot_memory(lengths, memory, args.output, params)
    
    print(f"\nFigure saved to: {args.output}")
    print(f"Figure saved to: {args.output.with_suffix('.pdf')}")
    
    print("\nStatistics:")
    print(f"  Sequences: {count:,}")
    print(f"  Length range: {length_range[0]} - {length_range[1]} nt")
    print(f"  Memory range: {memory_range[0]:.2f} - {memory_range[1]:.2f} MB")
    if slope:
        print(f"  Slope: {slope:.4f} MB/nt")
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
