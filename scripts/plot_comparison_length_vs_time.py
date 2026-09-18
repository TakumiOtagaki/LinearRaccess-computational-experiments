#!/usr/bin/env python3
"""
Plot Figure 2: bpRNA Length vs Runtime Comparison
X-axis: Sequence length (nt)
Y-axis: Runtime (seconds)
Compares: LinearRaccess, Raccess Exact, Raccess Window (W=100)
"""

import json
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["font.size"] = 14


def load_linearraccess_data(results_dir: Path) -> tuple[list[int], list[float], dict]:
    """Load LinearRaccess data from batch run_summary.json files."""
    lengths = []
    times = []
    params = {}
    
    # Try to load parameters from manifest.json
    manifest_file = results_dir / "manifest.json"
    if manifest_file.exists():
        with open(manifest_file) as f:
            manifest = json.load(f)
        if "settings" in manifest:
            params = manifest["settings"]
    
    batch_dirs = sorted(results_dir.glob("batch_*"))
    
    for batch_dir in batch_dirs:
        summary_file = batch_dir / "run_summary.json"
        if not summary_file.exists():
            continue
            
        with open(summary_file) as f:
            data = json.load(f)
        
        # Also try to get params from batch manifest if not already loaded
        if not params:
            batch_manifest = batch_dir / "manifest.json"
            if batch_manifest.exists():
                with open(batch_manifest) as f:
                    manifest = json.load(f)
                if "settings" in manifest:
                    params = manifest["settings"]
        
        for record in data.get("records", []):
            if record.get("status") == "SUCCESS":
                lengths.append(record["length"])
                times.append(record["elapsed_sec"])
    
    return lengths, times, params


def load_raccess_data(results_dir: Path, method_name: str) -> tuple[list[int], list[float]]:
    """Load Raccess data from individual JSON files."""
    lengths = []
    times = []
    
    json_files = sorted(results_dir.glob("*.json"))
    
    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)
        
        if data.get("status") == "SUCCESS":
            lengths.append(data["seq_length"])
            times.append(data["elapsed_sec"])
    
    return lengths, times


def plot_comparison(lengths_dict: dict, times_dict: dict, output_path: Path, linearraccess_params: dict = None):
    """Create comparison scatter plot of length vs time."""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    colors = {
        'LinearRaccess': '#d62728',            # Red
        'Raccess Exact': '#2ca02c',             # Green
        'Raccess Window (W=100)': '#1f77b4'     # Blue
    }
    
    markers = {
        'LinearRaccess': 'o',
        'Raccess Exact': 's',
        'Raccess Window (W=100)': '^'
    }
    
    for method in ['LinearRaccess', 'Raccess Exact', 'Raccess Window (W=100)']:
        if method not in lengths_dict:
            continue
            
        lengths_arr = np.array(lengths_dict[method])
        times_arr = np.array(times_dict[method])
        
        # Build label with parameters
        if method == 'LinearRaccess' and linearraccess_params:
            beam = linearraccess_params.get('beam', 'N/A')
            label = f'LinearRaccess (b={beam}) (n={len(lengths_arr):,})'
        elif method == 'Raccess Window (W=100)':
            label = f'Raccess Window (W=100) (n={len(lengths_arr):,})'
        else:
            label = f'{method} (n={len(lengths_arr):,})'
        
        # Scatter plot
        ax.scatter(lengths_arr, times_arr, alpha=0.4, s=15, 
                   c=colors[method], marker=markers[method],
                   edgecolors='none', label=label)
        
        # Linear fit (only for linear methods)
        if method != 'Raccess Exact' and len(lengths_arr) > 1:
            slope, intercept, r_value, p_value, std_err = stats.linregress(lengths_arr, times_arr)
            x_fit = np.linspace(lengths_arr.min(), lengths_arr.max(), 100)
            y_fit = slope * x_fit + intercept
            ax.plot(x_fit, y_fit, color=colors[method], linewidth=2, alpha=0.8,
                   label=f'{method} fit (R²={r_value**2:.3f})')
    
    # Linear scale for Y axis
    ax.set_xlabel('Sequence length (nt)', fontsize=18)
    ax.set_ylabel('Runtime (seconds)', fontsize=18)
    ax.set_title('Runtime Comparison: LinearRaccess vs Raccess', fontsize=18, fontweight='bold')
    
    # Grid
    ax.grid(True, alpha=0.3, which='both')
    
    # Legend
    ax.legend(loc='upper left', fontsize=13)
    ax.tick_params(axis='both', labelsize=14)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    
    print(f"Figure saved to: {output_path}")
    print(f"Figure saved to: {output_path.with_suffix('.pdf')}")
    
    # Print statistics
    print(f"\nStatistics:")
    for method in ['LinearRaccess', 'Raccess Exact', 'Raccess Window (W=100)']:
        if method not in lengths_dict:
            continue
        lengths = lengths_dict[method]
        times = times_dict[method]
        print(f"\n  {method}:")
        print(f"    Sequences: {len(lengths):,}")
        print(f"    Length range: {min(lengths)} - {max(lengths)} nt")
        print(f"    Time range: {min(times):.3f} - {max(times):.3f} sec")


def main():
    parser = argparse.ArgumentParser(description='Plot bpRNA Length vs Runtime Comparison')
    parser.add_argument('--linearraccess-dir', type=Path, 
                        default=Path('benchmark/results/bprna_runs/8666393'),
                        help='Directory containing LinearRaccess batch results')
    parser.add_argument('--raccess-exact-dir', type=Path,
                        default=Path('benchmark/results/raccess_exact'),
                        help='Directory containing Raccess exact results')
    parser.add_argument('--raccess-w100-dir', type=Path,
                        default=Path('benchmark/results/raccess_w100'),
                        help='Directory containing Raccess window (W=100) results')
    parser.add_argument('--output', type=Path, 
                        default=Path('manuscript/figures/figure2_length_vs_time_comparison.png'),
                        help='Output figure path')
    
    args = parser.parse_args()
    
    lengths_dict = {}
    times_dict = {}
    linearraccess_params = {}
    
    # Load LinearRaccess data
    if args.linearraccess_dir.exists():
        print(f"Loading LinearRaccess data from: {args.linearraccess_dir}")
        lengths, times, linearraccess_params = load_linearraccess_data(args.linearraccess_dir)
        if lengths:
            lengths_dict['LinearRaccess'] = lengths
            times_dict['LinearRaccess'] = times
            print(f"  Loaded {len(lengths):,} sequences")
            if linearraccess_params:
                print(f"  Parameters: {linearraccess_params}")
    
    # Load Raccess Exact data
    if args.raccess_exact_dir.exists():
        print(f"Loading Raccess Exact data from: {args.raccess_exact_dir}")
        lengths, times = load_raccess_data(args.raccess_exact_dir, 'raccess-exact')
        if lengths:
            lengths_dict['Raccess Exact'] = lengths
            times_dict['Raccess Exact'] = times
            print(f"  Loaded {len(lengths):,} sequences")
    
    # Load Raccess Window data
    if args.raccess_w100_dir.exists():
        print(f"Loading Raccess Window (W=100) data from: {args.raccess_w100_dir}")
        lengths, times = load_raccess_data(args.raccess_w100_dir, 'raccess-w100')
        if lengths:
            lengths_dict['Raccess Window (W=100)'] = lengths
            times_dict['Raccess Window (W=100)'] = times
            print(f"  Loaded {len(lengths):,} sequences")
    
    if not lengths_dict:
        print("Error: No data loaded")
        return
    
    plot_comparison(lengths_dict, times_dict, args.output, linearraccess_params)


if __name__ == '__main__':
    main()
