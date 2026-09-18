#!/usr/bin/env python3
"""Plot Figure 3: Sequence Length vs Memory Usage (RSS) comparison.

Compares three methods:
- LinearRaccess (beam=100)
- Raccess Exact
- Raccess Window (W=100)
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from scipy import stats

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams["font.size"] = 14


def load_memory_data(data_dir: Path, method_name: str):
    """Load memory usage data from JSON files."""
    lengths = []
    memory = []
    
    if not data_dir.exists():
        print(f"Warning: Directory not found: {data_dir}")
        return lengths, memory
    
    # Check if this is LinearRaccess batch data (run_summary.json format)
    for summary_file in data_dir.rglob("run_summary.json"):
        try:
            with open(summary_file, 'r') as f:
                data = json.load(f)
            
            for record in data.get("records", []):
                if record.get("status") == "SUCCESS" and record.get("peak_rss_mb") is not None:
                    lengths.append(record["length"])
                    memory.append(record["peak_rss_mb"])
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            continue
    
    # Also check for individual JSON files (Raccess format)
    if len(lengths) == 0:
        for json_file in data_dir.glob("*.json"):
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                
                if data.get("status") == "SUCCESS" and data.get("memory_kb") is not None:
                    lengths.append(data["seq_length"])
                    # Convert KB to MB
                    memory.append(data["memory_kb"] / 1024)
            except (json.JSONDecodeError, KeyError, TypeError) as e:
                continue
    
    return lengths, memory


def plot_memory_comparison(lengths_dict, memory_dict, output_path, output_pdf):
    """Create scatter plot comparing memory usage of three methods."""
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    colors = {
        'LinearRaccess (b=100)': '#d62728',  # Red
        'Raccess Exact': '#2ca02c',           # Green
        'Raccess Window (W=100)': '#1f77b4'   # Blue
    }
    
    markers = {
        'LinearRaccess (b=100)': 'o',
        'Raccess Exact': 's',
        'Raccess Window (W=100)': '^'
    }
    
    for method_name in ['LinearRaccess (b=100)', 'Raccess Exact', 'Raccess Window (W=100)']:
        if method_name not in lengths_dict or len(lengths_dict[method_name]) == 0:
            print(f"Warning: No data for {method_name}")
            continue
        
        lengths = np.array(lengths_dict[method_name])
        memory = np.array(memory_dict[method_name])
        
        # Plot scatter
        ax.scatter(lengths, memory, alpha=0.5, s=20, 
                  color=colors[method_name], marker=markers[method_name],
                  label=method_name)
        
        # Fit line (only if we have enough points)
        if len(lengths) > 1:
            # Use linear fit for LinearRaccess and Raccess Window
            # Use quadratic fit for Raccess Exact memory (O(N^2); time is
            # cubic but memory is quadratic for full-span exact folding)
            if 'Exact' in method_name:
                coeffs = np.polyfit(lengths, memory, 2)
                fit_line = np.poly1d(coeffs)
                x_fit = np.linspace(min(lengths), max(lengths), 100)
                ax.plot(x_fit, fit_line(x_fit), '--',
                       color=colors[method_name], alpha=0.7, linewidth=2)
            else:
                # Linear fit
                slope, intercept, r_value, p_value, std_err = stats.linregress(lengths, memory)
                x_fit = np.linspace(min(lengths), max(lengths), 100)
                y_fit = slope * x_fit + intercept
                ax.plot(x_fit, y_fit, '--', 
                       color=colors[method_name], alpha=0.7, linewidth=2)
    ax.set_xlabel('Sequence length (nt)', fontsize=18)
    ax.set_ylabel('Peak RSS memory (MB)', fontsize=18)
    ax.set_title('Memory Usage Comparison: LinearRaccess vs Raccess', fontsize=18, fontweight='bold')
    ax.legend(loc='upper left', fontsize=13)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    
    # Save
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
    print(f"\nFigure saved to: {output_path}")
    print(f"Figure saved to: {output_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot Figure 3: Memory Usage Comparison'
    )
    parser.add_argument('--linearraccess-dir', type=Path,
                       default=Path('benchmark/results/bprna_memory_runs/8677021'),
                       help='Directory containing LinearRaccess memory results')
    parser.add_argument('--raccess-exact-dir', type=Path,
                       default=Path('benchmark/results/raccess_exact'),
                       help='Directory containing Raccess exact results')
    parser.add_argument('--raccess-w100-dir', type=Path,
                       default=Path('benchmark/results/raccess_w100'),
                       help='Directory containing Raccess window (W=100) results')
    parser.add_argument('--output', type=Path,
                       default=Path('manuscript/figures/figure3_memory.png'),
                       help='Output PNG file path')
    parser.add_argument('--output-pdf', type=Path,
                       default=Path('manuscript/figures/figure3_memory.pdf'),
                       help='Output PDF file path')
    
    args = parser.parse_args()
    
    # Load data from all three methods
    lengths_dict = {}
    memory_dict = {}
    
    print("Loading memory data...")
    
    # LinearRaccess
    print(f"\n1. LinearRaccess (b=100):")
    lengths, memory = load_memory_data(args.linearraccess_dir, 'linearraccess')
    if lengths:
        lengths_dict['LinearRaccess (b=100)'] = lengths
        memory_dict['LinearRaccess (b=100)'] = memory
        print(f"   Loaded {len(lengths):,} sequences")
        print(f"   Length range: {min(lengths)} - {max(lengths)} nt")
        print(f"   Memory range: {min(memory):.2f} - {max(memory):.2f} MB")
    
    # Raccess Exact
    print(f"\n2. Raccess Exact:")
    lengths, memory = load_memory_data(args.raccess_exact_dir, 'raccess-exact')
    if lengths:
        lengths_dict['Raccess Exact'] = lengths
        memory_dict['Raccess Exact'] = memory
        print(f"   Loaded {len(lengths):,} sequences")
        print(f"   Length range: {min(lengths)} - {max(lengths)} nt")
        print(f"   Memory range: {min(memory):.2f} - {max(memory):.2f} MB")
    
    # Raccess Window (W=100)
    print(f"\n3. Raccess Window (W=100):")
    lengths, memory = load_memory_data(args.raccess_w100_dir, 'raccess-w100')
    if lengths:
        lengths_dict['Raccess Window (W=100)'] = lengths
        memory_dict['Raccess Window (W=100)'] = memory
        print(f"   Loaded {len(lengths):,} sequences")
        print(f"   Length range: {min(lengths)} - {max(lengths)} nt")
        print(f"   Memory range: {min(memory):.2f} - {max(memory):.2f} MB")
    
    # Plot
    print("\nGenerating plot...")
    plot_memory_comparison(lengths_dict, memory_dict, args.output, args.output_pdf)


if __name__ == '__main__':
    main()
