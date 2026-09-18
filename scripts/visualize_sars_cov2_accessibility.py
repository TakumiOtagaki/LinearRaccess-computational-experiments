#!/usr/bin/env python3
"""Visualize SARS-CoV-2 accessibility analysis results."""
import argparse
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np


def load_results(filepath):
    """Load analysis results from JSON."""
    with open(filepath) as f:
        return json.load(f)


def plot_accessibility_heatmap(data, output_path):
    """Create heatmap of accessibility across genome positions and window sizes."""
    top_regions = data['top_accessible_regions']
    seq_len = data['sequence_length']
    
    # Get window sizes
    window_sizes = sorted([int(w) for w in top_regions.keys()])
    
    # Create figure
    fig, axes = plt.subplots(len(window_sizes), 1, figsize=(16, 2*len(window_sizes)))
    if len(window_sizes) == 1:
        axes = [axes]
    
    for ax, w in zip(axes, window_sizes):
        regions = top_regions[str(w)]
        
        # Extract positions and values (filter out inf values)
        valid_regions = [r for r in regions if np.isfinite(r['accessibility'])]
        positions = [r['position'] for r in valid_regions]
        values = [r['accessibility'] for r in valid_regions]
        
        # Plot as scatter
        ax.scatter(positions, values, alpha=0.6, s=50, c='steelblue')
        
        ax.set_xlabel('Genome position (nt)', fontsize=10)
        ax.set_ylabel('Opening energy (-RT log P) [kcal/mol]', fontsize=9)
        ax.set_title(f'Window size {w} nt - Top accessible regions', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, seq_len)
        
        # Annotate top 3 (from valid regions)
        for i in range(min(3, len(valid_regions))):
            ax.annotate(f"{valid_regions[i]['position']:,}", 
                       xy=(positions[i], values[i]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)
    
    plt.suptitle(f"SARS-CoV-2 Genome Accessibility Analysis\n{data['sequence_id']}", 
                 fontsize=14, y=1.00)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {output_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close()


def plot_sliding_window_profile(data, output_path):
    """Plot sliding window accessibility profile."""
    sliding_data = data['sliding_window_analysis']
    windows = sliding_data['top_regions']
    seq_len = data['sequence_length']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8))
    
    # Plot 1: Top windows along genome
    starts = [w['start'] for w in windows]
    ends = [w['end'] for w in windows]
    means = [w['mean_accessibility'] for w in windows]
    maxs = [w['max_accessibility'] for w in windows]
    
    # Plot as horizontal bars
    for i, (start, end, mean_val, max_val) in enumerate(zip(starts, ends, means, maxs)):
        alpha = 0.7 if i < 20 else 0.3
        ax1.barh(i, end-start, left=start, height=0.8, 
                alpha=alpha, color='steelblue')
    
    ax1.set_xlabel('Genome position (nt)', fontsize=11)
    ax1.set_ylabel('Window rank', fontsize=11)
    ax1.set_title('High-Accessibility Sliding Windows', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='x')
    ax1.set_xlim(0, seq_len)
    
    # Plot 2: Mean accessibility values
    ax2.plot(range(len(means)), means, 'o-', color='steelblue', markersize=6)
    ax2.set_xlabel('Window rank', fontsize=11)
    ax2.set_ylabel('Mean opening energy [kcal/mol]', fontsize=11)
    ax2.set_title('Mean Opening Energy by Window Rank', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f"SARS-CoV-2 Sliding Window Analysis ({sliding_data['window_size']} nt windows)", 
                 fontsize=13)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {output_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close()


def plot_top_regions_comparison(data, output_path):
    """Compare top accessible regions across different window sizes."""
    top_regions = data['top_accessible_regions']
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    window_sizes = sorted([int(w) for w in top_regions.keys()])
    colors = plt.cm.viridis(np.linspace(0, 1, len(window_sizes)))
    
    for w, color in zip(window_sizes, colors):
        regions = top_regions[str(w)]
        positions = [r['position'] for r in regions]
        values = [r['accessibility'] for r in regions]
        
        ax.plot(range(1, len(values)+1), values, 'o-', 
               label=f'{w} nt', color=color, markersize=6, linewidth=2)
    
    ax.set_xlabel('Rank', fontsize=12)
    ax.set_ylabel('Opening energy (-RT log P) [kcal/mol]', fontsize=12)
    ax.set_title('Top Accessible Regions by Window Size', fontsize=13)
    ax.legend(title='Window size', loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {output_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize SARS-CoV-2 accessibility analysis'
    )
    parser.add_argument('--input', type=Path, required=True,
                        help='Input JSON file from analyze_sars_cov2_accessibility.py')
    parser.add_argument('--output-dir', type=Path,
                        default=Path('manuscript/figures'),
                        help='Output directory for figures')
    
    args = parser.parse_args()
    
    # Load data
    print(f"Loading results from {args.input}...")
    data = load_results(args.input)
    
    print(f"Sequence: {data['sequence_id']}")
    print(f"Length: {data['sequence_length']:,} nt")
    print(f"Window sizes: {list(data['top_accessible_regions'].keys())}")
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    print("\nGenerating visualizations...")
    
    # Plot 1: Heatmap
    print("1. Accessibility heatmap...")
    plot_accessibility_heatmap(
        data, 
        args.output_dir / 'sars_cov2_accessibility_heatmap.png'
    )
    
    # Plot 2: Sliding window profile
    print("2. Sliding window profile...")
    plot_sliding_window_profile(
        data,
        args.output_dir / 'sars_cov2_sliding_windows.png'
    )
    
    # Plot 3: Top regions comparison
    print("3. Top regions comparison...")
    plot_top_regions_comparison(
        data,
        args.output_dir / 'sars_cov2_top_regions.png'
    )
    
    print("\n✓ All visualizations complete!")
    
    # Print summary for paper
    print("\n" + "="*70)
    print("PAPER SUMMARY")
    print("="*70)
    print(f"\nDataset: {data['sequence_id']}")
    print(f"Genome length: {data['sequence_length']:,} nt")
    print(f"Beam width: {data['beam_width']}")
    print(f"\nWindow sizes analyzed: {', '.join(map(str, sorted([int(w) for w in data['top_accessible_regions'].keys()])))}")
    
    print("\nTop 5 most accessible regions by window size:")
    for w in sorted([int(w) for w in data['top_accessible_regions'].keys()]):
        regions = data['top_accessible_regions'][str(w)]
        if regions:
            top = regions[0]
            print(f"  {w:3d} nt: Position {top['position']:>6,} (energy: {top['accessibility']:.3f} kcal/mol)")
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
