#!/usr/bin/env python3
"""Create improved heatmap for SARS-CoV-2 accessibility with proper handling of inf values."""
import argparse
from pathlib import Path
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.patches as mpatches


def create_custom_colormap():
    """Create colormap: red for low energy (accessible), white for >=100 or invalid."""
    # Red gradient: low energy (accessible) -> higher energy (less accessible)
    # Darkest red (lowest energy) to white (highest energy before cutoff)
    colors = ['#67000d', '#a50f15', '#cb181d', '#ef3b2c', '#fb6a4a', 
              '#fc9272', '#fcbba1', '#fee0d2', '#fff5f0', 'white']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('accessibility', colors, N=n_bins)
    cmap.set_bad('white')  # NaN/inf will be white
    return cmap


def plot_genome_heatmap(data, output_path, max_energy=100, selected_windows=[3, 10, 30]):
    """Plot genome-wide accessibility heatmap without binning."""
    # Use full accessibility data (all positions)
    full_data = data.get('full_accessibility_data', data['top_accessible_regions'])
    seq_len = data['sequence_length']
    # Use only selected window sizes
    window_sizes = [w for w in selected_windows if str(w) in full_data]
    
    if not window_sizes:
        print(f"Warning: None of the selected windows {selected_windows} found in data")
        return
    
    # Wide figure: 1 pixel per nucleotide, but limited to reasonable width
    # Use aspect='auto' to stretch if needed
    fig_width = min(seq_len / 100, 40)  # Max 40 inches wide
    fig, axes = plt.subplots(len(window_sizes), 1, figsize=(fig_width, 2*len(window_sizes)))
    if len(window_sizes) == 1:
        axes = [axes]
    
    cmap = create_custom_colormap()
    
    for idx, (ax, w) in enumerate(zip(axes, window_sizes)):
        regions = full_data[str(w)]
        
        # Create genome-wide array (initialize with NaN for white background)
        genome_values = np.full(seq_len, np.nan)
        
        # Fill in values for all positions with accessibility data
        print(f"  Window {w} nt: {len(regions)} positions")
        for r in regions:
            pos = r['position']
            energy = r['accessibility']
            if pos < seq_len:
                # Show values <= max_energy as red gradient
                # Values > max_energy or inf become NaN (white)
                if np.isfinite(energy) and energy <= max_energy:
                    genome_values[pos] = energy
                else:
                    genome_values[pos] = np.nan  # White for inf or > max_energy
        
        # Reshape for imshow: (1, seq_len) - no binning
        genome_values_reshaped = genome_values.reshape(1, -1)
        
        # Plot heatmap without binning
        im = ax.imshow(genome_values_reshaped, aspect='auto', cmap=cmap, 
                      vmin=0, vmax=max_energy,
                      extent=[0, seq_len, 0, 1])
        
        ax.set_ylabel(f'{w} nt', fontsize=11, rotation=0, ha='right', va='center')
        ax.set_yticks([])
        
        if idx == len(window_sizes) - 1:
            ax.set_xlabel('Genome position (nt)', fontsize=12)
        
        # Add colorbar for the last subplot
        if idx == 0:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.01)
            cbar.set_label('Opening energy (kcal/mol)', fontsize=10)
            cbar.ax.text(1.2, 0, 'inf', transform=cbar.ax.transAxes, 
                        fontsize=9, va='bottom')
    
    plt.suptitle(f"SARS-CoV-2 Genome Accessibility Heatmap\n{data['sequence_id']}", 
                 fontsize=14, y=0.995)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {output_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close()


def plot_top_regions_bars(data, output_path):
    """Plot top accessible regions as horizontal bars with energy values."""
    top_regions = data['top_accessible_regions']
    seq_len = data['sequence_length']
    window_sizes = sorted([int(w) for w in top_regions.keys()])
    
    fig, axes = plt.subplots(len(window_sizes), 1, figsize=(16, 2*len(window_sizes)))
    if len(window_sizes) == 1:
        axes = [axes]
    
    for ax, w in zip(axes, window_sizes):
        regions = top_regions[str(w)]
        
        # Filter finite values and sort by energy
        valid = [(r['position'], r['accessibility']) for r in regions 
                if np.isfinite(r['accessibility'])]
        valid.sort(key=lambda x: x[1])  # Sort by energy (ascending)
        
        if not valid:
            ax.text(0.5, 0.5, f'Window {w} nt: All regions have inf energy', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_xlim(0, seq_len)
            continue
        
        positions = [p for p, _ in valid]
        energies = [e for _, e in valid]
        
        # Color by energy level
        colors = plt.cm.YlOrRd(np.linspace(0.2, 0.9, len(energies)))
        
        bars = ax.barh(range(len(positions)), [1]*len(positions), 
                      left=positions, color=colors, alpha=0.8)
        
        ax.set_xlim(0, seq_len)
        ax.set_ylim(-0.5, len(positions) - 0.5)
        ax.set_yticks(range(len(positions)))
        ax.set_yticklabels([f'{p:,}' for p in positions], fontsize=8)
        ax.set_ylabel(f'{w} nt', fontsize=11, rotation=0, ha='right')
        
        # Add energy values as text
        for i, (pos, energy) in enumerate(valid):
            ax.text(pos + seq_len*0.01, i, f'{energy:.1f}', 
                   va='center', fontsize=8)
        
        if w == window_sizes[-1]:
            ax.set_xlabel('Genome position (nt)', fontsize=11)
    
    plt.suptitle(f"Top Accessible Regions by Window Size\n{data['sequence_id']}", 
                 fontsize=13, y=1.00)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {output_path}")
    print(f"Saved: {pdf_path}")
    
    plt.close()


def generate_latex_table(data, output_path):
    """Generate LaTeX table for paper."""
    top_regions = data['top_accessible_regions']
    window_sizes = sorted([int(w) for w in top_regions.keys()])
    
    lines = [
        r'\begin{table}[ht]',
        r'\centering',
        r'\caption{Top 5 most accessible regions in SARS-CoV-2 genome by window size}',
        r'\label{tab:sars_cov2_accessibility}',
        r'\begin{tabular}{clcc}',
        r'\hline',
        r'Window size & Position & Region (\%) & Opening energy (kcal/mol) \\',
        r'\hline'
    ]
    
    for w in window_sizes:
        regions = top_regions[str(w)]
        valid = [(r['position'], r['accessibility']) for r in regions 
                if np.isfinite(r['accessibility'])]
        valid.sort(key=lambda x: x[1])
        
        for i, (pos, energy) in enumerate(valid[:5], 1):
            pct = (pos / data['sequence_length']) * 100
            region = '5\' UTR' if pct < 5 else ('3\' UTR' if pct > 95 else 'ORF')
            lines.append(f'{w} nt & {pos:,} & {region} ({pct:.1f}\%) & {energy:.2f} \\\\')
        
        if w != window_sizes[-1]:
            lines.append(r'\hline')
    
    lines.extend([
        r'\hline',
        r'\end{tabular}',
        r'\end{table}'
    ])
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Create improved SARS-CoV-2 visualizations')
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=Path('manuscript/figures'))
    
    args = parser.parse_args()
    
    with open(args.input) as f:
        data = json.load(f)
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Creating visualizations...")
    
    # Heatmap (no binning, only 3 window sizes)
    print("1. Genome heatmap (no binning, 3/10/30 nt)...")
    plot_genome_heatmap(data, args.output_dir / 'figure_sars_cov2_genome_accessibility.png', 
                       selected_windows=[3, 10, 30])
    
    # Bar chart
    print("2. Top regions bar chart...")
    plot_top_regions_bars(data, args.output_dir / 'figure_sars_cov2_top_regions.png')
    
    # LaTeX table
    print("3. LaTeX table...")
    generate_latex_table(data, args.output_dir / 'table_sars_cov2_accessibility.tex')
    
    print("\n✓ Done!")


if __name__ == '__main__':
    main()
