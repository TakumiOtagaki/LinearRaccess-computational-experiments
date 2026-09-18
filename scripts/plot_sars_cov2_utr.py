#!/usr/bin/env python3
"""
SARS-CoV-2 UTR accessibility visualization.
Focus on 5' UTR (1-265) and 3' UTR (29675-29903) regions.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import argparse

def parse_genbank_annotation(gb_file):
    """Parse GenBank file to extract UTR and stem-loop annotations."""
    annotations = {
        '5utr': (1, 265),
        '3utr': (29675, 29903),
        'stem_loops': []
    }
    
    with open(gb_file, 'r') as f:
        content = f.read()
    
    # Extract stem-loop features
    import re
    # Pattern: stem_loop <start>..<end>
    pattern = r'stem_loop\s+(\d+)\.\.(\d+)'
    matches = re.findall(pattern, content)
    
    for start, end in matches:
        annotations['stem_loops'].append((int(start), int(end)))
    
    return annotations

def plot_utr_region(data, region_name, start_pos, end_pos, annotations, output_dir, max_energy=100):
    """Plot accessibility heatmap for a specific UTR region."""
    
    full_data = data.get('full_accessibility_data', {})
    region_length = end_pos - start_pos + 1
    
    # Create figure with 3 subplots (one per window size)
    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)
    fig.suptitle(f'SARS-CoV-2 {region_name} Accessibility\n(Positions {start_pos}-{end_pos}, {region_length} nt)', 
                 fontsize=14, fontweight='bold')
    
    window_sizes = [3, 10, 30]
    colors = ['Reds', 'Reds', 'Reds']
    
    for idx, (ax, w) in enumerate(zip(axes, window_sizes)):
        if str(w) not in full_data:
            continue
        
        regions = full_data[str(w)]
        
        # Extract data for this UTR region only
        region_values = np.full(region_length, np.nan)
        
        for r in regions:
            pos = r['position']
            energy = r['accessibility']
            
            # Check if position is in UTR region
            if start_pos <= pos <= end_pos:
                local_pos = pos - start_pos
                if np.isfinite(energy) and energy <= max_energy:
                    region_values[local_pos] = energy
                else:
                    region_values[local_pos] = np.nan
        
        # Create custom colormap: white for NaN, red gradient for values
        cmap = plt.cm.Reds
        cmap.set_bad('white')
        
        # Plot as heatmap (1 row)
        im = ax.imshow(region_values.reshape(1, -1), aspect='auto', cmap=cmap,
                      vmin=0, vmax=max_energy,
                      extent=[start_pos, end_pos + 1, 0, 1])
        
        # Add stem-loop annotations
        for sl_start, sl_end in annotations['stem_loops']:
            if start_pos <= sl_start <= end_pos or start_pos <= sl_end <= end_pos:
                # Clip to region
                plot_start = max(sl_start, start_pos)
                plot_end = min(sl_end, end_pos)
                rect = Rectangle((plot_start, 0.1), plot_end - plot_start, 0.8,
                                linewidth=2, edgecolor='blue', facecolor='none',
                                linestyle='--', alpha=0.7)
                ax.add_patch(rect)
                # Add label for first subplot only
                if idx == 0:
                    ax.text((plot_start + plot_end) / 2, 1.1, f'SL\n{sl_start}-{sl_end}',
                           ha='center', va='bottom', fontsize=7, color='blue')
        
        ax.set_ylabel(f'{w} nt', fontsize=11, rotation=0, ha='right', va='center')
        ax.set_yticks([])
        ax.set_ylim(-0.2, 1.3)
        
        # Colorbar for first subplot
        if idx == 0:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.01, shrink=0.8)
            cbar.set_label('Opening energy\n(kcal/mol)', fontsize=9)
    
    axes[-1].set_xlabel('Genome position (nt)', fontsize=12)
    
    plt.tight_layout()
    
    # Save
    png_path = f'{output_dir}/figure_sars_cov2_{region_name.lower().replace(" ", "_")}_accessibility.png'
    pdf_path = f'{output_dir}/figure_sars_cov2_{region_name.lower().replace(" ", "_")}_accessibility.pdf'
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f'Saved: {png_path}')
    print(f'Saved: {pdf_path}')
    
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot SARS-CoV-2 UTR accessibility')
    parser.add_argument('--input', type=str, required=True,
                       help='Input JSON file from analyze_sars_cov2_accessibility.py')
    parser.add_argument('--annotation', type=str, 
                       default='benchmark/data/annotation/NC_045512.2.gb',
                       help='GenBank annotation file')
    parser.add_argument('--output-dir', type=str, default='manuscript/figures',
                       help='Output directory for figures')
    parser.add_argument('--max-energy', type=float, default=100,
                       help='Maximum energy value for colormap')
    
    args = parser.parse_args()
    
    # Load data
    with open(args.input, 'r') as f:
        data = json.load(f)
    
    # Parse annotation
    annotations = parse_genbank_annotation(args.annotation)
    print(f"Parsed annotation:")
    print(f"  5' UTR: {annotations['5utr']}")
    print(f"  3' UTR: {annotations['3utr']}")
    print(f"  Stem-loops: {len(annotations['stem_loops'])} found")
    for sl in annotations['stem_loops']:
        print(f"    {sl[0]}-{sl[1]}")
    
    # Plot 5' UTR
    print("\nPlotting 5' UTR...")
    plot_utr_region(data, "5' UTR", annotations['5utr'][0], annotations['5utr'][1],
                   annotations, args.output_dir, args.max_energy)
    
    # Plot 3' UTR
    print("\nPlotting 3' UTR...")
    plot_utr_region(data, "3' UTR", annotations['3utr'][0], annotations['3utr'][1],
                   annotations, args.output_dir, args.max_energy)
    
    print("\n✓ Done!")

if __name__ == '__main__':
    main()
