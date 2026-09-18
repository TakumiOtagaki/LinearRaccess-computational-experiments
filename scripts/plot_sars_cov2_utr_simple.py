#!/usr/bin/env python3
"""
SARS-CoV-2 UTR accessibility visualization - Simple version.
Shows top 3 accessible regions WITHIN UTRs for each window size.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import argparse

def find_top_regions_in_utr(full_data, window_size, start_pos, end_pos, max_energy=100, n_top=3):
    """Find top n accessible regions within a specific UTR."""
    if str(window_size) not in full_data:
        return []
    
    regions = full_data[str(window_size)]
    
    # Filter to UTR region
    utr_regions = []
    for r in regions:
        pos = r['position']
        energy = r['accessibility']
        # Region is [pos, pos + window_size - 1]
        region_end = pos + window_size - 1
        
        # Check if region overlaps with UTR
        if pos <= end_pos and region_end >= start_pos:
            if np.isfinite(energy) and energy <= max_energy:
                utr_regions.append((pos, energy))
    
    # Sort by energy (lower = more accessible) and take top n
    utr_regions.sort(key=lambda x: x[1])
    return utr_regions[:n_top]

def plot_utr_top_regions(data, region_name, start_pos, end_pos, output_dir, max_energy=100):
    """Plot top 3 accessible regions for each window size as rectangles."""
    
    full_data = data.get('full_accessibility_data', {})
    region_length = end_pos - start_pos + 1
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 5))
    
    window_sizes = [3, 10, 30]
    colors = ['#d62728', '#ff7f0e', '#2ca02c']  # Red, orange, green
    y_positions = [2.5, 1.5, 0.5]  # Vertical positions for each window size
    
    # Plot each window size
    for idx, (w, y_pos, color) in enumerate(zip(window_sizes, y_positions, colors)):
        # Get top 3 regions WITHIN this UTR
        top3 = find_top_regions_in_utr(full_data, w, start_pos, end_pos, max_energy, 3)
        
        print(f"  {w} nt window: found {len(top3)} regions")
        
        if len(top3) == 0:
            ax.text((start_pos + end_pos) / 2, y_pos, f'No accessible regions (<{max_energy} kcal/mol)',
                   ha='center', va='center', fontsize=10, style='italic', color='gray')
            continue
        
        # Draw rectangles for top 3
        for rank, (pos, energy) in enumerate(top3):
            # Rectangle spans from pos to pos + w - 1
            rect_start = max(pos, start_pos)
            rect_end = min(pos + w - 1, end_pos)
            rect_width = rect_end - rect_start + 1
            
            # Draw rectangle
            alpha = 0.5 + 0.15 * (2 - rank)  # Higher rank = more opaque
            rect = Rectangle((rect_start, y_pos - 0.3), rect_width, 0.6,
                           linewidth=2, edgecolor=color, facecolor=color,
                           alpha=alpha)
            ax.add_patch(rect)
            
            # Add opening energy value
            ax.text(rect_start + rect_width/2, y_pos, f'{energy:.1f}',
                   ha='center', va='center', fontsize=9, fontweight='bold',
                   color='white' if alpha > 0.6 else 'black')
            
            # Add position label above
            ax.text(rect_start + rect_width/2, y_pos + 0.45, f'{pos}',
                   ha='center', va='bottom', fontsize=8, color=color)
            
            print(f"    #{rank+1}: pos {pos}, energy={energy:.2f} kcal/mol")
    
    # Configure axes
    margin = max(10, region_length * 0.05)
    ax.set_xlim(start_pos - margin, end_pos + margin)
    ax.set_ylim(-0.5, 3.5)
    ax.set_xlabel('Genome position (nt)', fontsize=12)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([f'{w} nt window' for w in window_sizes])
    ax.set_title(f"SARS-CoV-2 {region_name} - Top 3 Accessible Regions\n(Positions {start_pos}-{end_pos}, {region_length} nt)", 
                fontsize=14, fontweight='bold')
    
    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor=colors[0], alpha=0.5, label='3 nt window'),
        mpatches.Patch(facecolor=colors[1], alpha=0.5, label='10 nt window'),
        mpatches.Patch(facecolor=colors[2], alpha=0.5, label='30 nt window')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    # Grid
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')
    ax.axhline(y=0, color='black', linewidth=0.5)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    # Save
    png_path = f'{output_dir}/figure_sars_cov2_{region_name.lower().replace(" ", "_")}_top3.png'
    pdf_path = f'{output_dir}/figure_sars_cov2_{region_name.lower().replace(" ", "_")}_top3.pdf'
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f'  Saved: {png_path}')
    
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot SARS-CoV-2 UTR top regions')
    parser.add_argument('--input', type=str, required=True,
                       help='Input JSON file from analyze_sars_cov2_accessibility.py')
    parser.add_argument('--output-dir', type=str, default='manuscript/figures',
                       help='Output directory for figures')
    parser.add_argument('--max-energy', type=float, default=100,
                       help='Maximum energy value to consider')
    
    args = parser.parse_args()
    
    # Load data
    with open(args.input, 'r') as f:
        data = json.load(f)
    
    # Plot 5' UTR
    print("\n5' UTR (1-265):")
    plot_utr_top_regions(data, "5' UTR", 1, 265, args.output_dir, args.max_energy)
    
    # Plot 3' UTR
    print("\n3' UTR (29675-29903):")
    plot_utr_top_regions(data, "3' UTR", 29675, 29903, args.output_dir, args.max_energy)
    
    print("\n✓ Done!")

if __name__ == '__main__':
    main()
