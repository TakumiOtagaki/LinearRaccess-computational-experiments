#!/usr/bin/env python3
"""Analyze SARS-CoV-2 genome accessibility and identify high-accessibility regions."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import json
import re
from collections import defaultdict
import numpy as np


def parse_fasta(filepath):
    """Parse single-sequence FASTA file."""
    with open(filepath) as f:
        lines = f.readlines()
    
    header = lines[0].strip()
    sequence = ''.join(line.strip() for line in lines[1:])
    seq_id = header[1:].split()[0] if header.startswith('>') else header
    return seq_id, sequence


def run_linearraccess(sequence, seq_id, linracc_bin, access_lengths, beam=100, 
                      c_multi=30, c_hairpin=30, timeout=600):
    """Run LinearRaccess and parse output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write input
        input_fa = Path(tmpdir) / 'input.fa'
        input_fa.write_text(f'>{seq_id}\n{sequence}\n')
        
        output_file = Path(tmpdir) / 'output.txt'
        
        # Build command
        cmd = [
            str(linracc_bin),
            f'-seqfile={input_fa}',
            f'-outfile={output_file}',
            f'-access_len={",".join(map(str, access_lengths))}',
            '-engine=lincapr',
            '-energy=turner2004',
            f'-beam={beam}',
            f'-c_multi={c_multi}',
            f'-c_hairpin={c_hairpin}',
            '-metadata'
            # Note: Default output is opening energy (-RT log P), not probability
        ]
        
        # Run
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        if result.returncode != 0:
            raise RuntimeError(f"LinearRaccess failed: {result.stderr}")
        
        # Parse output
        return parse_accessibility_output(output_file)


def parse_accessibility_output(filepath):
    """Parse LinearRaccess output file (opening energy format)."""
    accessibility = defaultdict(dict)  # access_len -> {position -> energy}
    
    with open(filepath) as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('>'):
            continue
        
        parts = line.split('\t')
        if len(parts) != 2:
            continue
        
        pos = int(parts[0])
        # Parse accessibility values: "len1,energy1;len2,energy2;..."
        # Energy is -RT log(P) in kcal/mol
        entries = parts[1].rstrip(';').split(';')
        for entry in entries:
            if ',' in entry:
                w_str, e_str = entry.split(',')
                w = int(w_str)
                e = float(e_str)
                accessibility[w][pos] = e
    
    return accessibility


def find_high_accessibility_regions(accessibility_data, top_n=10):
    """Find top N high accessibility regions for each window size."""
    results = {}
    
    for window_size, pos_probs in accessibility_data.items():
        # Sort by accessibility value descending
        sorted_regions = sorted(pos_probs.items(), key=lambda x: x[1], reverse=True)
        results[window_size] = sorted_regions[:top_n]
    
    return results


def analyze_sliding_windows(accessibility_data, sequence_length, window_size=100):
    """Analyze accessibility in sliding windows across the genome."""
    # For each position, get average accessibility across all window sizes
    avg_accessibility = {}
    
    for pos in range(sequence_length):
        values = []
        for w, pos_dict in accessibility_data.items():
            if pos in pos_dict:
                values.append(pos_dict[pos])
        if values:
            avg_accessibility[pos] = np.mean(values)
    
    # Find high-accessibility sliding windows
    sliding_windows = []
    step = window_size // 2
    
    for start in range(0, sequence_length - window_size + 1, step):
        end = start + window_size
        window_values = [avg_accessibility.get(i, 0) for i in range(start, end)]
        if window_values:
            avg_val = np.mean(window_values)
            max_val = max(window_values)
            sliding_windows.append({
                'start': start,
                'end': end,
                'mean_accessibility': avg_val,
                'max_accessibility': max_val,
                'length': window_size
            })
    
    # Sort by mean accessibility
    sliding_windows.sort(key=lambda x: x['mean_accessibility'], reverse=True)
    return sliding_windows


def main():
    parser = argparse.ArgumentParser(
        description='Analyze SARS-CoV-2 genome accessibility profiles'
    )
    parser.add_argument('--input', type=Path,
                        default=Path('benchmark/data/NC_045512.2.fa'),
                        help='Input SARS-CoV-2 FASTA file')
    parser.add_argument('--linracc', type=Path,
                        default=Path('vendor/LinearRaccess/build-raccess/LinRacc'),
                        help='Path to LinearRaccess binary')
    parser.add_argument('--output', type=Path,
                        default=Path('results/sars_cov2_accessibility_analysis.json'),
                        help='Output JSON file')
    parser.add_argument('--access-lens', type=int, nargs='+',
                        default=[1, 3, 7, 10, 20, 30, 50],
                        help='Window lengths to analyze')
    parser.add_argument('--beam', type=int, default=100,
                        help='Beam width')
    parser.add_argument('--top-n', type=int, default=10,
                        help='Number of top regions to report per window size')
    parser.add_argument('--sliding-window', type=int, default=100,
                        help='Sliding window size for region analysis')
    
    args = parser.parse_args()
    
    # Parse input
    print(f"Reading {args.input}...")
    seq_id, sequence = parse_fasta(args.input)
    seq_len = len(sequence)
    print(f"Sequence: {seq_id}")
    print(f"Length: {seq_len:,} nt")
    
    # Run LinearRaccess
    print(f"\nRunning LinearRaccess (beam={args.beam})...")
    print(f"Window sizes: {args.access_lens}")
    
    try:
        accessibility_data = run_linearraccess(
            sequence, seq_id, args.linracc, args.access_lens, 
            beam=args.beam
        )
        print(f"✓ Computed accessibility for {len(accessibility_data)} window sizes")
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    # Find high accessibility regions
    print(f"\nFinding top {args.top_n} high-accessibility regions...")
    top_regions = find_high_accessibility_regions(accessibility_data, args.top_n)
    
    # Analyze sliding windows
    print(f"Analyzing sliding windows (size={args.sliding_window})...")
    sliding_results = analyze_sliding_windows(
        accessibility_data, seq_len, args.sliding_window
    )
    
    # Prepare output
    results = {
        'sequence_id': seq_id,
        'sequence_length': seq_len,
        'beam_width': args.beam,
        'parameters': {
            'access_lengths': args.access_lens,
            'c_multi': 30,
            'c_hairpin': 30,
            'energy_model': 'turner2004'
        },
        'full_accessibility_data': {
            # Save ALL positions for each window size
            str(w): [
                {'position': pos, 'accessibility': energy}
                for pos, energy in sorted(accessibility_data[w].items())
            ]
            for w in accessibility_data
        },
        'top_accessible_regions': {
            str(w): [
                {'position': pos, 'accessibility': prob}
                for pos, prob in regions
            ]
            for w, regions in top_regions.items()
        },
        'sliding_window_analysis': {
            'window_size': args.sliding_window,
            'top_regions': sliding_results[:20]  # Top 20 sliding windows
        },
        'statistics': {
            'window_sizes_analyzed': len(accessibility_data),
            'total_positions': seq_len
        }
    }
    
    # Save results
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✓ Results saved to: {args.output}")
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY: Top Accessible Regions")
    print("="*60)
    
    for w in sorted(top_regions.keys()):
        regions = top_regions[w]
        print(f"\nWindow size {w} nt:")
        for i, (pos, prob) in enumerate(regions[:5], 1):
            print(f"  {i}. Position {pos:,}: accessibility = {prob:.6f}")
    
    print("\n" + "="*60)
    print("Sliding Window Analysis")
    print("="*60)
    print(f"\nTop 5 {args.sliding_window}-nt windows by mean accessibility:")
    for i, win in enumerate(sliding_results[:5], 1):
        print(f"  {i}. {win['start']:,}-{win['end']:,}: "
              f"mean={win['mean_accessibility']:.6f}, "
              f"max={win['max_accessibility']:.6f}")
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
