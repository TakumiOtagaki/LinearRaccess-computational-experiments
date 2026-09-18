#!/usr/bin/env python3
"""Filter bpRNA sequences by length (>= 50 nt)."""

from pathlib import Path
import sys

def filter_fasta(input_path: Path, output_path: Path, min_length: int = 50):
    """Filter FASTA sequences by minimum length."""
    kept = 0
    total = 0
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        current_header = None
        current_seq = []
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    seq = seq.replace('T', 'U').replace('t', 'u').upper()
                    total += 1
                    if len(seq) >= min_length:
                        outfile.write(f"{current_header}\n")
                        outfile.write(f"{seq}\n")
                        kept += 1
                
                # Start new sequence
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Process last sequence
        if current_header and current_seq:
            seq = ''.join(current_seq)
            seq = seq.replace('T', 'U').replace('t', 'u').upper()
            total += 1
            if len(seq) >= min_length:
                outfile.write(f"{current_header}\n")
                outfile.write(f"{seq}\n")
                kept += 1
    
    return total, kept

if __name__ == '__main__':
    input_file = Path('/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/data/bpRNA_1m_90.fasta')
    output_file = Path('/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/data/bpRNA_1m_90_filtered_50nt.fa')
    
    print(f"Filtering {input_file}...")
    total, kept = filter_fasta(input_file, output_file, min_length=50)
    
    print(f"Total sequences: {total}")
    print(f"Kept (>= 50 nt): {kept}")
    print(f"Filtered out: {total - kept}")
    print(f"Output written to: {output_file}")
