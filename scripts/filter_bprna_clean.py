#!/usr/bin/env python3
"""Filter bpRNA sequences: length >= 50 AND only ACGTU/acgtu."""

from pathlib import Path
import re

def filter_fasta(input_path: Path, output_path: Path, min_length: int = 50):
    """Filter FASTA sequences by minimum length and valid alphabet."""
    kept = 0
    total = 0
    invalid = 0
    too_short = 0
    
    valid_pattern = re.compile(r'^[ACGTU]+$', re.IGNORECASE)
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        current_header = None
        current_seq = []
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                # Process previous sequence
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    seq = seq.upper().replace('T', 'U')
                    total += 1
                    
                    if len(seq) < min_length:
                        too_short += 1
                    elif not valid_pattern.match(seq):
                        invalid += 1
                    else:
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
            seq = seq.upper().replace('T', 'U')
            total += 1
            
            if len(seq) < min_length:
                too_short += 1
            elif not valid_pattern.match(seq):
                invalid += 1
            else:
                outfile.write(f"{current_header}\n")
                outfile.write(f"{seq}\n")
                kept += 1
    
    return total, kept, too_short, invalid

if __name__ == '__main__':
    input_file = Path('/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/data/bpRNA_1m_90.fasta')
    output_file = Path('/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/data/bpRNA_1m_90_filtered_clean.fa')
    
    print(f"Filtering {input_file}...")
    total, kept, too_short, invalid = filter_fasta(input_file, output_file, min_length=50)
    
    print(f"\n=== Results ===")
    print(f"Total sequences: {total}")
    print(f"Kept (>= 50 nt, ACGTU only): {kept}")
    print(f"Filtered out (too short): {too_short}")
    print(f"Filtered out (invalid chars): {invalid}")
    print(f"\nOutput: {output_file}")
