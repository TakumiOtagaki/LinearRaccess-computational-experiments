#!/usr/bin/env python3
"""Split bpRNA FASTA into chunks for parallel processing."""

from pathlib import Path
import math

def split_fasta(input_path: Path, output_dir: Path, chunks: int = 27):
    """Split FASTA into N chunks."""
    # Read all sequences
    sequences = []
    current_id = None
    current_seq = []
    
    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id and current_seq:
                    sequences.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        # Last sequence
        if current_id and current_seq:
            sequences.append((current_id, ''.join(current_seq)))
    
    print(f"Total sequences: {len(sequences)}")
    
    # Calculate chunk size
    chunk_size = math.ceil(len(sequences) / chunks)
    print(f"Chunk size: ~{chunk_size} per file")
    
    # Write chunks
    output_dir.mkdir(parents=True, exist_ok=True)
    for i in range(chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, len(sequences))
        
        chunk_file = output_dir / f"chunk_{i+1:02d}.fa"
        with open(chunk_file, 'w') as f:
            for j in range(start, end):
                seq_id, seq = sequences[j]
                f.write(f">{seq_id}\n")
                # Write sequence in lines of 80 chars
                for k in range(0, len(seq), 80):
                    f.write(seq[k:k+80] + '\n')
        
        print(f"Chunk {i+1}: {end-start} sequences -> {chunk_file}")
    
    return len(sequences)

if __name__ == '__main__':
    input_file = Path('/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/data/bpRNA_1m_90_filtered_clean.fa')
    output_dir = Path('/gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/data/bprna_chunks')
    
    total = split_fasta(input_file, output_dir, chunks=27)
    print(f"\nSplit complete: {total} sequences into 27 chunks")
