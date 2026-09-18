#!/usr/bin/env python3
"""
Run Raccess on a single sequence and record time/memory.
"""
import argparse
import json
import subprocess
import tempfile
import os
import re
from pathlib import Path

RT = 0.61633008  # kcal/mol, from Raccess documentation

def energy_to_prob(energy):
    """Convert opening energy to probability."""
    import math
    return math.exp(-energy / RT)

def parse_raccess_output(output_file):
    """Parse Raccess output and return accessibility values."""
    results = {}
    current_seq = None
    
    with open(output_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq = line[1:]
                results[current_seq] = {}
            elif line and current_seq:
                # Format: pos\tlen1,ene1;len2,ene2;
                parts = line.split('\t')
                if len(parts) == 2:
                    pos = int(parts[0])
                    entries = parts[1].rstrip(';').split(';')
                    for entry in entries:
                        if ',' in entry:
                            length, energy = entry.split(',')
                            key = f"w{length}"
                            if key not in results[current_seq]:
                                results[current_seq][key] = {}
                            results[current_seq][key][pos] = energy_to_prob(float(energy))
    
    return results

def run_raccess(seq_file, seq_id, seq_length, access_len, max_span, output_dir, method_name):
    """Run Raccess and return timing/memory results."""
    
    output_file = os.path.join(output_dir, f"{seq_id}_{method_name}_out.txt")
    
    # Build command
    cmd = [
        "vendor/raccess/src/raccess/run_raccess",
        f"-seqfile={seq_file}",
        f"-access_len={access_len}",
        f"-max_span={max_span}",
        f"-outfile={output_file}"
    ]
    
    # Run with time
    time_cmd = ["/usr/bin/time", "-v"] + cmd
    
    try:
        result = subprocess.run(
            time_cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        # Parse time output
        stderr = result.stderr
        
        # Extract elapsed time
        elapsed_match = re.search(r'Elapsed \(wall clock\) time.*?:(\d+):([\d.]+)', stderr)
        if elapsed_match:
            minutes = int(elapsed_match.group(1))
            seconds = float(elapsed_match.group(2))
            elapsed_sec = minutes * 60 + seconds
        else:
            elapsed_match = re.search(r'Elapsed \(wall clock\) time.*?(\d+):(\d+):([\d.]+)', stderr)
            if elapsed_match:
                hours = int(elapsed_match.group(1))
                minutes = int(elapsed_match.group(2))
                seconds = float(elapsed_match.group(3))
                elapsed_sec = hours * 3600 + minutes * 60 + seconds
            else:
                elapsed_sec = None
        
        # Extract memory
        memory_match = re.search(r'Maximum resident set size.*?:(\d+)', stderr)
        memory_kb = int(memory_match.group(1)) if memory_match else None
        
        # Check output
        output_path = f"{output_file}/raccess_out.txt"
        if not os.path.exists(output_path):
            # Try alternative location
            output_path = f"{output_file}"
        
        status = "SUCCESS" if result.returncode == 0 and os.path.exists(output_path) else "FAILED"
        
        return {
            "seq_id": seq_id,
            "seq_length": seq_length,
            "method": method_name,
            "access_len": access_len,
            "max_span": max_span,
            "elapsed_sec": elapsed_sec,
            "memory_kb": memory_kb,
            "status": status,
            "exit_code": result.returncode
        }
        
    except subprocess.TimeoutExpired:
        return {
            "seq_id": seq_id,
            "seq_length": seq_length,
            "method": method_name,
            "access_len": access_len,
            "max_span": max_span,
            "elapsed_sec": None,
            "memory_kb": None,
            "status": "TIMEOUT",
            "exit_code": -1
        }
    except Exception as e:
        return {
            "seq_id": seq_id,
            "seq_length": seq_length,
            "method": method_name,
            "access_len": access_len,
            "max_span": max_span,
            "elapsed_sec": None,
            "memory_kb": None,
            "status": f"ERROR: {str(e)}",
            "exit_code": -1
        }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seq-file", required=True)
    parser.add_argument("--seq-id", required=True)
    parser.add_argument("--seq-length", type=int, required=True)
    parser.add_argument("--access-len", default="7")
    parser.add_argument("--max-span", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--output-json", required=True)
    
    args = parser.parse_args()
    
    result = run_raccess(
        args.seq_file,
        args.seq_id,
        args.seq_length,
        args.access_len,
        args.max_span,
        args.output_dir,
        args.method
    )
    
    with open(args.output_json, "w") as f:
        json.dump(result, f, indent=2)
    
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
