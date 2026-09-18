#!/usr/bin/env python3
"""
Run LinearRaccess on a single sequence and record time/memory.
"""
import argparse
import json
import subprocess
import os
import re
import tempfile

def run_linearraccess(seq_file, seq_id, seq_length, access_len, beam, output_dir, method_name):
    """Run LinearRaccess and return timing/memory results."""
    
    output_file = os.path.join(output_dir, f"{seq_id}_{method_name}.txt")
    
    # Build command
    cmd = [
        "vendor/LinearRaccess/build-raccess/LinRacc",
        "--input", seq_file,
        "--access-len", str(access_len),
        "--beam", str(beam),
        "--output", output_file,
        "--log-temp"
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
        
        # Parse time output from stderr
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
        status = "SUCCESS" if result.returncode == 0 and os.path.exists(output_file) else "FAILED"
        
        return {
            "seq_id": seq_id,
            "seq_length": seq_length,
            "method": method_name,
            "access_len": access_len,
            "beam": beam,
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
            "beam": beam,
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
            "beam": beam,
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
    parser.add_argument("--access-len", type=int, default=7)
    parser.add_argument("--beam", type=int, default=200)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--output-json", required=True)
    
    args = parser.parse_args()
    
    result = run_linearraccess(
        args.seq_file,
        args.seq_id,
        args.seq_length,
        args.access_len,
        args.beam,
        args.output_dir,
        args.method
    )
    
    with open(args.output_json, "w") as f:
        json.dump(result, f, indent=2)
    
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
