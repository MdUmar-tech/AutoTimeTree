#!/usr/bin/env python3
"""
AutoTimeTree - Step 1: Sequence Cleaning & MAFFT Alignment
Reads an input FASTA file, cleans sequence headers, and executes MAFFT multiple sequence alignment.
Outputs are saved to the results/ folder.
"""
import sys
import os
import argparse
import subprocess

def parse_fasta(filename):
    """Parse FASTA records from a file."""
    records = []
    header = None
    seq_lines = []
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    records.append((header, ''.join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            records.append((header, ''.join(seq_lines)))
    return records

def clean_header(header):
    """Sanitize FASTA header to avoid issues with phylogenetic tools."""
    clean = header.replace(' ', '_').replace('.', '_').replace(':', '_').replace(',', '_')
    while '__' in clean:
        clean = clean.replace('__', '_')
    return clean

def main():
    parser = argparse.ArgumentParser(description="AutoTimeTree Step 1: Clean headers and run MAFFT alignment.")
    parser.add_argument("--input", "-i", default="data/example_input.fasta", help="Path to input FASTA file (including sequences + outgroup)")
    parser.add_argument("--output", "-o", default="results/aligned_sequences.fasta", help="Path for aligned output FASTA file")
    parser.add_argument("--cleaned", "-c", default="results/cleaned_sequences.fasta", help="Path for cleaned FASTA file prior to alignment")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[Error] Input FASTA file '{args.input}' not found!")
        sys.exit(1)

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    cleaned_dir = os.path.dirname(args.cleaned)
    if cleaned_dir and not os.path.exists(cleaned_dir):
        os.makedirs(cleaned_dir, exist_ok=True)

    print(f"[*] Reading sequences from '{args.input}'...")
    raw_records = parse_fasta(args.input)
    print(f"[+] Loaded {len(raw_records)} sequences.")

    cleaned_records = []
    for header, seq in raw_records:
        cid = clean_header(header)
        cleaned_records.append((cid, seq.upper()))

    with open(args.cleaned, "w") as out_f:
        for cid, seq in cleaned_records:
            out_f.write(f">{cid}\n{seq}\n")

    print(f"[+] Cleaned sequences saved to '{args.cleaned}'.")

    # Run MAFFT alignment
    print(f"[*] Running MAFFT alignment (--auto) on '{args.cleaned}'...")
    cmd = f"mafft --auto '{args.cleaned}' > '{args.output}'"
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"[+] Multiple Sequence Alignment successfully saved to '{args.output}'.")
    except subprocess.CalledProcessError as e:
        print(f"[Error] MAFFT alignment failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
