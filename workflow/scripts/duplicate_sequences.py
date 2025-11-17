#!/usr/bin/env python3
"""
Duplicate sequences to account for circular nature of circRNAs.
Concatenates each sequence with itself (seq + seq) to handle circularity.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import multiprocessing as mp
from functools import partial


def duplicate_sequence(record):
    """Duplicate a single sequence record."""
    # Create new sequence by concatenating with itself
    duplicated_seq = record.seq + record.seq
    # Create new record with duplicated sequence
    new_record = record.__class__(
        duplicated_seq,
        id=record.id,
        description=record.description
    )
    return new_record


def process_sequences_parallel(input_file, output_file, threads):
    """Process sequences in parallel."""
    print(f"Reading sequences from {input_file}")
    records = list(SeqIO.parse(input_file, "fasta"))
    print(f"Found {len(records)} sequences")
    
    # Use multiprocessing to duplicate sequences
    with mp.Pool(threads) as pool:
        duplicated_records = list(tqdm(
            pool.imap(duplicate_sequence, records),
            total=len(records),
            desc="Duplicating sequences"
        ))
    
    # Write duplicated sequences to output file
    print(f"Writing duplicated sequences to {output_file}")
    SeqIO.write(duplicated_records, output_file, "fasta")
    print(f"Successfully wrote {len(duplicated_records)} duplicated sequences")


def main():
    parser = argparse.ArgumentParser(
        description="Duplicate circRNA sequences to account for circular nature"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file with circRNA sequences"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file with duplicated sequences"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=mp.cpu_count(),
        help="Number of threads to use (default: all available cores)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    try:
        process_sequences_parallel(args.input, args.output, args.threads)
        if args.verbose:
            print("Sequence duplication completed successfully")
    except Exception as e:
        print(f"Error processing sequences: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
