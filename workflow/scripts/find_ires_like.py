#!/usr/bin/env python3
"""
find_ires_like.py - Find IRES-like hexamers in circRNA sequences

This script reads circRNA sequences from a FASTA file, duplicates each sequence
to handle the circular nature of circRNAs, and searches for IRES-like hexamers.
It uses multithreading for better performance.

Usage:
    python find_ires_like.py --input circRNA_sequences.fasta --output results.csv --hexamers "AATATA,ATATAT"

Author: Norbert Moldovan
"""

import argparse
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Default IRES-like hexamers list
DEFAULT_HEXAMERS = [
    "AATATA", "AAAATA", "AAATAT", "AAATAA", "AATAAA", "ATATAA", "AAAAAA", 
    "AAATTA", "ATATAT", "AAAAAT", "ATAAAT", "ATATTA", "AATAAT", "TATATA", 
    "ATAATA", "AATTAA", "ATAAAA", "TATAAT", "ATAAGA", "AATATT", "AAAATT", 
    "AAATAC", "TAATAT", "CATATA", "ATTAAT", "ACATAT", "AATACA", "ATACAA", 
    "AATTAT", "ATATAG", "TAAATA", "TTATAA", "AGAAGA", "ATTATA", "TATAAA", 
    "AACATA", "TTAATA", "TATATT", "AATACT", "ATATAC", "AATAAG", "AAATTC", 
    "GAGATA", "CAAAAA", "AATTTA", "AAGATA", "AACATT", "ATTAGG", "ACATAA", 
    "GAAGAA", "TATACT", "AGATAT", "TCAAGC", "AAGAAT", "AAACAT", "ATTATT", 
    "ACAAAA", "AAAAGA", "AATCAA", "AAAGAC", "TAAGAA", "ATAAAC", "TAGATT", 
    "ATAAAG", "AATATC", "TAATAA", "ATTCGA", "TATTTT", "TAATTA", "TATATG", 
    "GGAGAT", "TAATCT", "TAAAAA", "AAATCC", "ATCAAG", "ATACTG", "CATTAG", 
    "TGACAT", "ATTTAA", "AGATTA", "TAAACA", "CGAAAC", "TATTAA", "AATAGA", 
    "AATTCA", "ATAAGT", "AAACAA", "ATACTA", "ATATCT", "AAGAAG", "TATACA", 
    "GACATA", "TGAATA", "TAAGAC", "AACTGA", "TTATAT", "TTTAAA"
]

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Find IRES-like hexamers in circRNA sequences')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file with circRNA sequences')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file for results')
    parser.add_argument('-x', '--hexamers', help='Comma-separated list of IRES-like hexamers (e.g., "AATATA,ATATAT")')
    parser.add_argument('-f', '--hexamers-file', help='File containing IRES-like hexamers (one per line)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to use (default: 4)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    return parser.parse_args()

def read_hexamers_from_file(file_path):
    """Read hexamers from a file, one per line."""
    try:
        with open(file_path, 'r') as f:
            hexamers = [line.strip().upper() for line in f if line.strip()]
        logger.info(f"Read {len(hexamers)} hexamers from {file_path}")
        return hexamers
    except Exception as e:
        logger.error(f"Error reading hexamers from {file_path}: {e}")
        return DEFAULT_HEXAMERS

def parse_hexamers_string(hexamers_str):
    """Parse a comma-separated string of hexamers."""
    if not hexamers_str:
        return DEFAULT_HEXAMERS
    
    try:
        # Split by comma and strip whitespace
        hexamers = [h.strip().upper() for h in hexamers_str.split(',') if h.strip()]
        logger.info(f"Parsed {len(hexamers)} hexamers from input string")
        return hexamers
    except Exception as e:
        logger.error(f"Error parsing hexamers string: {e}")
        return DEFAULT_HEXAMERS

def find_hexamer_matches(record, hexamers, min_length=6):
    """
    Find all hexamer matches in a circRNA sequence.
    
    Args:
        record: SeqIO record object
        hexamers: List of hexamer sequences to search for
        min_length: Minimum length of hexamers (default: 6)
    
    Returns:
        List of dictionaries with match information
    """
    seq_id = record.id
    seq = str(record.seq).upper()
    seq_len = len(seq)
    
    # Duplicate the sequence to handle circular nature
    extended_seq = seq + seq
    
    results = []
    
    # Create a compiled regex pattern for each hexamer for efficiency
    for hexamer in hexamers:
        if len(hexamer) < min_length:
            continue
            
        pattern = re.compile(hexamer)
        
        # Find all matches in the extended sequence
        for match in pattern.finditer(extended_seq):
            start_pos = match.start()
            end_pos = match.end() - 1  # Convert to inclusive end
            
            # Determine if the match spans the BSJ
            spans_bsj = start_pos < seq_len and end_pos >= seq_len
            
            # Adjust positions for circular representation
            if start_pos >= seq_len:
                start_pos -= seq_len
            if end_pos >= seq_len:
                end_pos -= seq_len
                
            # Create result entry
            result = {
                'name': seq_id,
                'hexamer': hexamer,
                'start': start_pos + 1,  # Convert to 1-based
                'end': end_pos + 1,      # Convert to 1-based
                'spans_bsj': spans_bsj
            }
            results.append(result)
    
    return results

def process_sequence(record, hexamers):
    """Process a single sequence record (for parallel execution)."""
    try:
        return find_hexamer_matches(record, hexamers)
    except Exception as e:
        logger.error(f"Error processing sequence {record.id}: {e}")
        return []

def main():
    """Main function."""
    args = parse_arguments()
    
    # Set up logging verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Load hexamers
    if args.hexamers:
        hexamers = parse_hexamers_string(args.hexamers)
    elif args.hexamers_file:
        hexamers = read_hexamers_from_file(args.hexamers_file)
    else:
        logger.info("Using default hexamer list")
        hexamers = DEFAULT_HEXAMERS
    
    logger.info(f"Loaded {len(hexamers)} IRES-like hexamers")
    
    # Check if input file exists
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Read sequences
    try:
        records = list(SeqIO.parse(args.input, "fasta"))
        logger.info(f"Read {len(records)} sequences from {args.input}")
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        sys.exit(1)
    
    # Process sequences in parallel
    all_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_sequence, record, hexamers) for record in records]
        
        for future in futures:
            results = future.result()
            all_results.extend(results)
    
    logger.info(f"Found {len(all_results)} hexamer matches")
    
    # Convert results to DataFrame for easier manipulation
    if all_results:
        df = pd.DataFrame(all_results)
        
        # Count hexamers per sequence
        hexamer_counts = df.groupby('name').size().reset_index(name='hexamer_count')
        
        # Add sequence length information
        seq_lengths = {record.id: len(record.seq) for record in records}
        hexamer_counts['seq_length'] = hexamer_counts['name'].map(seq_lengths)
        
        # Calculate hexamer density (hexamers per kb)
        hexamer_counts['hexamer_density'] = hexamer_counts['hexamer_count'] * 1000 / hexamer_counts['seq_length']
        
        # Write results to CSV
        df.to_csv(args.output, index=False)
        
        # Write summary to a separate file
        summary_file = args.output.replace('.csv', '_summary.csv')
        hexamer_counts.to_csv(summary_file, index=False)
        
        logger.info(f"Results written to {args.output}")
        logger.info(f"Summary written to {summary_file}")
    else:
        # Create empty output file
        with open(args.output, 'w') as f:
            f.write("name,hexamer,start,end,spans_bsj\n")
        logger.warning(f"No hexamer matches found. Empty output file created: {args.output}")

if __name__ == "__main__":
    main() 