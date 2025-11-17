#!/usr/bin/env python3
"""
Extract ORFs from circRNA sequences and translate them to amino acid sequences.
Takes CPC2 results and duplicated sequences as input, outputs peptides.fasta.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm


def extract_and_translate_orf(sequence, orf_start, peptide_length):
    """
    Extract ORF sequence and translate to amino acids based on peptide length.
    
    Args:
        sequence: DNA sequence (BioPython Seq or string)
        orf_start: Starting position of ORF (0-based index)
        peptide_length: Peptide length in amino acids
    
    Returns:
        Amino acid sequence, or None if extraction fails
    """
    seq_str = str(sequence)
    
    # Calculate DNA length: peptide_length amino acids * 3 nucleotides per amino acid
    # Plus 3 nucleotides for the stop codon
    dna_length = peptide_length * 3 + 3
    
    # Calculate end position
    orf_end = orf_start + dna_length
    
    # Check if we have enough sequence
    if orf_end > len(seq_str):
        return None
    
    # Extract ORF sequence (inclusive of stop codon)
    orf_seq = seq_str[orf_start:orf_end]
    
    # Translate to amino acids
    seq_obj = Seq(orf_seq)
    aa_seq = seq_obj.translate(to_stop=False)  # Don't stop early, translate all codons
    
    # Remove the stop codon (*) if present
    aa_str = str(aa_seq).rstrip('*')
    
    # Verify we got the expected length
    if len(aa_str) != peptide_length:
        return None
    
    return aa_str


def main():
    parser = argparse.ArgumentParser(
        description="Extract ORFs from circRNA sequences and translate to amino acids"
    )
    parser.add_argument(
        "--cpc2_results",
        required=True,
        help="Input CPC2 results TSV file"
    )
    parser.add_argument(
        "--sequences",
        required=True,
        help="Input duplicated sequences FASTA file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output peptides FASTA file"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not Path(args.cpc2_results).exists():
        print(f"Error: CPC2 results file {args.cpc2_results} does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.sequences).exists():
        print(f"Error: Sequences file {args.sequences} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Read CPC2 results
    print(f"Reading CPC2 results from {args.cpc2_results}")
    try:
        cpc2_df = pd.read_csv(args.cpc2_results, sep='\t')
        print(f"Found {len(cpc2_df)} total results")
    except Exception as e:
        print(f"Error reading CPC2 results: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Filter for coding sequences
    if 'label' not in cpc2_df.columns:
        print("Error: 'label' column not found in CPC2 results", file=sys.stderr)
        sys.exit(1)
    
    coding_df = cpc2_df[cpc2_df['label'] == 'coding'].copy()
    print(f"Found {len(coding_df)} coding sequences")
    
    if len(coding_df) == 0:
        print("Warning: No coding sequences found. Creating empty output file.", file=sys.stderr)
        with open(args.output, 'w') as f:
            pass
        sys.exit(0)
    
    # Read sequences into a dictionary
    print(f"Reading sequences from {args.sequences}")
    try:
        sequence_dict = {}
        for record in SeqIO.parse(args.sequences, "fasta"):
            sequence_dict[record.id] = record.seq
        print(f"Found {len(sequence_dict)} sequences")
    except Exception as e:
        print(f"Error reading sequences: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check if peptide_length column exists
    if 'peptide_length' not in coding_df.columns:
        print("Error: 'peptide_length' column not found in CPC2 results", file=sys.stderr)
        sys.exit(1)
    
    # Process each coding sequence
    print("Extracting and translating ORFs...")
    peptides = []
    errors = []
    
    for _, row in tqdm(coding_df.iterrows(), total=len(coding_df), desc="Processing"):
        seq_id = row['#ID']
        orf_start = int(row['ORF_Start'])
        peptide_length = int(row['peptide_length'])
        
        # Skip if peptide_length is 0 or invalid
        if peptide_length <= 0:
            if args.verbose:
                print(f"Warning: Invalid peptide_length for sequence {seq_id}", file=sys.stderr)
            continue
        
        # Get sequence
        if seq_id not in sequence_dict:
            if args.verbose:
                print(f"Warning: Sequence {seq_id} not found in FASTA file", file=sys.stderr)
            continue
        
        sequence = sequence_dict[seq_id]
        
        # Extract and translate ORF
        aa_seq = extract_and_translate_orf(sequence, orf_start, peptide_length)
        
        if aa_seq is None:
            if args.verbose:
                print(f"Warning: Could not extract ORF for sequence {seq_id}", file=sys.stderr)
            errors.append(seq_id)
            continue
        
        peptides.append((seq_id, aa_seq))
    
    # Write peptides to output file
    print(f"Writing {len(peptides)} peptides to {args.output}")
    try:
        with open(args.output, 'w') as f:
            for seq_id, aa_seq in peptides:
                f.write(f">{seq_id}\n{aa_seq}\n")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Successfully wrote {len(peptides)} peptides")
    
    # Report errors if any
    if errors:
        print(f"\nWarning: {len(errors)} sequences could not be processed")
        if args.verbose:
            for error in errors[:10]:  # Show first 10 errors
                print(f"  {error}", file=sys.stderr)
            if len(errors) > 10:
                print(f"  ... and {len(errors) - 10} more", file=sys.stderr)


if __name__ == "__main__":
    main()
