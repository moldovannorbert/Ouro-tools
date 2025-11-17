#!/usr/bin/env python3
"""
merge_ires_predictions.py - Merge IRES and IRES-like predictions with CPC2 ORF results

This script reads IRES and IRES-like predictions along with circRNA sequences and CPC2 results,
filters IRES/IRES-like elements that are upstream of CPC2-predicted ORFs,
and calculates the distance between IRES/IRES-like elements and the ORF start.

Usage:
    python merge_ires_predictions.py --ires ires_results.csv --ireslike ireslike_results.csv 
                                    --fasta circRNA_sequences.fasta --cpc2 cpc2_results.tsv
                                    --output merged_results.csv --summary orf_summary.csv --threads 8

Author: Norbert Moldovan
"""

import argparse
import re
import sys
import logging
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Merge IRES and IRES-like predictions with CPC2 ORF results')
    parser.add_argument('--ires', required=True, help='Input CSV file with IRES predictions')
    parser.add_argument('--ireslike', required=True, help='Input CSV file with IRES-like predictions')
    parser.add_argument('--fasta', required=True, help='Input FASTA file with circRNA sequences')
    parser.add_argument('--cpc2', required=True, help='Input TSV file with CPC2 results')
    parser.add_argument('--output', required=True, help='Output CSV file for merged results')
    parser.add_argument('--summary', required=True, help='Output CSV file for CDS summary')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    return parser.parse_args()

def read_ires_predictions(ires_file):
    """Read IRES predictions from CSV file."""
    try:
        df = pd.read_csv(ires_file)
        logger.info(f"Read {len(df)} IRES predictions from {ires_file}")
        return df
    except Exception as e:
        logger.error(f"Error reading IRES predictions from {ires_file}: {e}")
        return pd.DataFrame()

def read_ireslike_predictions(ireslike_file):
    """Read IRES-like predictions from CSV file."""
    try:
        df = pd.read_csv(ireslike_file)
        logger.info(f"Read {len(df)} IRES-like predictions from {ireslike_file}")
        return df
    except Exception as e:
        logger.error(f"Error reading IRES-like predictions from {ireslike_file}: {e}")
        return pd.DataFrame()

def read_sequences(fasta_file):
    """Read circRNA sequences from FASTA file."""
    try:
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq).upper()
        logger.info(f"Read {len(sequences)} sequences from {fasta_file}")
        return sequences
    except Exception as e:
        logger.error(f"Error reading sequences from {fasta_file}: {e}")
        return {}

def read_cpc2_results(cpc2_file):
    """Read CPC2 results from TSV file."""
    try:
        df = pd.read_csv(cpc2_file, sep='\t')
        logger.info(f"Read {len(df)} CPC2 results from {cpc2_file}")
        
        # Create a dictionary mapping circRNA names to ORF_Start positions
        cpc2_dict = {}
        for _, row in df.iterrows():
            circrna_name = row['#ID']
            orf_start = row['ORF_Start']
            cpc2_dict[circrna_name] = orf_start
            
        logger.info(f"Created CPC2 dictionary with {len(cpc2_dict)} entries")
        return cpc2_dict
    except Exception as e:
        logger.error(f"Error reading CPC2 results from {cpc2_file}: {e}")
        return {}


def get_ires_sequence(seq, start, end):
    """Extract IRES sequence from duplicated circRNA sequence."""
    # For duplicated sequences, we can directly extract the sequence
    # since the duplication already handles the circular nature
    return seq[start-1:end]

def process_circrna(circrna_id, sequence, ires_df, ireslike_df, cpc2_dict):
    """
    Process a single circRNA to find IRES elements upstream of CPC2-predicted ORFs.
    
    Args:
        circrna_id: The circRNA identifier
        sequence: The circRNA sequence
        ires_df: DataFrame with IRES predictions
        ireslike_df: DataFrame with IRES-like predictions
        cpc2_dict: Dictionary mapping circRNA names to ORF_Start positions
    
    Returns:
        List of dictionaries with merged results
    """
    try:
        # Get ORF start position from CPC2 results
        orf_start = cpc2_dict.get(circrna_id, 0)
        
        # Skip circRNAs with no ORF found (ORF_Start == 0)
        if orf_start == 0:
            return []
            
        # Filter IRES predictions for this circRNA
        ires_preds = ires_df[ires_df['name'] == circrna_id] if not ires_df.empty else pd.DataFrame()
        
        # Filter IRES-like predictions for this circRNA
        ireslike_preds = ireslike_df[ireslike_df['name'] == circrna_id] if not ireslike_df.empty else pd.DataFrame()
        
        results = []
        
        # Process IRES predictions
        for _, ires in ires_preds.iterrows():
            ires_start = int(ires['start'])
            ires_end = int(ires['stop'])
            ires_seq = get_ires_sequence(sequence, ires_start, ires_end)
            
            # Check if IRES is upstream of ORF start in the duplicated sequence
            # Since we're using duplicated sequences, we can directly compare positions
            # IRES is upstream if it ends before ORF start in the duplicated sequence
            if ires_end < orf_start:
                distance = orf_start - ires_end
            else:
                # IRES is after ORF in the duplicated sequence, so not upstream
                continue
            
            # Only keep IRES elements that are upstream (end before ORF start)
            result = {
                'name': circrna_id,
                'type': 'IRES',
                'start': ires_start,
                'end': ires_end,
                'sequence': ires_seq,
                'score': ires['score'] if 'score' in ires.index else None,
                'orf_start': orf_start,
                'distance_to_orf': distance,
                'spans_bsj': ires_start > ires_end or (ires_start <= len(sequence)//2 and ires_end > len(sequence)//2)
            }
            results.append(result)
        
        # Process IRES-like predictions
        for _, ireslike in ireslike_preds.iterrows():
            ireslike_start = int(ireslike['start'])
            ireslike_end = int(ireslike['end'])
            ireslike_seq = get_ires_sequence(sequence, ireslike_start, ireslike_end)
            
            # Check if IRES-like is upstream of ORF start in the duplicated sequence
            # Since we're using duplicated sequences, we can directly compare positions
            # IRES-like is upstream if it ends before ORF start in the duplicated sequence
            if ireslike_end < orf_start:
                distance = orf_start - ireslike_end
            else:
                # IRES-like is after ORF in the duplicated sequence, so not upstream
                continue
            
            # Only keep IRES-like elements that are upstream (end before ORF start)
            result = {
                'name': circrna_id,
                'type': 'IRES-like',
                'start': ireslike_start,
                'end': ireslike_end,
                'sequence': ireslike_seq,
                'hexamer': ireslike['hexamer'] if 'hexamer' in ireslike.index else None,
                'orf_start': orf_start,
                'distance_to_orf': distance,
                'spans_bsj': ireslike_start > ireslike_end or (ireslike_start <= len(sequence)//2 and ireslike_end > len(sequence)//2)
            }
            results.append(result)
        
        return results
    
    except Exception as e:
        logger.error(f"Error processing circRNA {circrna_id}: {e}")
        return []

def main():
    """Main function."""
    args = parse_arguments()
    
    # Set up logging verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Read input files
    ires_df = read_ires_predictions(args.ires)
    ireslike_df = read_ireslike_predictions(args.ireslike)
    sequences = read_sequences(args.fasta)
    cpc2_dict = read_cpc2_results(args.cpc2)
    
    if not sequences:
        logger.error("No sequences found. Exiting.")
        sys.exit(1)
    
    if not cpc2_dict:
        logger.error("No CPC2 results found. Exiting.")
        sys.exit(1)
    
    # Process each circRNA in parallel
    all_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_circrna, circrna_id, sequence, ires_df, ireslike_df, cpc2_dict) 
                  for circrna_id, sequence in sequences.items()]
        
        for future in futures:
            results = future.result()
            all_results.extend(results)
    
    logger.info(f"Found {len(all_results)} IRES/IRES-like elements upstream of ORFs")
    
    # Convert results to DataFrame
    if all_results:
        df = pd.DataFrame(all_results)
        
        # Write results to CSV
        df.to_csv(args.output, index=False)
        logger.info(f"Results written to {args.output}")
        
        # Create summary statistics
        summary = df.groupby(['name', 'type']).agg({
            'distance_to_orf': ['mean', 'min', 'max', 'count']
        }).reset_index()
        
        # Flatten multi-level columns
        summary.columns = ['_'.join(col).strip('_') for col in summary.columns.values]
        
        # Write summary to CSV
        summary.to_csv(args.summary, index=False)
        logger.info(f"Summary written to {args.summary}")
    else:
        # Create empty output files
        with open(args.output, 'w') as f:
            f.write("name,type,start,end,sequence,score,orf_start,distance_to_orf,spans_bsj\n")
        with open(args.summary, 'w') as f:
            f.write("name,type,distance_to_orf_mean,distance_to_orf_min,distance_to_orf_max,distance_to_orf_count\n")
        logger.warning(f"No IRES/IRES-like elements upstream of ORFs found. Empty output files created.")

if __name__ == "__main__":
    main() 