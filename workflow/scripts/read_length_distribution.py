#!/usr/bin/env python3

import argparse
import pandas as pd
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate read length distribution')
    parser.add_argument('-i', 
                        '--input', 
                        required=True, 
                        help='Input file with read information (can be gzipped)')
    parser.add_argument('-o', 
                        '--output', 
                        required=True, 
                        help='Output file for length distribution')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Read input file (handles both .tsv and .tsv.gz)
    if args.input.endswith('.gz'):
        with gzip.open(args.input, 'rt') as f:
            df = pd.read_csv(f, sep='\t')
    else:
        df = pd.read_csv(args.input, sep='\t')
    
    # For duplicate readIDs, keep the one with highest percentIdentity
    df = df.sort_values('percentIdentity', ascending=False)
    df = df.drop_duplicates(subset='readIDs', keep='first')
    
    # Calculate length distribution
    length_counts = df['lengths'].value_counts().reset_index()
    length_counts.columns = ['read_length', 'count']
    
    # Calculate proportions
    total_reads = length_counts['count'].sum()
    length_counts['proportion'] = length_counts['count'] / total_reads
    
    # Sort by read length
    length_counts = length_counts.sort_values('read_length')
    
    # Save to output file
    length_counts.to_csv(args.output,
                         sep='\t', 
                         index=False)

if __name__ == '__main__':
    main()



