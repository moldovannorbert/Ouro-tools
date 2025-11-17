#!/usr/bin/env python3
"""
Script to extract circRNA sequences from a reference genome based on circRNA info file.
"""

import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Extract circRNA sequences from reference genome.')
    parser.add_argument('--circRNA_info', required=True, help='Path to circRNA info file')
    parser.add_argument('--ref_genome', required=True, help='Path to reference genome FASTA file')
    parser.add_argument('--output', required=True, help='Path to output FASTA file')
    return parser.parse_args()

def parse_isoform(isoform_str):
    """Extract exon coordinates from isoform string."""
    exons = []
    for exon in isoform_str.split(','):
        exon = exon.strip()
        if '-' in exon:
            parts = exon.split('-')
        elif '|' in exon:
            parts = exon.split('|')
        else:
            continue

        if len(parts) != 2:
            continue

        try:
            start, end = map(int, parts)
            exons.append((start, end))
        except ValueError:
            continue

    return exons

def extract_sequence(chrom, exons, genome_dict, strand):
    """Extract sequence from reference genome based on exon coordinates."""
    seq = ""
    for start, end in exons:
        # Extract sequence for this exon
        exon_seq = str(genome_dict[chrom].seq[start:end])
        seq += exon_seq
    
    # If on negative strand, return reverse complement
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    
    return seq

def main():
    """Main function to extract circRNA sequences."""
    args = parse_args()
    
    # Read reference genome
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.ref_genome, "fasta"))
    
    # Read circRNA info file and extract sequences
    output_records = []
    with open(args.circRNA_info, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            
            # Extract chromosome and strand
            chrom = fields[0]
            strand = fields[6]
            
            # Extract circ_id and isoform information using regex
            circ_id_match = re.search(r'circ_id "([^"]+)"', fields[8])
            isoform_match = re.search(r'isoform "([^"]+)"', fields[8])
            gene_name_match = re.search(r'gene_name "([^"]+)"', fields[8])
            circ_type_match = re.search(r'circ_type "([^"]+)"', fields[8])
            
            if circ_id_match and isoform_match:
                circ_id = circ_id_match.group(1)
                isoform_str = isoform_match.group(1)
                gene_name = gene_name_match.group(1) if gene_name_match else "Unknown"
                circ_type = circ_type_match.group(1) if circ_type_match else "Unknown"
                # Check if there are alternative isoform definitions separated by '|'
                if '|' in isoform_str:
                    isoform_defs = [alt.strip() for alt in isoform_str.split('|') if alt.strip()]
                else:
                    isoform_defs = [isoform_str]
                
                # Process each isoform definition separately
                for idx, iso_def in enumerate(isoform_defs, start=1):
                    exons = parse_isoform(iso_def)
                    if not exons:
                        continue
                    try:
                        seq = extract_sequence(chrom, exons, genome_dict, strand)
                        # Append isoform index if multiple isoforms
                        if len(isoform_defs) > 1:
                            record_id = f"{circ_id}|{gene_name}|{circ_type}_iso{idx}"
                        else:
                            record_id = f"{circ_id}|{gene_name}|{circ_type}"
                        record_description = f"isoform={iso_def} strand={strand} circ_type={circ_type}"
                        record = SeqRecord(Seq(seq), id=record_id, description=record_description)
                        output_records.append(record)
                    except KeyError:
                        print(f"Warning: Chromosome {chrom} not found in reference genome")
                        continue
    
    # Write sequences to FASTA file
    SeqIO.write(output_records, args.output, "fasta")
    print(f"Extracted {len(output_records)} circRNA sequences to {args.output}")

if __name__ == "__main__":
    main() 