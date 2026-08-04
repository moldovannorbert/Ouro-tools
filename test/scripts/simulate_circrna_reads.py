#!/usr/bin/env python3

import argparse
import gzip
import random
import os
import sys
import multiprocessing
from dataclasses import dataclass
from typing import Tuple, List, Optional
import pysam
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor


@dataclass
class Adapter:
    """Class to store adapter information"""
    name: str
    start_sequence: Tuple[str, str] = None
    end_sequence: Tuple[str, str] = None


@dataclass
class Barcode:
    """Class to store barcode information"""
    name: str
    start_sequence: Tuple[str, str] = None
    end_sequence: Tuple[str, str] = None


# Define adapters
ADAPTERS = [
    Adapter('SQK-NSK007',
                    start_sequence=('SQK-NSK007_Y_Top', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),

            Adapter('Rapid',
                    start_sequence=('Rapid_adapter',
                                    'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA')),

            Adapter('RBK004_upstream',
                    start_sequence=('RBK004_upstream', 'AATGTACTTCGTTCAGTTACGGCTTGGGTGTTTAACC')),

            Adapter('SQK-MAP006',
            start_sequence=('SQK-MAP006_Y_Top_SK63', 'GGTTGTTTCTGTTGGTGCTGATATTGCT'),
            end_sequence=('SQK-MAP006_Y_Bottom_SK64', 'GCAATATCAGCACCAACAGAAA')),

            Adapter('SQK-MAP006 short',
            start_sequence=('SQK-MAP006_Short_Y_Top_LI32', 'CGGCGTCTGCTTGGGTGTTTAACCT'),
            end_sequence=('SQK-MAP006_Short_Y_Bottom_LI33', 'GGTTAAACACCCAAGCAGACGCCG')),

            # The PCR adapters are used both in PCR DNA kits and some cDNA kits.
            Adapter('PCR adapters 1',
                    start_sequence=('PCR_1_start', 'ACTTGCCTGTCGCTCTATCTTC'),
            end_sequence=('PCR_1_end', 'GAAGATAGAGCGACAGGCAAGT')),

            Adapter('PCR adapters 2',
                    start_sequence=('PCR_2_start', 'TTTCTGTTGGTGCTGATATTGC'),
            end_sequence=('PCR_2_end', 'GCAATATCAGCACCAACAGAAA')),

            Adapter('PCR adapters 3',
                    start_sequence=('PCR_3_start', 'TACTTGCCTGTCGCTCTATCTTC'),
            end_sequence=('PCR_3_end', 'GAAGATAGAGCGACAGGCAAGTA')),

    # 1D^2 kit adapters
            Adapter('1D^2 part 1',
                    start_sequence=('1D2_part_1_start', 'GAGAGGTTCCAAGTCAGAGAGGTTCCT'),
            end_sequence=('1D2_part_1_end', 'AGGAACCTCTCTGACTTGGAACCTCTC')),

            Adapter('1D^2 part 2',
                    start_sequence=('1D2_part_2_start', 'CTTCGTTCAGTTACGTATTGCTGGCGTCTGCTT'),
            end_sequence=('1D2_part_2_end', 'CACCCAAGCAGACGCCAGCAATACGTAACT')),

            Adapter('cDNA SSP',
            start_sequence=('cDNA_SSP', 'TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG'),
            end_sequence=('cDNA_SSP_rev', 'CCCGGCCGTAATGGCAGCAATATCAGCACCAACAGAAA')),

    # Custom adapters for the Ouro-seq protocol.
    Adapter('Ouro-seq',
            start_sequence=('P-SSP', 'AAGCAGTGGTATCAACGCAGAGTGGTTT'),
            end_sequence=('RH-P-CRTA', 'AAGCAGTGGTATCAACGCAGAGTATGCAACGCAACT')),
]


# Define barcodes
BARCODES = [
                # Some barcoding kits (like the native barcodes) use the rev comp barcode at the start
            # of the read and the forward barcode at the end of the read.
    Barcode('Barcode 1 (reverse)',
                    start_sequence=('BC01_rev', 'CACAAAGACACCGACAACTTTCTT'),
                    end_sequence=('BC01', 'AAGAAAGTTGTCGGTGTCTTTGTG')),
    Barcode('Barcode 2 (reverse)',
                    start_sequence=('BC02_rev', 'ACAGACGACTACAAACGGAATCGA'),
                    end_sequence=('BC02', 'TCGATTCCGTTTGTAGTCGTCTGT')),
    Barcode('Barcode 3 (reverse)',
                    start_sequence=('BC03_rev', 'CCTGGTAACTGGGACACAAGACTC'),
                    end_sequence=('BC03', 'GAGTCTTGTGTCCCAGTTACCAGG')),
    Barcode('Barcode 4 (reverse)',
                    start_sequence=('BC04_rev', 'TAGGGAAACACGATAGAATCCGAA'),
                    end_sequence=('BC04', 'TTCGGATTCTATCGTGTTTCCCTA')),
    Barcode('Barcode 5 (reverse)',
                    start_sequence=('BC05_rev', 'AAGGTTACACAAACCCTGGACAAG'),
                    end_sequence=('BC05', 'CTTGTCCAGGGTTTGTGTAACCTT')),
    Barcode('Barcode 6 (reverse)',
                    start_sequence=('BC06_rev', 'GACTACTTTCTGCCTTTGCGAGAA'),
                    end_sequence=('BC06', 'TTCTCGCAAAGGCAGAAAGTAGTC')),
    Barcode('Barcode 7 (reverse)',
                    start_sequence=('BC07_rev', 'AAGGATTCATTCCCACGGTAACAC'),
                    end_sequence=('BC07', 'GTGTTACCGTGGGAATGAATCCTT')),
    Barcode('Barcode 8 (reverse)',
                    start_sequence=('BC08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA'),
                    end_sequence=('BC08', 'TTCAGGGAACAAACCAAGTTACGT')),
    Barcode('Barcode 9 (reverse)',
                    start_sequence=('BC09_rev', 'AACCAAGACTCGCTGTGCCTAGTT'),
                    end_sequence=('BC09', 'AACTAGGCACAGCGAGTCTTGGTT')),
    Barcode('Barcode 10 (reverse)',
                    start_sequence=('BC10_rev', 'GAGAGGACAAAGGTTTCAACGCTT'),
                    end_sequence=('BC10', 'AAGCGTTGAAACCTTTGTCCTCTC')),
    Barcode('Barcode 11 (reverse)',
                    start_sequence=('BC11_rev', 'TCCATTCCCTCCGATAGATGAAAC'),
                    end_sequence=('BC11', 'GTTTCATCTATCGGAGGGAATGGA')),
    Barcode('Barcode 12 (reverse)',
                    start_sequence=('BC12_rev', 'TCCGATTCTGCTTCTTTCTACCTG'),
                    end_sequence=('BC12', 'CAGGTAGAAAGAAGCAGAATCGGA')),

            # Other barcoding kits (like the PCR and rapid barcodes) use the forward barcode at the
            # start of the read and the rev comp barcode at the end of the read.
    Barcode('Barcode 1 (forward)',
                    start_sequence=('BC01', 'AAGAAAGTTGTCGGTGTCTTTGTG'),
                    end_sequence=('BC01_rev', 'CACAAAGACACCGACAACTTTCTT')),
    Barcode('Barcode 2 (forward)',
                    start_sequence=('BC02', 'TCGATTCCGTTTGTAGTCGTCTGT'),
                    end_sequence=('BC02_rev', 'ACAGACGACTACAAACGGAATCGA')),
    Barcode('Barcode 3 (forward)',
                    start_sequence=('BC03', 'GAGTCTTGTGTCCCAGTTACCAGG'),
                    end_sequence=('BC03_rev', 'CCTGGTAACTGGGACACAAGACTC')),
    Barcode('Barcode 4 (forward)',
                    start_sequence=('BC04', 'TTCGGATTCTATCGTGTTTCCCTA'),
                    end_sequence=('BC04_rev', 'TAGGGAAACACGATAGAATCCGAA')),
    Barcode('Barcode 5 (forward)',
                    start_sequence=('BC05', 'CTTGTCCAGGGTTTGTGTAACCTT'),
                    end_sequence=('BC05_rev', 'AAGGTTACACAAACCCTGGACAAG')),
    Barcode('Barcode 6 (forward)',
                    start_sequence=('BC06', 'TTCTCGCAAAGGCAGAAAGTAGTC'),
                    end_sequence=('BC06_rev', 'GACTACTTTCTGCCTTTGCGAGAA')),
    Barcode('Barcode 7 (forward)',
                    start_sequence=('BC07', 'GTGTTACCGTGGGAATGAATCCTT'),
                    end_sequence=('BC07_rev', 'AAGGATTCATTCCCACGGTAACAC')),
    Barcode('Barcode 8 (forward)',
                    start_sequence=('BC08', 'TTCAGGGAACAAACCAAGTTACGT'),
                    end_sequence=('BC08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA')),
    Barcode('Barcode 9 (forward)',
                    start_sequence=('BC09', 'AACTAGGCACAGCGAGTCTTGGTT'),
                    end_sequence=('BC09_rev', 'AACCAAGACTCGCTGTGCCTAGTT')),
    Barcode('Barcode 10 (forward)',
                    start_sequence=('BC10', 'AAGCGTTGAAACCTTTGTCCTCTC'),
                    end_sequence=('BC10_rev', 'GAGAGGACAAAGGTTTCAACGCTT')),
    Barcode('Barcode 11 (forward)',
                    start_sequence=('BC11', 'GTTTCATCTATCGGAGGGAATGGA'),
                    end_sequence=('BC11_rev', 'TCCATTCCCTCCGATAGATGAAAC')),
    Barcode('Barcode 12 (forward)',
                    start_sequence=('BC12', 'CAGGTAGAAAGAAGCAGAATCGGA'),
            end_sequence=('BC12_rev', 'TCCGATTCTGCTTCTTTCTACCTG'))
]


def make_full_native_barcode_adapter(barcode_num):
    """Create a full native barcode adapter with the complete sequence."""
    barcode = next((x for x in BARCODES if x.name == f'Barcode {barcode_num} (reverse)'), None)
    if not barcode:
        return None
    
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter(f'Native barcoding {barcode_num} (full sequence)',
                   start_sequence=(f'NB{barcode_num:02d}_start', start_full_seq),
                   end_sequence=(f'NB{barcode_num:02d}_end', end_full_seq))


def make_old_full_rapid_barcode_adapter(barcode_num):
    """Create a full rapid barcode adapter (SQK-RBK001) with the complete sequence."""
    barcode = next((x for x in BARCODES if x.name == f'Barcode {barcode_num} (forward)'), None)
    if not barcode:
        return None
    
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'TATTGCT' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter(f'Rapid barcoding {barcode_num} (full sequence, old)',
                   start_sequence=(f'RB{barcode_num:02d}_full', start_full_seq))


def make_new_full_rapid_barcode_adapter(barcode_num):
    """Create a full rapid barcode adapter (SQK-RBK004) with the complete sequence."""
    barcode = next((x for x in BARCODES if x.name == f'Barcode {barcode_num} (forward)'), None)
    if not barcode:
        return None
    
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'GCTTGGGTGTTTAACC' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter(f'Rapid barcoding {barcode_num} (full sequence, new)',
                   start_sequence=(f'RB{barcode_num:02d}_full', start_full_seq))


# Create full adapter sequences for all barcodes
FULL_ADAPTERS = []

# Add native barcoding full adapters (for barcodes 1-12)
for i in range(1, 13):
    adapter = make_full_native_barcode_adapter(i)
    if adapter:
        FULL_ADAPTERS.append(adapter)

# Add old rapid barcoding full adapters (for barcodes 1-12)
for i in range(1, 13):
    adapter = make_old_full_rapid_barcode_adapter(i)
    if adapter:
        FULL_ADAPTERS.append(adapter)

# Add new rapid barcoding full adapters (for barcodes 1-12)
for i in range(1, 13):
    adapter = make_new_full_rapid_barcode_adapter(i)
    if adapter:
        FULL_ADAPTERS.append(adapter)

# Add full adapters to the main adapter list
ADAPTERS.extend(FULL_ADAPTERS)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Simulate circular RNA reads for Nanopore sequencing.')
    
    # Add option to list all adapters and barcodes
    parser.add_argument('--list_options', action='store_true',
                        help='List all available adapters and barcodes and exit.')
    
    parser.add_argument('--genome', required=False, help='Path to the reference genome (FASTA format).')
    parser.add_argument('--gtf', required=False, help='Path to the genome annotation (GTF format).')
    parser.add_argument('--output', required=False, help='Output FASTQ.gz file path.')
    parser.add_argument('--num_circrnas', type=int, default=100, help='Number of circRNAs to simulate.')
    parser.add_argument('--max_exons', type=int, default=5, help='Maximum number of exons per circRNA.')
    parser.add_argument('--max_copies', type=int, default=5, help='Maximum number of copies (concatemers) per read.')
    parser.add_argument('--adapter', choices=[a.name for a in ADAPTERS], default='SQK-NSK007', 
                        help='Adapter to use for simulation.')
    parser.add_argument('--barcodes', nargs='+', choices=[b.name for b in BARCODES],
                        default=['Barcode 1 (forward)'], help='Barcodes to use for simulation.')
    parser.add_argument('--read_length', type=int, default=1000, 
                        help='Average read length (before concatemerization).')
    parser.add_argument('--read_length_std', type=int, default=200, 
                        help='Standard deviation of read length.')
    parser.add_argument('--num_reads', type=int, default=1000, 
                        help='Number of reads to simulate.')
    parser.add_argument('--min_detections', type=int, default=5,
                        help='Minimum number of times each circRNA should be detected/sequenced.')
    parser.add_argument('--threads', type=int, default=1, 
                        help='Number of CPU threads to use.')
    parser.add_argument('--seed', type=int, default=42, 
                        help='Random seed for reproducibility.')
    
    args = parser.parse_args()
    
    # If --list_options is specified, print all adapters and barcodes and exit
    if args.list_options:
        print("Available Adapters:")
        print("-------------------")
        for adapter in ADAPTERS:
            print(f"  {adapter.name}")
            if adapter.start_sequence:
                print(f"    Start: {adapter.start_sequence[0]} - {adapter.start_sequence[1][:20]}...")
            if adapter.end_sequence:
                print(f"    End:   {adapter.end_sequence[0]} - {adapter.end_sequence[1][:20]}...")
        
        print("\nAvailable Barcodes:")
        print("------------------")
        for barcode in BARCODES:
            print(f"  {barcode.name}")
            if barcode.start_sequence:
                print(f"    Start: {barcode.start_sequence[0]} - {barcode.start_sequence[1]}")
            if barcode.end_sequence:
                print(f"    End:   {barcode.end_sequence[0]} - {barcode.end_sequence[1]}")
        
        sys.exit(0)
    
    # Check required arguments if not just listing options
    if not args.genome:
        parser.error("--genome is required unless --list_options is specified")
    if not args.gtf:
        parser.error("--gtf is required unless --list_options is specified")
    if not args.output:
        parser.error("--output is required unless --list_options is specified")
    
    return args


def parse_gtf(gtf_path):
    """Parse GTF file to extract exon information."""
    print(f"Parsing GTF file: {gtf_path}")
    exons = []
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
                
            chrom = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            strand = fields[6]
            
            # Parse attributes
            attr_dict = {}
            attrs = fields[8].split(';')
            for attr in attrs:
                attr = attr.strip()
                if not attr:
                    continue
                key_value = attr.split(' ', 1)
                if len(key_value) != 2:
                    continue
                key, value = key_value
                attr_dict[key] = value.strip('"')
            
            # Get gene_id and transcript_id
            gene_id = attr_dict.get('gene_id', 'unknown')
            transcript_id = attr_dict.get('transcript_id', 'unknown')
            
            exons.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'transcript_id': transcript_id
            })
    
    print(f"Extracted {len(exons)} exons from GTF")
    return exons


def get_exon_sequence(genome, exon):
    """Extract exon sequence from the genome."""
    try:
        # Try to get the sequence directly
        seq = genome[exon['chrom']][exon['start']:exon['end']]
        if exon['strand'] == '-':
            seq = seq.reverse_complement()
        return seq
    except KeyError:
        # Handle case where chromosome name doesn't exist in genome
        # Try common chromosome name variations
        chrom = exon['chrom']
        alt_names = []
        
        # Try with/without 'chr' prefix
        if chrom.startswith('chr'):
            alt_names.append(chrom[3:])  # Remove 'chr' prefix
        else:
            alt_names.append('chr' + chrom)  # Add 'chr' prefix
            
        # Try removing version suffixes (e.g., _KZ208906v1_fix)
        if '_' in chrom:
            base_chrom = chrom.split('_')[0]
            alt_names.append(base_chrom)
            if not base_chrom.startswith('chr'):
                alt_names.append('chr' + base_chrom)
        
        # Try the alternative names
        for alt_name in alt_names:
            if alt_name in genome:
                seq = genome[alt_name][exon['start']:exon['end']]
                if exon['strand'] == '-':
                    seq = seq.reverse_complement()
                return seq
        
        # If we get here, we couldn't find a matching chromosome
        print(f"Warning: Chromosome {exon['chrom']} not found in genome. Skipping exon.")
        return Seq("")


def generate_circrna(genome, exons, max_exons):
    """Generate a circular RNA by selecting random exons."""
    # We'll try multiple genes until we find one with proper splice sites
    max_attempts = 20  # Maximum number of genes to try
    
    for _ in range(max_attempts):
        # Select a random gene
        gene_ids = list(set(exon['gene_id'] for exon in exons))
        gene_id = random.choice(gene_ids)
        
        # Get all exons from this gene
        gene_exons = [e for e in exons if e['gene_id'] == gene_id]
        
        if not gene_exons:
            continue
        
        # Select a random number of exons (at least 1, at most max_exons)
        num_exons = random.randint(1, min(max_exons, len(gene_exons)))
        selected_exons = random.sample(gene_exons, num_exons)
        
        # Sort exons by their position in the transcript
        selected_exons.sort(key=lambda e: e['start'])
        
        # Get the strand (all exons should have the same strand)
        strand = selected_exons[0]['strand']
        
        # Check for canonical splice sites
        if len(selected_exons) > 0:
            first_exon = selected_exons[0]
            last_exon = selected_exons[-1]
            
            # Get chromosomes for both exons
            first_chrom = get_matching_chrom(genome, first_exon['chrom'])
            last_chrom = get_matching_chrom(genome, last_exon['chrom'])
            
            if first_chrom is None or last_chrom is None:
                continue
                
            # Check for GT/AG splice site signals based on strand
            if strand == '+':
                # Get the 2 nucleotides after the last exon (donor site, should be GT)
                if last_exon['end'] + 2 <= len(genome[last_chrom]):
                    donor_site = genome[last_chrom][last_exon['end']:last_exon['end']+2]
                else:
                    continue
                
                # Get the 2 nucleotides before the first exon (acceptor site, should be AG)
                if first_exon['start'] >= 2:
                    acceptor_site = genome[first_chrom][first_exon['start']-2:first_exon['start']]
                else:
                    continue
            else:  # strand == '-'
                # On the - strand, we need to take the reverse complement
                # Get the 2 nucleotides before the last exon (donor site on - strand, should be AC [complement of GT])
                if last_exon['start'] >= 2:
                    donor_site = genome[last_chrom][last_exon['start']-2:last_exon['start']]
                    donor_site = donor_site.reverse_complement()
                else:
                    continue
                
                # Get the 2 nucleotides after the first exon (acceptor site on - strand, should be CT [complement of AG])
                if first_exon['end'] + 2 <= len(genome[first_chrom]):
                    acceptor_site = genome[first_chrom][first_exon['end']:first_exon['end']+2]
                    acceptor_site = acceptor_site.reverse_complement()
                else:
                    continue
            
            # Check if the splice sites match the GT/AG pattern
            if str(donor_site).upper() != "GT" or str(acceptor_site).upper() != "AG":
                # Doesn't have canonical splice sites, try another gene
                continue
        
        # Extract and concatenate exon sequences
        circrna_seq = Seq("")
        for exon in selected_exons:
            exon_seq = get_exon_sequence(genome, exon)
            circrna_seq += exon_seq
        
        # If we couldn't get any sequences, try another gene
        if len(circrna_seq) == 0:
            continue
        
        return circrna_seq, strand
    
    # If we tried max_attempts genes and couldn't find one with canonical splice sites
    print("Warning: Could not find a suitable circRNA with canonical GT/AG splice sites")
    return None, "unknown"


def get_matching_chrom(genome, chrom):
    """Find a matching chromosome name in the genome."""
    if chrom in genome:
        return chrom
    
    alt_names = []
    
    # Try with/without 'chr' prefix
    if chrom.startswith('chr'):
        alt_names.append(chrom[3:])  # Remove 'chr' prefix
    else:
        alt_names.append('chr' + chrom)  # Add 'chr' prefix
        
    # Try removing version suffixes (e.g., _KZ208906v1_fix)
    if '_' in chrom:
        base_chrom = chrom.split('_')[0]
        alt_names.append(base_chrom)
        if not base_chrom.startswith('chr'):
            alt_names.append('chr' + base_chrom)
    
    # Try the alternative names
    for alt_name in alt_names:
        if alt_name in genome:
            return alt_name
    
    return None


def simulate_read(circrna_seq, max_copies, strand, adapter, barcode):
    """Simulate a nanopore read from a circular RNA."""
    if not circrna_seq:
        return None
    
    # Determine number of copies (concatemers)
    num_copies = random.randint(1, max_copies)
    
    # Apply random priming (start from random position)
    start_pos = random.randint(0, len(circrna_seq) - 1)
    read_seq = circrna_seq[start_pos:] + circrna_seq[:start_pos]
    
    # Create concatemeric read
    concatemer_seq = read_seq * num_copies
    
    # Check if we should use a full adapter sequence instead
    barcode_num = None
    if barcode.name.startswith('Barcode '):
        try:
            barcode_num = int(barcode.name.split()[1])
        except (ValueError, IndexError):
            pass
    
    # Use full adapter if available and appropriate
    if barcode_num and adapter.name == 'SQK-NSK007' and '(reverse)' in barcode.name:
        # Use native barcoding full adapter
        full_adapter_name = f'Native barcoding {barcode_num} (full sequence)'
        full_adapter = next((a for a in ADAPTERS if a.name == full_adapter_name), None)
        if full_adapter:
            start_adapter = full_adapter.start_sequence[1] if full_adapter.start_sequence else ""
            end_adapter = full_adapter.end_sequence[1] if full_adapter.end_sequence else ""
            return start_adapter + concatemer_seq + end_adapter
    
    elif barcode_num and adapter.name == 'Rapid' and '(forward)' in barcode.name:
        # Use old rapid barcoding full adapter
        full_adapter_name = f'Rapid barcoding {barcode_num} (full sequence, old)'
        full_adapter = next((a for a in ADAPTERS if a.name == full_adapter_name), None)
        if full_adapter:
            start_adapter = full_adapter.start_sequence[1] if full_adapter.start_sequence else ""
            return start_adapter + concatemer_seq
    
    elif barcode_num and adapter.name == 'RBK004_upstream' and '(forward)' in barcode.name:
        # Use new rapid barcoding full adapter
        full_adapter_name = f'Rapid barcoding {barcode_num} (full sequence, new)'
        full_adapter = next((a for a in ADAPTERS if a.name == full_adapter_name), None)
        if full_adapter:
            start_adapter = full_adapter.start_sequence[1] if full_adapter.start_sequence else ""
            return start_adapter + concatemer_seq
    
    # If no full adapter is used, proceed with the standard approach
    # Add adapter and barcode
    start_adapter = ""
    if adapter.start_sequence:
        start_adapter = adapter.start_sequence[1]
    
    end_adapter = ""
    if adapter.end_sequence:
        end_adapter = adapter.end_sequence[1]
    
    # Get barcode sequences
    start_barcode = ""
    if barcode.start_sequence:
        start_barcode = barcode.start_sequence[1]
    
    end_barcode = ""
    if barcode.end_sequence:
        end_barcode = barcode.end_sequence[1]
    
    # Construct the final sequence with adapters and barcodes
    # 5' adapter + 5' barcode + concatemeric sequence + 3' barcode + 3' adapter
    final_seq = start_adapter + start_barcode + concatemer_seq + end_barcode + end_adapter
    
    # Introduce random errors to simulate real nanopore reads (simplified)
    # In a real simulation, you'd want a more sophisticated error model
    error_rate = 0.05  # 5% error rate
    final_seq_list = list(final_seq)
    for i in range(len(final_seq_list)):
        if random.random() < error_rate:
            error_type = random.choice(['sub', 'ins', 'del'])
            if error_type == 'sub':
                final_seq_list[i] = random.choice('ACGT')
            elif error_type == 'ins':
                final_seq_list[i] = final_seq_list[i] + random.choice('ACGT')
            elif error_type == 'del' and i < len(final_seq_list) - 1:
                final_seq_list[i] = ''
    
    final_seq = ''.join(final_seq_list)
    
    return final_seq


def generate_read_batch(args, genome, exons, read_batch_size, read_id_offset):
    """Generate a batch of reads."""
    adapter = next(a for a in ADAPTERS if a.name == args.adapter)
    barcode_objects = [next(b for b in BARCODES if b.name == bc_name) for bc_name in args.barcodes]
    
    reads = []
    
    # First, determine how many unique circRNAs we need
    num_circrnas = read_batch_size // args.min_detections
    if num_circrnas == 0:
        num_circrnas = 1
    
    # Generate the unique circRNAs
    circrnas = []
    for _ in range(num_circrnas):
        circrna_seq, strand = generate_circrna(genome, exons, args.max_exons)
        if circrna_seq is not None:
            circrnas.append((circrna_seq, strand))
    
    # If we couldn't generate any valid circRNAs, return empty list
    if not circrnas:
        return []
    
    # Ensure each circRNA appears at least min_detections times
    read_id = read_id_offset
    for circrna_seq, strand in circrnas:
        # Generate min_detections reads for this circRNA
        for _ in range(args.min_detections):
            if read_id - read_id_offset >= read_batch_size:
                break
                
            # Select a random barcode
            barcode = random.choice(barcode_objects)
            
            # Simulate read
            read_seq = simulate_read(circrna_seq, args.max_copies, strand, adapter, barcode)
            if read_seq is None:
                continue
            
            # Generate quality scores (simplified)
            qual = "".join(chr(random.randint(33, 73)) for _ in range(len(read_seq)))
            
            # Create record
            record = SeqRecord(
                Seq(read_seq),
                id=f"simulated_read_{read_id}",
                description=f"circRNA_len={len(circrna_seq)} copies={args.max_copies} adapter={args.adapter} barcode={barcode.name}",
                letter_annotations={"phred_quality": [ord(q) - 33 for q in qual]}
            )
            
            reads.append(record)
            read_id += 1
    
    # If we still need more reads, add them by randomly selecting from the generated circRNAs
    while len(reads) < read_batch_size:
        circrna_seq, strand = random.choice(circrnas)
        
        # Select a random barcode
        barcode = random.choice(barcode_objects)
        
        # Simulate read
        read_seq = simulate_read(circrna_seq, args.max_copies, strand, adapter, barcode)
        if read_seq is None:
            continue
        
        # Generate quality scores (simplified)
        qual = "".join(chr(random.randint(33, 73)) for _ in range(len(read_seq)))
        
        # Create record
        record = SeqRecord(
            Seq(read_seq),
            id=f"simulated_read_{read_id}",
            description=f"circRNA_len={len(circrna_seq)} copies={args.max_copies} adapter={args.adapter} barcode={barcode.name}",
            letter_annotations={"phred_quality": [ord(q) - 33 for q in qual]}
        )
        
        reads.append(record)
        read_id += 1
    
    return reads


def main():
    args = parse_arguments()
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    print(f"Loading reference genome: {args.genome}")
    genome = {}
    with open(args.genome, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            genome[record.id] = record.seq
    
    print(f"Loaded genome with {len(genome)} chromosomes/contigs")
    
    # Parse GTF file to get exon information
    exons = parse_gtf(args.gtf)
    
    # Check for potential chromosome name mismatches
    gtf_chroms = set(exon['chrom'] for exon in exons)
    genome_chroms = set(genome.keys())
    
    # Check for common prefixes
    has_chr_prefix_gtf = any(chrom.startswith('chr') for chrom in gtf_chroms)
    has_chr_prefix_genome = any(chrom.startswith('chr') for chrom in genome_chroms)
    
    if has_chr_prefix_gtf != has_chr_prefix_genome:
        print("WARNING: Detected potential chromosome name format mismatch between GTF and genome:")
        print(f"  - GTF chromosomes {'have' if has_chr_prefix_gtf else 'do not have'} 'chr' prefix")
        print(f"  - Genome chromosomes {'have' if has_chr_prefix_genome else 'do not have'} 'chr' prefix")
        print("The script will attempt to handle this automatically, but you may get better results by")
        print("ensuring consistent chromosome naming between your GTF and genome files.")
    
    # Calculate batch size for parallel processing
    reads_per_worker = args.num_reads // args.threads
    
    # Generate reads in parallel
    all_reads = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for i in range(args.threads):
            read_batch_size = reads_per_worker
            if i == args.threads - 1:  # Last worker takes remaining reads
                read_batch_size += args.num_reads % args.threads
            
            read_id_offset = i * reads_per_worker
            future = executor.submit(
                generate_read_batch, 
                args, 
                genome, 
                exons, 
                read_batch_size, 
                read_id_offset
            )
            futures.append(future)
        
        for future in futures:
            try:
                batch_reads = future.result()
                all_reads.extend(batch_reads)
            except Exception as e:
                print(f"ERROR in worker process: {e}")
                raise
    
    print(f"Writing {len(all_reads)} simulated reads to {args.output}")
    
    # Write to FASTQ.gz
    with gzip.open(args.output, 'wt') as f:
        SeqIO.write(all_reads, f, 'fastq')


if __name__ == "__main__":
    main() 