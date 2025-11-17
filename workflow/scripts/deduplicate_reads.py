#!/usr/bin/env python3
"""
Deduplicate circRNA reads based on start and end sequence similarity.
Identifies PCR duplicates from randomly-primed circRNA long reads and merges them
into polished consensus sequences using multiple sequence alignment.
"""

import argparse
import sys
import gzip
import hashlib
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, Union
import multiprocessing as mp
from functools import partial
import tempfile
import os

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from tqdm import tqdm
    import edlib
    import spoa
except ImportError as e:
    print(f"Error: Missing required dependency: {e}", file=sys.stderr)
    print("Please install required packages: biopython, tqdm, edlib, spoa", file=sys.stderr)
    sys.exit(1)


class ReadRecord:
    """Container for read data with end sequences."""
    
    def __init__(self, record: SeqRecord, start_seq: str, end_seq: str, 
                 quality: Optional[str] = None):
        self.record = record
        self.start_seq = start_seq
        self.end_seq = end_seq
        self.quality = quality
        self.original_id = record.id
        self.sequence = str(record.seq)
    
    def __len__(self):
        return len(self.sequence)


def extract_end_sequences(seq: str, end_length: int) -> Tuple[str, str]:
    """Extract start and end sequences from a read."""
    if len(seq) < 2 * end_length:
        # If read is too short, use the whole sequence for both ends
        return seq, seq
    
    start_seq = seq[:end_length]
    end_seq = seq[-end_length:]
    return start_seq, end_seq


def generate_kmer_hash(start_seq: str, end_seq: str, k: int = 10) -> str:
    """Generate a hash key from k-mers of start and end sequences."""
    def get_kmers(sequence: str, k: int) -> set:
        if len(sequence) < k:
            return {sequence}
        return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}
    
    start_kmers = get_kmers(start_seq, k)
    end_kmers = get_kmers(end_seq, k)
    
    # Combine k-mers and create hash
    combined_kmers = sorted(start_kmers | end_kmers)
    kmer_string = "|".join(combined_kmers)
    return hashlib.md5(kmer_string.encode()).hexdigest()


def calculate_similarity(seq1: str, seq2: str) -> float:
    """Calculate similarity between two sequences using edlib alignment."""
    if seq1 == seq2:
        return 1.0
    
    # Use edlib for fast alignment
    result = edlib.align(seq1, seq2, mode="HW", task="distance")
    edit_distance = result["editDistance"]
    max_len = max(len(seq1), len(seq2))
    
    if max_len == 0:
        return 1.0
    
    similarity = 1.0 - (edit_distance / max_len)
    return similarity


def validate_similarity(reads_in_group: List[ReadRecord], 
                       similarity_threshold: float, 
                       end_length: int) -> List[List[ReadRecord]]:
    """Validate similarity within a group and split into true clusters."""
    if len(reads_in_group) <= 1:
        return [reads_in_group]
    
    clusters = []
    remaining_reads = reads_in_group.copy()
    
    while remaining_reads:
        # Start a new cluster with the first remaining read
        current_cluster = [remaining_reads.pop(0)]
        
        # Find all reads similar to any read in current cluster
        i = 0
        while i < len(remaining_reads):
            read = remaining_reads[i]
            is_similar = False
            
            # Check similarity to any read in current cluster
            for cluster_read in current_cluster:
                # Check both start and end similarity
                start_sim = calculate_similarity(read.start_seq, cluster_read.start_seq)
                end_sim = calculate_similarity(read.end_seq, cluster_read.end_seq)
                
                # Both ends must be similar
                if start_sim >= similarity_threshold and end_sim >= similarity_threshold:
                    is_similar = True
                    break
            
            if is_similar:
                current_cluster.append(read)
                remaining_reads.pop(i)
            else:
                i += 1
        
        clusters.append(current_cluster)
    
    return clusters


def generate_consensus(cluster_reads: List[ReadRecord]) -> ReadRecord:
    """Generate consensus sequence from a cluster of reads using spoa."""
    if len(cluster_reads) == 1:
        return cluster_reads[0]
    
    # Extract sequences for alignment
    sequences = [read.sequence for read in cluster_reads]
    
    try:
        # Use spoa for multiple sequence alignment and consensus
        # algorithm=0 is local alignment, algorithm=1 is global alignment
        result = spoa.poa(sequences, algorithm=0)
        consensus_seq = result[0]  # First element is the consensus sequence
        
        # Create new record with consensus
        # Use the longest original read as template for metadata
        template_read = max(cluster_reads, key=len)
        
        # Update ID to indicate this is a consensus
        new_id = f"{template_read.original_id}_consensus_{len(cluster_reads)}dups"
        
        consensus_record = SeqRecord(
            Seq(consensus_seq),
            id=new_id,
            description=f"consensus from {len(cluster_reads)} reads"
        )
        
        # For FASTQ, generate average quality scores
        quality = None
        if template_read.quality:
            # Simple approach: use the quality from the longest read
            # In a more sophisticated implementation, you'd calculate consensus quality
            quality = template_read.quality
        
        return ReadRecord(consensus_record, 
                         extract_end_sequences(consensus_seq, len(template_read.start_seq))[0],
                         extract_end_sequences(consensus_seq, len(template_read.start_seq))[1],
                         quality)
    
    except Exception as e:
        print(f"Warning: spoa consensus failed for cluster, using longest read: {e}", file=sys.stderr)
        # Fallback: return the longest read
        return max(cluster_reads, key=len)


def process_cluster(cluster: List[ReadRecord]) -> Tuple[ReadRecord, Dict[str, int]]:
    """Process a single cluster and return consensus + statistics."""
    consensus = generate_consensus(cluster)
    
    stats = {
        'total_reads': len(cluster),
        'unique_reads': 1 if len(cluster) == 1 else 1,
        'duplicate_reads': len(cluster) - 1 if len(cluster) > 1 else 0
    }
    
    return consensus, stats


def group_by_kmers(reads: List[ReadRecord], end_length: int, k: int = 10) -> Dict[str, List[ReadRecord]]:
    """Group reads by k-mer hash of their end sequences."""
    groups = defaultdict(list)
    
    for read in reads:
        hash_key = generate_kmer_hash(read.start_seq, read.end_seq, k)
        groups[hash_key].append(read)
    
    return dict(groups)


def load_reads(input_file: str, end_length: int) -> List[ReadRecord]:
    """Load reads from FASTQ or FASTA file and extract end sequences."""
    reads = []
    
    # Determine file format and compression
    file_path = Path(input_file)
    is_gzipped = file_path.suffix == '.gz'
    
    # Open file with appropriate handler
    if is_gzipped:
        file_handle = gzip.open(input_file, 'rt')
    else:
        file_handle = open(input_file, 'r')
    
    try:
        # Determine format from file extension
        if file_path.suffixes[-1] in ['.fastq', '.fq']:
            format_type = 'fastq'
        elif file_path.suffixes[-1] in ['.fasta', '.fa']:
            format_type = 'fasta'
        else:
            # Try to auto-detect from content
            first_line = file_handle.readline().strip()
            file_handle.seek(0)
            if first_line.startswith('@'):
                format_type = 'fastq'
            elif first_line.startswith('>'):
                format_type = 'fasta'
            else:
                raise ValueError(f"Cannot determine file format for {input_file}")
        
        # Parse sequences
        for record in SeqIO.parse(file_handle, format_type):
            start_seq, end_seq = extract_end_sequences(str(record.seq), end_length)
            
            # Extract quality if available
            quality = None
            if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
                quality = ''.join(chr(q + 33) for q in record.letter_annotations['phred_quality'])
            
            read_record = ReadRecord(record, start_seq, end_seq, quality)
            reads.append(read_record)
    
    finally:
        file_handle.close()
    
    return reads


def write_output(consensus_reads: List[ReadRecord], output_file: str, 
                input_format: str, is_gzipped: bool):
    """Write consensus reads to output file."""
    output_path = Path(output_file)
    
    # Determine output format
    if input_format == 'fastq':
        output_format = 'fastq'
    else:
        output_format = 'fasta'
    
    # Open output file
    if is_gzipped:
        file_handle = gzip.open(output_file, 'wt')
    else:
        file_handle = open(output_file, 'w')
    
    try:
        # Convert ReadRecord objects back to SeqRecord for writing
        records = []
        for read_record in consensus_reads:
            record = read_record.record
            # Update sequence if consensus was generated
            if str(record.seq) != read_record.sequence:
                record.seq = Seq(read_record.sequence)
            
            # Add quality scores for FASTQ
            if output_format == 'fastq' and read_record.quality:
                record.letter_annotations = {'phred_quality': [ord(q) - 33 for q in read_record.quality]}
            
            records.append(record)
        
        SeqIO.write(records, file_handle, output_format)
    
    finally:
        file_handle.close()


def write_statistics(all_stats: List[Dict[str, int]], output_file: str):
    """Write deduplication statistics to TSV file."""
    # Aggregate statistics
    total_reads = sum(stats['total_reads'] for stats in all_stats)
    unique_reads = sum(stats['unique_reads'] for stats in all_stats)
    duplicate_reads = sum(stats['duplicate_reads'] for stats in all_stats)
    
    # Create output directory if needed
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("Metric\tCount\n")
        f.write(f"Total_input_reads\t{total_reads}\n")
        f.write(f"Unique_reads\t{unique_reads}\n")
        f.write(f"Duplicate_reads\t{duplicate_reads}\n")
        f.write(f"Deduplication_rate\t{duplicate_reads/total_reads:.4f}\n")


def process_reads_parallel(input_file: str, output_file: str, stats_file: str,
                          similarity_threshold: float, end_length: int, 
                          threads: int) -> None:
    """Main processing function with parallel execution."""
    print(f"Loading reads from {input_file}")
    reads = load_reads(input_file, end_length)
    print(f"Loaded {len(reads)} reads")
    
    if len(reads) == 0:
        print("No reads found in input file")
        # Create empty output files
        write_output([], output_file, 'fasta', False)
        write_statistics([], stats_file)
        return
    
    # Group reads by k-mer hashing
    print("Grouping reads by k-mer similarity...")
    kmer_groups = group_by_kmers(reads, end_length)
    print(f"Created {len(kmer_groups)} k-mer groups")
    
    # Validate similarity within groups
    print("Validating similarity within groups...")
    all_clusters = []
    for group_reads in tqdm(kmer_groups.values(), desc="Validating groups"):
        clusters = validate_similarity(group_reads, similarity_threshold, end_length)
        all_clusters.extend(clusters)
    
    print(f"Created {len(all_clusters)} final clusters")
    
    # Process clusters in parallel
    print("Generating consensus sequences...")
    with mp.Pool(threads) as pool:
        results = list(tqdm(
            pool.imap(process_cluster, all_clusters),
            total=len(all_clusters),
            desc="Processing clusters"
        ))
    
    # Extract consensus reads and statistics
    consensus_reads = [result[0] for result in results]
    all_stats = [result[1] for result in results]
    
    # Write output
    print(f"Writing deduplicated reads to {output_file}")
    input_format = 'fastq' if Path(input_file).suffixes[-1] in ['.fastq', '.fq'] else 'fasta'
    is_gzipped = Path(output_file).suffix == '.gz'
    write_output(consensus_reads, output_file, input_format, is_gzipped)
    
    # Write statistics
    print(f"Writing statistics to {stats_file}")
    write_statistics(all_stats, stats_file)
    
    # Print summary
    total_input = sum(stats['total_reads'] for stats in all_stats)
    total_output = len(consensus_reads)
    duplicates = total_input - total_output
    print(f"Summary: {total_input} input reads -> {total_output} unique reads ({duplicates} duplicates removed)")


def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate circRNA reads based on start/end sequence similarity"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTQ/FASTA file (supports .gz compression)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file (format matches input)"
    )
    parser.add_argument(
        "--stats",
        required=True,
        help="Output TSV file for deduplication statistics"
    )
    parser.add_argument(
        "--similarity",
        type=float,
        default=0.8,
        help="Similarity threshold for duplicate detection (0.0-1.0, default: 0.8)"
    )
    parser.add_argument(
        "--end-length",
        type=int,
        default=20,
        help="Number of bases to check on each end (default: 20)"
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
    
    # Validate input file
    if not Path(args.input).exists():
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Validate parameters
    if not 0.0 <= args.similarity <= 1.0:
        print("Error: Similarity threshold must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)
    
    if args.end_length < 1:
        print("Error: End length must be at least 1", file=sys.stderr)
        sys.exit(1)
    
    if args.threads < 1:
        print("Error: Thread count must be at least 1", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    try:
        process_reads_parallel(
            args.input, args.output, args.stats,
            args.similarity, args.end_length, args.threads
        )
        if args.verbose:
            print("Deduplication completed successfully")
    except Exception as e:
        print(f"Error during deduplication: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
