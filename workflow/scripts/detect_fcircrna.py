#!/usr/bin/env python3
import argparse
import gzip
import mappy as mp
import concurrent.futures

# ============================================================================
# Helper functions to load fasta files (reference genome and input reads)
# ============================================================================

def load_reference(ref_filepath):
    """
    Load a reference genome fasta file (or gzipped fasta) and return a dictionary:
         {chromosome_name: sequence_string, ...}
    """
    ref = {}
    open_func = gzip.open if ref_filepath.endswith(".gz") else open
    with open_func(ref_filepath, "rt") as fh:
        current_chrom = None
        seq_lines = []
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_chrom is not None:
                    ref[current_chrom] = "".join(seq_lines)
                # get first token as chrom name
                current_chrom = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_chrom is not None:
            ref[current_chrom] = "".join(seq_lines)
    return ref

def parse_fasta(fasta_path):
    """
    Generator to parse a fasta file (plain text or gzipped) yielding (header, sequence).
    """
    open_func = gzip.open if fasta_path.endswith(".gz") else open
    with open_func(fasta_path, "rt") as fh:
        header = None
        seq_lines = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield (header, "".join(seq_lines))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            yield (header, "".join(seq_lines))

# ============================================================================
# Helper functions for splice signal search within the reference genome region
# ============================================================================

def search_splice_signal_positive_upstream(seq, pos, window, motif):
    """
    Search for motif in the region [pos-window, pos] (upstream of pos)
    Return the (genomic) index of the motif if found, else None.
    We use rfind so that the match closest to pos is returned.
    """
    start = max(0, pos - window)
    region = seq[start:pos]
    i = region.rfind(motif)
    if i != -1:
        return start + i
    return None

def search_splice_signal_positive_downstream(seq, pos, window, motif):
    """
    Search for motif in the region [pos, pos+window] (downstream of pos)
    Return the index of the motif if found, else None.
    """
    region = seq[pos: pos+window]
    i = region.find(motif)
    if i != -1:
        return pos + i
    return None

def search_splice_signal_negative_downstream(seq, pos, window, motif):
    """
    For negative strand candidates: search downstream of pos
    """
    region = seq[pos: pos+window]
    i = region.find(motif)
    if i != -1:
        return pos + i
    return None

def search_splice_signal_negative_upstream(seq, pos, window, motif):
    """
    For negative strand candidates: search upstream of pos.
    Use rfind so that the match closest to pos is returned.
    """
    start = max(0, pos - window)
    region = seq[start: pos]
    i = region.rfind(motif)
    if i != -1:
        return start + i
    return None

# ============================================================================
# Main function
# ============================================================================

def init_worker(ref_file, min_seglen, min_distance, window):
    """
    Worker initializer to load the mappy aligner and reference genome into global variables.
    """
    global ALN, REF_GENOME, MIN_SEGLEN, MIN_DISTANCE, WINDOW
    ALN = mp.Aligner(ref_file)
    if ALN is None:
        raise Exception("Failed to load Aligner with reference " + ref_file)
    REF_GENOME = load_reference(ref_file)
    MIN_SEGLEN = min_seglen
    MIN_DISTANCE = min_distance
    WINDOW = window


def process_read(read_tuple):
    """
    Process a single read using the mappy aligner.
    If at least two mapping segments (of sufficient length) are found, 
    check for candidate fusion circRNA by assessing genomic separation and splice signals.
    Returns a tuple (header_out, seq) if a candidate is found, else None.
    """
    header, seq = read_tuple
    read_id = header.split()[0]
    hits = list(ALN.map(seq, cs=True))
    if len(hits) < 2:
        return None
    hits = sorted(hits, key=lambda h: h.q_st)
    hits = [h for h in hits if (h.q_en - h.q_st) >= MIN_SEGLEN]
    if len(hits) < 2:
        return None
    for i in range(len(hits) - 1):
        h1 = hits[i]
        h2 = hits[i+1]
        # Ensure the query segments do not overlap
        if h2.q_st - h1.q_en < 0:
            continue
        # If on the same reference sequence, check the minimum genomic distance
        if h1.ctg == h2.ctg and abs(h2.r_st - h1.r_en) < MIN_DISTANCE:
            continue

        left_ref = REF_GENOME.get(h1.ctg)
        right_ref = REF_GENOME.get(h2.ctg)
        if left_ref is None or right_ref is None:
            continue

        # Try positive strand configuration: search for "AG" upstream and "GT" downstream.
        pos_left_signal = search_splice_signal_positive_upstream(left_ref, h1.r_en, WINDOW, "AG")
        pos_right_signal = search_splice_signal_positive_downstream(right_ref, h2.r_st, WINDOW, "GT")
        if pos_left_signal is not None and pos_right_signal is not None:
            strand = "+"
            splice_signal = "AG/GT"
            left_signal_dist = h1.r_en - (pos_left_signal + 2)
            right_signal_dist = pos_right_signal - h2.r_st
        else:
            # Try negative strand configuration: search for "CT" downstream and "AC" upstream.
            neg_left_signal = search_splice_signal_negative_downstream(left_ref, h1.r_en, WINDOW, "CT")
            neg_right_signal = search_splice_signal_negative_upstream(right_ref, h2.r_st, WINDOW, "AC")
            if neg_left_signal is not None and neg_right_signal is not None:
                strand = "-"
                splice_signal = "CT/AC"
                left_signal_dist = neg_left_signal - h1.r_en
                right_signal_dist = h2.r_st - (neg_right_signal + 2)
            else:
                strand = "NA"
                splice_signal = "NA"
                left_signal_dist = "NA"
                right_signal_dist = "NA"

        left_mapping_str = f"{h1.ctg}:{h1.r_st}-{h1.r_en}"
        right_mapping_str = f"{h2.ctg}:{h2.r_st}-{h2.r_en}"
        left_exon_str = f"{h1.ctg}:{h1.r_st}-{h1.r_en}|{h1.r_en - h1.r_st}"
        right_exon_str = f"{h2.ctg}:{h2.r_st}-{h2.r_en}|{h2.r_en - h2.r_st}"
        header_out = (f">{read_id} {left_mapping_str}|{right_mapping_str} "
                      f"{strand} {left_exon_str},{right_exon_str} "
                      f"{splice_signal}|{left_signal_dist}-{right_signal_dist}")
        return (header_out, seq)
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Identify fusion circRNA from nanopore long-read sequencing data using approximate mapping."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input fasta file with pre-processed circular consensus reads.")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference genome fasta or gzipped fasta. Used for both mapping and splice signal search.")
    parser.add_argument("-o", "--output_prefix", required=True,
                        help="Output file prefix; results are written to prefix.cand_fcirc.fa")
    parser.add_argument("-d", "--min_distance", type=int, default=1_000_000,
                        help="Minimum genomic distance between fused parts if on the same chromosome (default: 1000000 bp)")
    parser.add_argument("-w", "--window", type=int, default=50,
                        help="Window size for splice signal search around mapping boundaries (default: 50 bp)")
    parser.add_argument("--min_seglen", type=int, default=150,
                        help="Minimum segment length for each mapping (default: 150 bp)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads for parallel processing (default: 1)")
    args = parser.parse_args()

    # Initialize globals once; this ensures only one copy of the reference and aligner is loaded.
    init_worker(args.reference, args.min_seglen, args.min_distance, args.window)
    candidates = []
    # Process reads in parallel using threads to share memory and reduce overhead.
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        # Using executor.map avoids storing all futures in memory and automatically yields results.
        for result in executor.map(process_read, parse_fasta(args.input)):
            if result is not None:
                candidates.append(result)

    out_filename = f"{args.output_prefix}.cand_fcirc.fa"
    with open(out_filename, "w") as out_fa:
        for header_out, seq in candidates:
            out_fa.write(header_out + "\n")
            out_fa.write(seq + "\n")
    print(f"Finished! Fusion circRNA candidates written to {out_filename}")

if __name__ == "__main__":
    main()
