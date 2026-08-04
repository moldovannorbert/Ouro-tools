# Graph Report - Ouro-tools  (2026-08-03)

## Corpus Check
- 15 files · ~18,369 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 197 nodes · 278 edges · 15 communities (13 shown, 2 thin omitted)
- Extraction: 96% EXTRACTED · 4% INFERRED · 0% AMBIGUOUS · INFERRED: 11 edges (avg confidence: 0.55)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `544bf7da`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- deduplicate_reads.py
- Circr.py
- calculate_potential
- simulate_circrna_reads.py
- detect_fcircrna.py
- merge_ires_predictions.py
- find_ires_like.py
- convert_to_gtf
- main
- Ouro-tools: Pipeline for the Analysis of circRNA from Long-Read Nanopore Sequencing
- process_sequences_parallel
- extract_and_translate_orf
- read_length_distribution.py
- run_circrna_simulation.sh

## God Nodes (most connected - your core abstractions)
1. `calculate_potential()` - 10 edges
2. `ReadRecord` - 10 edges
3. `main()` - 9 edges
4. `process_reads_parallel()` - 9 edges
5. `FindCDS` - 8 edges
6. `Fickett` - 7 edges
7. `third_party_processes()` - 7 edges
8. `process_read()` - 7 edges
9. `main()` - 7 edges
10. `generate_read_batch()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `generate_read_batch()` --calls--> `SeqRecord`  [INFERRED]
  test/scripts/simulate_circrna_reads.py →   _Bridges community 3 → community 0_
- `main()` --calls--> `SeqRecord`  [INFERRED]
  workflow/scripts/extract_circRNA_sequences.py →   _Bridges community 0 → community 8_

## Import Cycles
- None detected.

## Communities (15 total, 2 thin omitted)

### Community 0 - "deduplicate_reads.py"
Cohesion: 0.12
Nodes (26): SeqRecord, calculate_similarity(), extract_end_sequences(), generate_consensus(), generate_kmer_hash(), group_by_kmers(), load_reads(), main() (+18 more)

### Community 1 - "Circr.py"
Cohesion: 0.13
Nodes (25): add_CircBase_annotation(), calculate_exons(), compare_predicted(), create_circular_DB(), examine(), filter_interactions(), fix_bed_scores(), getCircrPath() (+17 more)

### Community 2 - "calculate_potential"
Cohesion: 0.11
Nodes (17): calculate_potential(), Fickett, FindCDS, getstatusoutput(), __main(), mRNA_translate(), protein_param(), calculate Fickett TESTCODE for full sequence 	NAR 10(17) 5303-531 	modified from (+9 more)

### Community 3 - "simulate_circrna_reads.py"
Cohesion: 0.11
Nodes (25): Adapter, Barcode, generate_circrna(), generate_read_batch(), get_exon_sequence(), get_matching_chrom(), main(), make_full_native_barcode_adapter() (+17 more)

### Community 4 - "detect_fcircrna.py"
Cohesion: 0.16
Nodes (17): init_worker(), load_reference(), main(), parse_fasta(), process_read(), Worker initializer to load the mappy aligner and reference genome into global va, Load a reference genome fasta file (or gzipped fasta) and return a dictionary:, Process a single read using the mappy aligner.     If at least two mapping segme (+9 more)

### Community 5 - "merge_ires_predictions.py"
Cohesion: 0.18
Nodes (15): get_ires_sequence(), main(), parse_arguments(), process_circrna(), Extract IRES sequence from duplicated circRNA sequence., Process a single circRNA to find IRES elements upstream of CPC2-predicted ORFs., Parse command line arguments., Read IRES predictions from CSV file. (+7 more)

### Community 6 - "find_ires_like.py"
Cohesion: 0.24
Nodes (11): find_hexamer_matches(), main(), parse_arguments(), parse_hexamers_string(), process_sequence(), Process a single sequence record (for parallel execution)., Parse command line arguments., Read hexamers from a file, one per line. (+3 more)

### Community 7 - "convert_to_gtf"
Cohesion: 0.29
Nodes (9): convert_to_gtf(), format_gtf_attributes(), main(), parse_attributes(), process_isoform(), Process isoform string to get exon coordinates.     Returns list of (start, end), Format attributes for GTF format., Convert CIRI-long .info file to GTF3 format. (+1 more)

### Community 8 - "main"
Cohesion: 0.31
Nodes (8): extract_sequence(), main(), parse_args(), parse_isoform(), Parse command line arguments., Extract exon coordinates from isoform string., Extract sequence from reference genome based on exon coordinates., Main function to extract circRNA sequences.

### Community 9 - "Ouro-tools: Pipeline for the Analysis of circRNA from Long-Read Nanopore Sequencing"
Cohesion: 0.29
Nodes (6): Installation Guide, Ouro-tools: Pipeline for the Analysis of circRNA from Long-Read Nanopore Sequencing, Overview, Repo Contents, System Requirements, User Guide

### Community 10 - "process_sequences_parallel"
Cohesion: 0.47
Nodes (5): duplicate_sequence(), main(), process_sequences_parallel(), Duplicate a single sequence record., Process sequences in parallel.

### Community 11 - "extract_and_translate_orf"
Cohesion: 0.67
Nodes (3): extract_and_translate_orf(), main(), Extract ORF sequence and translate to amino acids based on peptide length.

## Knowledge Gaps
- **6 isolated node(s):** `run_circrna_simulation.sh script`, `Overview`, `Repo Contents`, `System Requirements`, `Installation Guide` (+1 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **2 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `generate_read_batch()` connect `simulate_circrna_reads.py` to `deduplicate_reads.py`?**
  _High betweenness centrality (0.051) - this node is a cross-community bridge._
- **Why does `main()` connect `main` to `deduplicate_reads.py`?**
  _High betweenness centrality (0.024) - this node is a cross-community bridge._
- **What connects `run_circrna_simulation.sh script`, `Overview`, `Repo Contents` to the rest of the system?**
  _6 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `deduplicate_reads.py` be split into smaller, more focused modules?**
  _Cohesion score 0.11576354679802955 - nodes in this community are weakly interconnected._
- **Should `Circr.py` be split into smaller, more focused modules?**
  _Cohesion score 0.1339031339031339 - nodes in this community are weakly interconnected._
- **Should `calculate_potential` be split into smaller, more focused modules?**
  _Cohesion score 0.1111111111111111 - nodes in this community are weakly interconnected._
- **Should `simulate_circrna_reads.py` be split into smaller, more focused modules?**
  _Cohesion score 0.1076923076923077 - nodes in this community are weakly interconnected._