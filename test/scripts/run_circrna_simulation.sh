#!/bin/bash

# Exit on error
set -e

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH"
    exit 1
fi

# Activate the ouro-tools-testing environment
echo "Activating ouro-tools-testing conda environment..."
eval "$(conda shell.bash hook)"
conda activate ouro-tools-testing || { echo "Error: Failed to activate ouro-tools-testing environment. Make sure it exists."; exit 1; }

# Define parameters
GENOME_PATH="$1"
GTF_PATH="$2"
OUTPUT_PATH="$3"
MAX_EXONS="${4:-3}"  # Default to 3 if not provided
MIN_DETECTIONS="${5:-5}"  # Default to 5 if not provided

# Check if required arguments are provided
if [ -z "$GENOME_PATH" ] || [ -z "$GTF_PATH" ] || [ -z "$OUTPUT_PATH" ]; then
    echo "Usage: $0 <genome_path> <gtf_path> <output_path> [max_exons] [min_detections]"
    echo "Example: $0 /path/to/genome.fa /path/to/annotation.gtf /path/to/output.fastq.gz 5 3"
    exit 1
fi

# Check if the files exist
if [ ! -f "$GENOME_PATH" ]; then
    echo "Error: Genome file not found: $GENOME_PATH"
    exit 1
fi

if [ ! -f "$GTF_PATH" ]; then
    echo "Error: GTF file not found: $GTF_PATH"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_PATH")
mkdir -p "$OUTPUT_DIR"

echo "Note: Make sure your genome FASTA and GTF files use consistent chromosome naming."
echo "      The script will try to handle mismatches, but results may be better with consistent naming."

# Run the simulation
echo "Running circular RNA simulation..."
python "$(dirname "$0")/simulate_circrna_reads.py" \
    --genome "$GENOME_PATH" \
    --gtf "$GTF_PATH" \
    --output "$OUTPUT_PATH" \
    --num_reads 100 \
    --max_exons "$MAX_EXONS" \
    --min_detections "$MIN_DETECTIONS" \
    --threads 8 \
    --adapter "Ouro-seq" \
    --barcodes "Barcode 9 (reverse)"

echo "Simulation completed successfully!"
echo "Output file: $OUTPUT_PATH" 