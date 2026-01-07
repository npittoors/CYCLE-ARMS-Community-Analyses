#!/bin/bash

# Script to run BLASTn against the NCBI nt database for 18S sequences
# Usage: ./ncbi_blast_18S.sh <query_fasta> <output_dir> <db_path> [identity_threshold] [min_alignment_length] [threads]

# Parse arguments with defaults
QUERY_FASTA=$1
OUTPUT_DIR=$2
DB_PATH=$3

# Check if parameters are provided, otherwise use defaults
if [ "$#" -ge 4 ]; then
    IDENTITY_THRESHOLD=$4
else
    IDENTITY_THRESHOLD=85
fi
if [ "$#" -ge 5 ]; then
    MIN_ALIGNMENT=$5
else
    MIN_ALIGNMENT=300
fi
if [ "$#" -ge 6 ]; then
    THREADS=$6
else
    THREADS=40
fi

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Print parameters
echo "Running 18S NCBI BLAST with parameters:"
echo "  Query: $QUERY_FASTA"
echo "  Database: $DB_PATH"
echo "  Identity threshold: $IDENTITY_THRESHOLD%"
echo "  Minimum alignment: $MIN_ALIGNMENT bp"
echo "  Threads: $THREADS"

# Run BLASTn
echo "Running BLASTn against NCBI nt database..."
blastn -query $QUERY_FASTA \
    -db $DB_PATH \
    -out $OUTPUT_DIR/ncbi_blast_results.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -max_target_seqs 200 \
    -num_threads $THREADS \
    -perc_identity $IDENTITY_THRESHOLD

# Filter results based on criteria (identity threshold, >85% coverage, min_alignment length)
echo "Filtering results based on criteria..."
awk -v id="$IDENTITY_THRESHOLD" -v len="$MIN_ALIGNMENT" '$3 >= id && $13 >= 85 && $4 >= len' \
    $OUTPUT_DIR/ncbi_blast_results.tsv > $OUTPUT_DIR/ncbi_filtered_results.tsv

echo "BLAST search completed. Results saved to:"
echo "  Raw results: $OUTPUT_DIR/ncbi_blast_results.tsv"
echo "  Filtered results: $OUTPUT_DIR/ncbi_filtered_results.tsv"

# Count results
TOTAL_HITS=$(wc -l < $OUTPUT_DIR/ncbi_blast_results.tsv)
FILTERED_HITS=$(wc -l < $OUTPUT_DIR/ncbi_filtered_results.tsv)

echo "Total hits: $TOTAL_HITS"
echo "Hits passing criteria: $FILTERED_HITS"