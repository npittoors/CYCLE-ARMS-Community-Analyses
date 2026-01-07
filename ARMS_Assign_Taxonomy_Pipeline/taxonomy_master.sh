#!/bin/bash
# Configurable Taxonomy Pipeline
# Author: Based on Nicole Pittoors' scripts
# Date: April 22, 2025
# 
# This script provides a configurable pipeline for taxonomic assignment
# that can handle different genetic markers (COI or 18S) and database modes
# and can set prioirtizations for each database (e.g. MZG then NCBI)
# (full: NMNH+MZG+NCBI or reduced: NMNH+MZG only)

# Set strict error handling
set -e  # Exit on any error
set -u  # Treat unset variables as errors

# Default configuration
MARKER="COI"                   # COI or 18S
DB_MODE="full"                 # full or reduced
ASV_FASTA=""                   # Required input
OUTPUT_DIR="taxonomy_results"  # Output directory
THREADS=8                      # Number of CPU threads
EMAIL=""                       # Email for NCBI Entrez
BACKFILL="none"                # Taxonomy backfilling (none, worms, gbif)
DEBUG=false                    # Enable debug mode
PRIORITIZE_DB=true             # Whether to prioritize database source over percent identity
DB_ORDER="MZG,NMNH,NCBI"       # Database priority order (comma-separated, highest priority first)
# Marker-specific settings (will be set based on MARKER value)
IDENTITY_THRESHOLD=85          # Minimum percent identity threshold
MIN_ALIGNMENT=200              # Minimum alignment length
REFERENCE_PATHS=()             # Array to store reference paths
SEPARATE_BY_PHYLUM="false"  # Whether to separate metazoans from non-metazoans

# Function to display help
function show_help {
    echo "Configurable Taxonomy Pipeline"
    echo "Usage: $0 [options]"
    echo ""
    echo "Required options:"
    echo "  --asv-fasta FILE         ASV sequences FASTA file"
    echo ""
    echo "Basic options:"
    echo "  --marker TYPE            Marker type: COI or 18S (default: $MARKER)"
    echo "  --db-mode MODE           Database mode: full or reduced (default: $DB_MODE)"
    echo "                           full = NMNH + MZG + NCBI, reduced = NMNH + MZG only"
    echo "  --output-dir DIR         Output directory (default: $OUTPUT_DIR)"
    echo "  --threads NUM            Number of threads (default: $THREADS)"
    echo ""
    echo "Reference database options:"
    echo "  --nmnh-ref FILE          NMNH reference FASTA file (override auto-selection)"
    echo "  --mzg-ref FILE           MetaZooGene reference FASTA file (override auto-selection)"
    echo "  --ncbi-db PATH           NCBI nt database path (override auto-selection)"
    echo ""
    echo "Taxonomy processing options:"
    echo "  --email EMAIL            Email for NCBI Entrez (required for NCBI processing)"
    echo "  --backfill TYPE          Taxonomy backfilling: none, worms, or gbif (default: none)"
    echo "  --prioritize-db            Prioritize database source (MZG > NMNH > NCBI) over percent identity (default: true)"
    echo "  --no-prioritize-db         Sort hits only by percent identity, not by database source"
    echo "  --db-order LIST            Custom database priority order (comma-separated, highest priority first, default: MZG,NMNH,NCBI)"
    echo "  --identity NUM           Minimum percent identity threshold (default: auto by marker)"
    echo "  --alignment NUM          Minimum alignment length (default: auto by marker)"
    echo ""
    echo "  This pipeline automatically generates:"
    echo "  - Taxonomy assignments for sequences meeting criteria"
    echo "  - List of ASVs with no hits meeting criteria"
    echo "  - Complete list of ASVs missing from final taxonomy (missing_asvs.tsv)"
    echo "  - FASTA file with sequences for all missing ASVs (missing_asvs.fasta)"
    echo "Other options:"
    echo "  --debug                  Enable debug mode"
    echo "  --help                   Show this help message"
    exit 1
}

# Function to validate marker type
function validate_marker {
    if [[ ! "$MARKER" =~ ^(COI|18S)$ ]]; then
        echo "Error: Invalid marker type. Must be COI or 18S."
        exit 1
    fi
}

# Function to validate database mode
function validate_db_mode {
    if [[ ! "$DB_MODE" =~ ^(full|reduced)$ ]]; then
        echo "Error: Invalid database mode. Must be 'full' or 'reduced'."
        exit 1
    fi
}

# Function to set marker-specific settings
function set_marker_settings {
    # Set default parameters based on marker
    if [[ "$MARKER" == "COI" ]]; then
        # COI settings
        IDENTITY_THRESHOLD=${IDENTITY_THRESHOLD:-85}
        MIN_ALIGNMENT=${MIN_ALIGNMENT:-200}
        
        # Default reference paths for COI
        NMNH_REF=${NMNH_REF:-"references/COI/CYCLE_COI_refs.fasta"}
        MZG_REF=${MZG_REF:-"references/COI/MZGfasta-coi__T4000000__o00__C_world.dada2.fasta"}
        NCBI_DB=${NCBI_DB:-"/media/other/databases/BLAST/NCBI_nt_v4/nt"}
        
    elif [[ "$MARKER" == "18S" ]]; then
        # 18S settings
        IDENTITY_THRESHOLD=${IDENTITY_THRESHOLD:-85}
        MIN_ALIGNMENT=${MIN_ALIGNMENT:-300}
        
        # Default reference paths for 18S
        NMNH_REF=${NMNH_REF:-"references/18S/CYCLE_18S_refs.fasta"}
        MZG_REF=${MZG_REF:-"references/18S/MZGfasta-18S__T4000000__o00__C_world.dada2.fasta"}
        NCBI_DB=${NCBI_DB:-"/media/other/databases/BLAST/NCBI_nt_v4/nt"}
    fi
    
    # Store them in an array for easier access
    REFERENCE_PATHS=("$NMNH_REF" "$MZG_REF" "$NCBI_DB")
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --marker)
            MARKER="$2"
            shift 2
            ;;
        --db-mode)
            DB_MODE="$2"
            shift 2
            ;;
        --asv-fasta)
            ASV_FASTA="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --nmnh-ref)
            NMNH_REF="$2"
            shift 2
            ;;
        --mzg-ref)
            MZG_REF="$2"
            shift 2
            ;;
        --ncbi-db)
            NCBI_DB="$2"
            shift 2
            ;;
        --email)
            EMAIL="$2"
            shift 2
            ;;
        --backfill)
            BACKFILL="$2"
            shift 2
            ;;
        --identity)
            IDENTITY_THRESHOLD="$2"
            shift 2
            ;;
        --alignment)
            MIN_ALIGNMENT="$2"
            shift 2
            ;;
        --debug)
            DEBUG=true
            shift
            ;;
        --prioritize-db)
            PRIORITIZE_DB=true
            shift
            ;;
        --no-prioritize-db)
            PRIORITIZE_DB=false
            shift
            ;;
        --db-order)
            DB_ORDER="$2"
            shift 2
            ;;
        --separate-by-phylum)
            SEPARATE_BY_PHYLUM="true"
             shift
             ;;

        --help)
            show_help
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Validate inputs
if [ -z "$ASV_FASTA" ]; then
    echo "Error: ASV FASTA file is required (--asv-fasta)"
    show_help
fi

if [ ! -f "$ASV_FASTA" ]; then
    echo "Error: ASV FASTA file not found: $ASV_FASTA"
    exit 1
fi

# Validate marker and database mode
validate_marker
validate_db_mode

# Set marker-specific settings (if not provided)
set_marker_settings

# Check if email is provided for NCBI processing
if [[ "$DB_MODE" == "full" && -z "$EMAIL" ]]; then
    echo "Error: Email is required for NCBI processing in full database mode (--email)"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Set paths for output files
LOGS_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOGS_DIR"

NMNH_BLAST_OUT="${OUTPUT_DIR}/NMNH_blast_results.txt"
NMNH_TAXONOMY="${OUTPUT_DIR}/nmnh_taxonomy_results.tsv"
MZG_TAXONOMY="${OUTPUT_DIR}/MZG_taxonomy_results.tsv"
NCBI_FILTERED="${OUTPUT_DIR}/ncbi_filtered_results.tsv"
COMBINED_TAXONOMY="${OUTPUT_DIR}/combined_taxonomy_results.tsv"
FINAL_TAXONOMY="${OUTPUT_DIR}/final_taxonomy_${MARKER}.tsv"
NO_HITS_FILE="${OUTPUT_DIR}/no_hits.tsv"

# Debug information
if [ "$DEBUG" = true ]; then
    echo "=== Debug Information ==="
    echo "Marker type: $MARKER"
    echo "Database mode: $DB_MODE"
    echo "ASV FASTA: $ASV_FASTA"
    echo "Output directory: $OUTPUT_DIR"
    echo "NMNH reference: $NMNH_REF"
    echo "MZG reference: $MZG_REF"
    if [[ "$DB_MODE" == "full" ]]; then
        echo "NCBI database: $NCBI_DB"
    fi
    echo "Identity threshold: $IDENTITY_THRESHOLD"
    echo "Min alignment length: $MIN_ALIGNMENT"
    echo "Threads: $THREADS"
    echo "Email: $EMAIL"
    echo "Backfill mode: $BACKFILL"
    echo "========================="
    echo "Database prioritization: $PRIORITIZE_DB"
    echo "Database priority order: $DB_ORDER"
fi

# Check for required Python packages
echo "Checking Python requirements..."
python -c "import Bio; import pandas; print('Basic packages imported successfully')"

if [[ "$BACKFILL" == "gbif" ]]; then
    python -c "
try:
    import pygbif
    print('pygbif package imported successfully')
except ImportError:
    print('Warning: pygbif package not found. Please install with pip install pygbif')
"
fi

# Function to run NMNH BLAST and taxonomy processing
function run_nmnh_processing {
    echo "===== Processing NMNH voucher sequences ($MARKER) ====="
    
    # Run BLAST against NMNH voucher sequences
    echo "Running BLAST against NMNH voucher sequences..."
    
    # Create BLAST database if needed
    if [ ! -f "${NMNH_REF}.nin" ]; then
        echo "Creating BLAST database for NMNH references..."
        makeblastdb -in "$NMNH_REF" -dbtype nucl
    fi
    
    # Run BLAST with marker-specific parameters
    blastn -query "$ASV_FASTA" \
        -db "$NMNH_REF" \
        -out "$NMNH_BLAST_OUT" \
        -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
        -perc_identity "$IDENTITY_THRESHOLD" \
        -max_target_seqs 5 \
        -num_threads "$THREADS"
    
    echo "BLAST complete. Results saved to: $NMNH_BLAST_OUT"
    
    # Process NMNH BLAST results to extract taxonomy
    echo "Processing NMNH BLAST results for taxonomy..."
    
    if [ -f "nmnh_taxonomy_processor.py" ]; then
        python nmnh_taxonomy_processor.py \
            --blast "$NMNH_BLAST_OUT" \
            --output "$NMNH_TAXONOMY" \
            --marker "$MARKER" \
            --min-alignment "$MIN_ALIGNMENT"
        
        echo "NMNH taxonomy processing complete. Results saved to: $NMNH_TAXONOMY"
    else
        echo "Error: nmnh_taxonomy_processor.py not found. Cannot continue."
        exit 1
    fi
}

# Function to run MetaZooGene BLAST
function run_mzg_processing {
    echo "===== Processing MetaZooGene database ($MARKER) ====="
    
    # Run BLAST against MetaZooGene database
    echo "Running BLAST against MetaZooGene database..."
    
    # Use the MZG_BLAST.py script with marker-specific parameters
    python MZG_BLAST.py \
        --query "$ASV_FASTA" \
        --reference "$MZG_REF" \
        --output "${OUTPUT_DIR}/MZG_taxonomy" \
        --threads "$THREADS" \
        --marker "$MARKER" \
        --min-identity "$IDENTITY_THRESHOLD" \
        --min-alignment "$MIN_ALIGNMENT"
    
    # Rename the output file to be consistent with our naming
    mv "${OUTPUT_DIR}/MZG_taxonomy.tsv" "$MZG_TAXONOMY"
    
    echo "MetaZooGene processing complete. Results saved to: $MZG_TAXONOMY"
}

# Function to run NCBI BLAST (only in full mode)
function run_ncbi_processing {
    if [[ "$DB_MODE" != "full" ]]; then
        echo "Skipping NCBI processing (reduced database mode selected)"
        return
    fi
    
    echo "===== Processing NCBI nt database ($MARKER) ====="
    
    # Create directory for NCBI results
    NCBI_OUT_DIR="${OUTPUT_DIR}/ncbi_results"
    mkdir -p "$NCBI_OUT_DIR"
    
    # Run BLAST against NCBI nt database
    echo "Running BLAST against NCBI nt database..."
    
    # Check if running on server or local machine
    if [ -f "ncbi_blast.sh" ]; then
        # Use existing script
        echo "Using ncbi_blast.sh script..."
        bash ncbi_blast.sh "$ASV_FASTA" "$NCBI_OUT_DIR" "$NCBI_DB" "$IDENTITY_THRESHOLD" "$MIN_ALIGNMENT" "$THREADS"
        
        # Link output file
        if [ -f "${NCBI_OUT_DIR}/ncbi_filtered_results.tsv" ]; then
            ln -sf "$(readlink -f ${NCBI_OUT_DIR}/ncbi_filtered_results.tsv)" "$(readlink -f $NCBI_FILTERED)"
        else
            echo "Warning: NCBI filtered results not found. Using default BLAST approach."
            
            # Run BLAST with marker-specific parameters
            blastn -query "$ASV_FASTA" \
                -db "$NCBI_DB" \
                -out "${NCBI_OUT_DIR}/ncbi_blast_results.tsv" \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
                -max_target_seqs 200 \
                -num_threads "$THREADS" \
                -perc_identity "$IDENTITY_THRESHOLD"
            
            # Filter results based on criteria
            awk -v min_pident="$IDENTITY_THRESHOLD" -v min_align="$MIN_ALIGNMENT" \
                '$3 >= min_pident && $13 >= 85 && $4 >= min_align' \
                "${NCBI_OUT_DIR}/ncbi_blast_results.tsv" > "$NCBI_FILTERED"
            
            echo "NCBI BLAST complete. Filtered results saved to: $NCBI_FILTERED"
        fi
    else
        # Run BLAST with marker-specific parameters
        blastn -query "$ASV_FASTA" \
            -db "$NCBI_DB" \
            -out "${NCBI_OUT_DIR}/ncbi_blast_results.tsv" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -max_target_seqs 200 \
            -num_threads "$THREADS" \
            -perc_identity "$IDENTITY_THRESHOLD"
        
        # Filter results based on criteria
        awk -v min_pident="$IDENTITY_THRESHOLD" -v min_align="$MIN_ALIGNMENT" \
            '$3 >= min_pident && $13 >= 85 && $4 >= min_align' \
            "${NCBI_OUT_DIR}/ncbi_blast_results.tsv" > "$NCBI_FILTERED"
        
        echo "NCBI BLAST complete. Filtered results saved to: $NCBI_FILTERED"
    fi
}

# Function to select best hits and create combined taxonomy
function run_best_hit_selector {
    echo "===== Selecting best taxonomy hits ====="
    
    # Gather all taxonomy files
    TAXONOMY_FILES=()
    
    if [ -f "$NMNH_TAXONOMY" ]; then
        TAXONOMY_FILES+=("--nmnh" "$NMNH_TAXONOMY")
    else
        echo "Warning: NMNH taxonomy file not found. Continuing without it."
    fi
    
    if [ -f "$MZG_TAXONOMY" ]; then
        TAXONOMY_FILES+=("--mzg" "$MZG_TAXONOMY")
    else
        echo "Warning: MetaZooGene taxonomy file not found. Continuing without it."
    fi
    
    if [[ "$DB_MODE" == "full" && -f "$NCBI_FILTERED" ]]; then
        TAXONOMY_FILES+=("--ncbi" "$NCBI_FILTERED")
    elif [[ "$DB_MODE" == "full" ]]; then
        echo "Warning: NCBI filtered results not found. Continuing without it."
    fi
    
    # Exit if no taxonomy files found
    if [ ${#TAXONOMY_FILES[@]} -eq 0 ]; then
        echo "Error: No taxonomy files found. Exiting."
        exit 1
    fi
    
    # Run best taxonomy selector
    echo "Running best taxonomy selector..."
    
    if [ -f "best_taxonomy_selector.py" ]; then
        # Build command as an array
        CMD=("python" "best_taxonomy_selector.py")
        
        # Add all taxonomy files
        for arg in "${TAXONOMY_FILES[@]}"; do
            CMD+=("$arg")
        done
        
        # Add other required parameters
        CMD+=("--output" "$COMBINED_TAXONOMY")
        CMD+=("--email" "$EMAIL")
        CMD+=("--marker" "$MARKER")
        CMD+=("--db-mode" "$DB_MODE")
        CMD+=("--min-identity" "$IDENTITY_THRESHOLD")
        CMD+=("--min-alignment" "$MIN_ALIGNMENT")
        
        # Add database prioritization options
        if [[ "$PRIORITIZE_DB" == "true" ]]; then
            CMD+=("--prioritize-db")
        else
            CMD+=("--no-prioritize-db")
        fi
        
        # Add database priority order
        if [[ -n "$DB_ORDER" ]]; then
            CMD+=("--db-order" "$DB_ORDER")
        fi
        
        # Execute the command
        echo "Running: ${CMD[*]}"
        "${CMD[@]}"
        
        echo "Best taxonomy selection complete. Results saved to: $COMBINED_TAXONOMY"
    else
        echo "Error: best_taxonomy_selector.py not found. Cannot continue."
        exit 1
    fi
}

# Function to process final taxonomy
function run_taxonomy_processing {
    echo "===== Processing final taxonomy ====="
    
    if [ -f "final_process_taxonomy.py" ]; then
        echo "Running final taxonomy processing..."
        
        # Build command based on options
        CMD=("python" "final_process_taxonomy.py" "$COMBINED_TAXONOMY" "$FINAL_TAXONOMY")
        CMD+=("--marker" "$MARKER")
        CMD+=("--min-pident" "$IDENTITY_THRESHOLD")
        CMD+=("--min-alignment" "$MIN_ALIGNMENT")
        
        # Add backfill option if specified
        if [[ "$BACKFILL" != "none" ]]; then
            CMD+=("--source" "$BACKFILL")
        else
            CMD+=("--no-backfill")
        fi
        
        # Add no-hits output option
        CMD+=("--no-hits-file" "$NO_HITS_FILE")
        
        # Add database filter if specified
        if [[ -n "${DATABASE_FILTER:-}" ]]; then
            CMD+=("--database-filter" "$DATABASE_FILTER")
        fi
        
        # Add the separate-by-phylum option
        if [[ "$SEPARATE_BY_PHYLUM" == "true" ]]; then
            CMD+=("--separate-by-phylum")
        fi
      
        # Add changes log path
        CMD+=("--changes-log" "${OUTPUT_DIR}/taxonomy_changes.log")
        
        # Execute the command
        echo "Running: ${CMD[*]}"
        "${CMD[@]}"
        
        echo "Final taxonomy processing complete. Results saved to: $FINAL_TAXONOMY"
        echo "ASVs with no acceptable hits saved to: $NO_HITS_FILE"
    else
        echo "Warning: final_process_taxonomy.py not found. Skipping final processing."
        cp "$COMBINED_TAXONOMY" "$FINAL_TAXONOMY"
    fi
}

function find_missing_asvs {
    echo "===== Finding ASVs without taxonomy hits ====="
    
    # Create a temporary Python script
    TEMP_SCRIPT="${OUTPUT_DIR}/find_missing_asvs.py"
    
    cat > "$TEMP_SCRIPT" << 'EOF'
#!/usr/bin/env python3
from Bio import SeqIO
import sys
import pandas as pd

fasta_file = sys.argv[1]  # Original ASV FASTA file
results_file = sys.argv[2]  # Final taxonomy results
output_file = sys.argv[3]  # Output file for ASVs with no hits

# Get all ASV IDs from FASTA
all_asvs = set()
seq_records = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    all_asvs.add(record.id)
    seq_records[record.id] = record

print(f"Found {len(all_asvs)} total ASVs in FASTA file")

# Get ASVs with hits from results file
results_df = pd.read_csv(results_file, sep='\t')
hit_asvs = set(results_df['ASV_ID'])

print(f"Found {len(hit_asvs)} ASVs with taxonomy hits")

# Find ASVs with no hits
missing_asvs = all_asvs - hit_asvs
print(f"Found {len(missing_asvs)} ASVs with no taxonomy hits")

# Write the no-hit ASVs to a file
with open(output_file, 'w') as f:
    f.write("ASV_ID\n")
    for asv_id in sorted(missing_asvs):
        f.write(f"{asv_id}\n")

print(f"Saved {len(missing_asvs)} ASVs with no hits to {output_file}")

# Extract sequences for no-hit ASVs
fasta_output = output_file.replace('.tsv', '.fasta')
with open(fasta_output, 'w') as f:
    for asv_id in sorted(missing_asvs):
        if asv_id in seq_records:
            SeqIO.write(seq_records[asv_id], f, "fasta")

print(f"Saved sequences for {len(missing_asvs)} ASVs with no hits to {fasta_output}")
EOF

    # Make the script executable
    chmod +x "$TEMP_SCRIPT"
    
    # Run the script
    python "$TEMP_SCRIPT" "$ASV_FASTA" "$FINAL_TAXONOMY" "${OUTPUT_DIR}/missing_asvs.tsv"
    
    echo "Missing ASVs analysis complete. Results saved to:"
    echo "  List of ASVs: ${OUTPUT_DIR}/missing_asvs.tsv"
    echo "  Sequences: ${OUTPUT_DIR}/missing_asvs.fasta"
}

# Main execution
echo "===== Starting $MARKER Taxonomy Pipeline ($DB_MODE mode) ====="
echo "Using settings:"
echo "  ASV FASTA: $ASV_FASTA"
echo "  Identity threshold: $IDENTITY_THRESHOLD%"
echo "  Min alignment length: $MIN_ALIGNMENT bp"
echo "  Threads: $THREADS"
echo ""

# Run all steps in sequence
run_nmnh_processing
run_mzg_processing

if [[ "$DB_MODE" == "full" ]]; then
    run_ncbi_processing
fi

run_best_hit_selector
run_taxonomy_processing
find_missing_asvs

echo "===== $MARKER Taxonomy Pipeline Complete ====="
echo "Final results saved to: $FINAL_TAXONOMY"
echo "ASVs with no hits meeting criteria saved to: $NO_HITS_FILE"
echo "Complete list of missing ASVs saved to: ${OUTPUT_DIR}/missing_asvs.tsv"
echo "Sequences for missing ASVs saved to: ${OUTPUT_DIR}/missing_asvs.fasta"