# Configurable Taxonomy Pipeline
# Author: Nicole Pittoors (with the assistance claude AI)
# Date Last Modified: April 25, 2025
A flexible pipeline for taxonomic assignment of amplicon sequence variants (ASVs) from marine metabarcoding studies, with support for COI and 18S markers.

## Overview

This pipeline enables taxonomic assignment of ASVs using a combination of reference databases:
- NMNH voucher database (curated reference specimens)
- MetaZooGene database (marine metazoan references)
- NCBI nt database (optional)

The pipeline supports both COI and 18S markers with appropriate threshold settings and can be run in either full mode (using all three databases) or reduced mode (using only NMNH and MetaZooGene databases).

## Features

- Marker-specific settings (COI or 18S)
- Configurable identity and alignment length thresholds
- Full or reduced database modes
- Prioritization of more reliable reference databases over raw percent identity (optional)
- Taxonomic filtering to exclude non-target organisms
- Comprehensive output including no-hit ASVs
- Taxonomy backfilling using WoRMS or GBIF

## Requirements

- Python 3.6+
- BioPython (`pip install biopython`)
- pandas (`pip install pandas`)
- requests (for WoRMS API calls)
- pygbif (optional, for GBIF taxonomy backfilling)
- BLAST+ command-line tools
- Reference databases:
  - NMNH voucher database
  - MetaZooGene reference database (COI/18S)
  - NCBI nt database (for full mode)

## Installation

1. Clone or download all scripts to a directory:
   ```bash
   git clone https://github.com/yourusername/taxonomy-pipeline.git
   cd taxonomy-pipeline
   ```

2. Make the scripts executable:
   ```bash
   chmod +x *.py *.sh
   ```

3. Install required Python packages:
   ```bash
   pip install biopython pandas requests
   
   # For GBIF taxonomy backfilling (optional but recommended):
   pip install pygbif
   ```

4. Ensure BLAST+ is installed and available in your PATH:
   ```bash
   # Check installation
   blastn -version
   
   # If using a conda environment:
   conda activate blast_env
   ```

## Directory Structure

Set up your working directory with the following structure:
```
taxonomy-pipeline/
├── taxonomy_master.sh           # Master script to run the pipeline
├── MZG_BLAST.py                 # MetaZooGene BLAST and taxonomy assignment
├── nmnh_taxonomy_processor.py   # NMNH voucher processor
├── ncbi_blast.sh                # NCBI BLAST script
├── best_taxonomy_selector.py    # Combine and select best hits
├── final_process_taxonomy.py    # Post-process and refine taxonomy
├── README.md                    # This file
├── references/
│   ├── COI/
│   │   ├── NMNH_COI_refs.fasta  # NMNH reference database for COI
│   │   └── MZG_COI_refs.fasta   # MetaZooGene reference for COI
│   └── 18S/
│       ├── NMNH_18S_refs.fasta  # NMNH reference database for 18S
│       └── MZG_18S_refs.fasta   # MetaZooGene reference for 18S
└── data/
    ├── COI_ASVs.fasta           # Your COI ASV sequences
    └── 18S_ASVs.fasta           # Your 18S ASV sequences
```

## Pipeline Components

### 1. Master Pipeline Script (taxonomy_master.sh)

Coordinates the entire workflow:
- Handles marker-specific settings (COI or 18S)
- Supports both full and reduced database modes
- Sets appropriate thresholds for each marker type

### 2. MetaZooGene BLAST (MZG_BLAST.py)

Processes MetaZooGene reference database:
- Handles both COI and 18S marker formats
- Parses taxonomy from MZG headers
- Applies appropriate filtering thresholds

### 3. NMNH Taxonomy Processor (nmnh_taxonomy_processor.py)

Processes NMNH voucher database:
- Parses taxonomic information from BLAST results
- Applies marker-specific thresholds
- Formats output for compatibility with other tools

### 4. Best Taxonomy Selector (best_taxonomy_selector.py)

Combines results from all reference databases:
- Prioritizes hits by database source and/or percent identity
- Handles both full and reduced database modes
- Produces standardized taxonomy output

### 5. Final Taxonomy Processor (final_process_taxonomy.py)

Post-processes combined taxonomy results:
- Filters results based on marker-specific thresholds
- Removes environmental sequences
- Backfills missing taxonomy using WoRMS or GBIF
- Forward-fills taxonomy ranks with appropriate symbols

## Usage

### Running the Master Pipeline

The `taxonomy_master.sh` script manages the entire pipeline. Here are examples for both markers:
Note: Make sure env with pygbif is activated for GBIF backfilling (edna env on Herrera Lab Server)
#### Example 1: COI with Full Database Mode

```bash
./taxonomy_master.sh \
  --marker COI \
  --db-mode full \
  --asv-fasta data/COI_ASVs.fasta \
  --output-dir results_COI \
  --nmnh-ref references/COI/CYCLE_COI_refs.fasta \
  --mzg-ref references/COI/MZGfasta-coi__T4000000__o00__C_world.dada2.fasta \
  --ncbi-db /media/other/databases/BLAST/NCBI_nt_v4/nt \
  --email your@email.com \
  --threads 8 \
  --prioritize-db \
  --backfill gbif \
  --separate-by-phylum
```

#### Example 2: 18S with Reduced Database Mode

```bash
./taxonomy_master.sh \
  --marker 18S \
  --db-mode reduced \
  --asv-fasta data/18S_ASVs.fasta \
  --output-dir results_18S \
  --nmnh-ref references/18S/CYCLE_18S_refs.fasta \
  --mzg-ref references/18S/MZGfasta-18S__T4000000__o00__C_world.dada2.fasta \
  --threads 8 \
  --backfill worms \
  --alignment 300
```

### Command-Line Options

#### Required Parameters
- `--asv-fasta`: Path to your ASV sequences in FASTA format

#### Important Options
- `--marker`: Marker type, either `COI` or `18S` (default: `COI`)
- `--db-mode`: Database mode, either `full` or `reduced` (default: `full`)
- `--output-dir`: Directory for output files (default: `taxonomy_results`)
- `--email`: Your email address (required for NCBI API access in full mode)
- `--threads`: Number of CPU threads to use (default: `8`)

#### Database Parameters
- `--nmnh-ref`: Custom path to NMNH reference FASTA
- `--mzg-ref`: Custom path to MetaZooGene reference FASTA
- `--ncbi-db`: Custom path to NCBI nt database

#### Filtering and Processing Parameters
- `--identity`: Override default percent identity threshold
- `--alignment`: Override default alignment length threshold
- `--prioritize-db`: Prioritize database source over percent identity (default: enabled)
- `--no-prioritize-db`: Disable database prioritization
- `--db-order`: Custom database priority order (comma-separated, default: "MZG,NMNH,NCBI")
- `--backfill`: Taxonomy backfilling using external databases (`none`, `worms`, or `gbif`)
- `--separate-by-phylum`: Separate metazoan and non-metazoan taxa into separate files
- `--debug`: Enable debug mode with verbose output

## Running Individual Scripts

While the master script is recommended, you can also run each component separately:

### 1. MZG BLAST

```bash
python MZG_BLAST.py \
  --query data/COI_ASVs.fasta \
  --reference references/COI/MZG_COI_refs.fasta \
  --output results/MZG_results \
  --marker COI \
  --min-identity 85 \
  --min-alignment 200 \
  --threads 8 \
  --no-hits-file results/MZG_no_hits.tsv
```

### 2. NMNH Taxonomy Processing

```bash
python nmnh_taxonomy_processor.py \
  --blast results/NMNH_blast_results.txt \
  --output results/nmnh_taxonomy_results.tsv \
  --marker COI \
  --min-identity 85 \
  --min-alignment 200 \
  --no-hits-file results/NMNH_no_hits.tsv
```

### 3. NCBI BLAST (Bash Script)

```bash
./ncbi_blast.sh \
  data/COI_ASVs.fasta \
  results/ \
  path/to/ncbi/nt \
  85 \
  200 \
  8
```

### 4. Best Taxonomy Selector

```bash
python best_taxonomy_selector.py \
  --nmnh results/nmnh_taxonomy_results.tsv \
  --mzg results/MZG_results.tsv \
  --ncbi results/ncbi_filtered_results.tsv \
  --output results/combined_taxonomy_results.tsv \
  --marker COI \
  --db-mode full \
  --min-identity 85 \
  --min-alignment 200 \
  --email your@email.com \
  --prioritize-db \
  --db-order "MZG,NMNH,NCBI" \
  --no-hits-file results/no_hits.tsv
```

### 5. Final Taxonomy Processing

```bash
python final_process_taxonomy.py \
  results/combined_taxonomy_results.tsv \
  results/final_taxonomy_COI.tsv \
  --marker COI \
  --min-pident 85 \
  --min-alignment 200 \
  --source gbif \
  --no-hits-file results/no_hits.tsv \
  --separate-by-phylum \
  --changes-log results/taxonomy_changes.log \
  --cache-file results/taxonomy_cache.json
```

## Marker-Specific Thresholds

Default thresholds for each marker type:

| Marker | Min Identity | Min Alignment |
|--------|--------------|--------------|
| COI    | 85%          | 200bp        |
| 18S    | 85%          | 300bp        |

## Output Files

The pipeline produces the following primary output files:

- `nmnh_taxonomy_results.tsv`: NMNH voucher taxonomy results
- `MZG_taxonomy_results.tsv`: MetaZooGene taxonomy results
- `ncbi_filtered_results.tsv`: NCBI BLAST results (full mode only)
- `combined_taxonomy_results.tsv`: Combined best hits across all databases
- `final_taxonomy_COI.tsv` or `final_taxonomy_18S.tsv`: Final taxonomy assignments with post-processing
- `no_hits.tsv`: ASVs with no hits meeting criteria
- `missing_asvs.tsv`: Complete list of ASVs missing from the final taxonomy
- `missing_asvs.fasta`: FASTA file with sequences for all missing ASVs
- `taxonomy_changes.log`: Log of taxonomy changes made during backfilling
- When using `--separate-by-phylum`:
  - `final_taxonomy_metazoan_COI.tsv`: Metazoan and Rhodophyta ASVs
  - `final_taxonomy_non_metazoan_COI.tsv`: Other phyla ASVs

## Taxonomy Backfilling

The pipeline supports two sources for taxonomy backfilling:

### WoRMS (World Register of Marine Species)
- More specialized for marine organisms
- Better for marine metabarcoding studies
- Slower API (requires longer sleep times between queries)
- More authoritative for marine taxonomy

### GBIF (Global Biodiversity Information Facility)
- Broader taxonomic coverage
- Faster API (can process more queries quickly)
- May have less marine-specific information
- Requires the `pygbif` package

## Troubleshooting

### Common Issues

- **BLAST errors**: Ensure BLAST+ is properly installed and available in PATH
- **Database errors**: Check that reference database paths are correct
- **No hits found**: Try reducing identity or alignment thresholds
- **NCBI API errors**: Ensure you've provided a valid email address

### FASTA Header Length Issues

BLAST has a 50-character limit for FASTA headers. If your reference headers are too long, use this Python script to fix them:

```python
#!/usr/bin/env python3
from Bio import SeqIO
import sys

input_file = "references/18S/NMNH_18S_refs.fasta"
output_file = "references/18S/NMNH_18S_refs_fixed.fasta"

with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        if len(record.id) > 50:
            parts = record.id.split('|')
            if len(parts) >= 3:
                new_id = f"{parts[0]}|{parts[1]}"
                remaining_space = 50 - len(new_id) - 1
                if remaining_space > 0:
                    new_id += f"|{parts[2][:remaining_space]}"
            else:
                new_id = record.id[:50]
            record.id = new_id
            record.description = ""
        SeqIO.write(record, output_handle, "fasta")
```

## Database Prioritization

The pipeline can prioritize more reliable reference databases over raw percent identity:

- **Default Priority Order**: MetaZooGene > NMNH > NCBI
- This ensures higher quality, manually curated marine references are used over more general references
- To disable, use the `--no-prioritize-db` flag
- To customize priority order, use `--db-order "DB1,DB2,DB3"`

## Example Advanced Workflow

### Running Both Markers with Comparison

```bash
# Run COI pipeline
./taxonomy_master.sh --marker COI --db-mode full \
  --asv-fasta data/COI_ASVs.fasta --output-dir results_COI \
  --email your@email.com --backfill gbif --prioritize-db \
  --separate-by-phylum

# Run 18S pipeline
./taxonomy_master.sh --marker 18S --db-mode full \
  --asv-fasta data/18S_ASVs.fasta --output-dir results_18S \
  --email your@email.com --backfill gbif --prioritize-db \
  --separate-by-phylum

# Compare phylum distribution
python -c "
import pandas as pd
coi = pd.read_csv('results_COI/final_taxonomy_COI.tsv', sep='\t')
s18 = pd.read_csv('results_18S/final_taxonomy_18S.tsv', sep='\t')
print('COI phyla distribution:')
print(coi['phylum'].value_counts())
print('\n18S phyla distribution:')
print(s18['phylum'].value_counts())
"
```

### Find Missing ASVs and Check Coverage

The pipeline automatically generates files with missing ASVs, but you can also run a custom analysis:

```bash
# Extract specific marker data
python -c "
import pandas as pd
from Bio import SeqIO
import sys

# Load final taxonomy results
results = pd.read_csv('results_COI/final_taxonomy_COI.tsv', sep='\t')

# Calculate percent identity statistics
print(f'Average percent identity: {results[\"pident\"].mean():.2f}%')
print(f'Minimum percent identity: {results[\"pident\"].min():.2f}%')
print(f'Maximum percent identity: {results[\"pident\"].max():.2f}%')

# Count by database source
print('\\nDatabase source distribution:')
print(results['database'].value_counts())

# Count by phylum
print('\\nPhylum distribution:')
print(results['phylum'].value_counts().head(10))
"
```

## Citation

If using this pipeline in your research, please cite:
- MetaZooGene: [https://metazoogene.org/](https://metazoogene.org/)
- BLAST: Altschul SF, et al. (1990). Basic local alignment search tool. J Mol Biol. 215:403-410.
- [ARMS Manuscript]
- [PRIMER CITATIONS]

## Complete Pipeline Workflow

The taxonomy pipeline follows these steps:

1. **NMNH Processing**:
   - Runs BLAST against NMNH voucher database
   - Processes BLAST results to extract taxonomy
   - Applies marker-specific thresholds

2. **MetaZooGene Processing**:
   - Runs BLAST against MetaZooGene database
   - Parses taxonomy from MZG headers
   - Filters based on identity and alignment thresholds

3. **NCBI Processing** (full mode only):
   - Runs BLAST against NCBI nt database
   - Filters results based on criteria
   - Creates a standardized output format

4. **Best Hit Selection**:
   - Combines results from all databases
   - Selects best hit for each ASV based on configurable criteria
   - Prioritizes by database source and/or percent identity

5. **Final Taxonomy Processing**:
   - Filters based on marker-specific thresholds
   - Removes environmental sequences
   - Backfills missing taxonomy using WoRMS or GBIF
   - Forward-fills taxonomy ranks with appropriate symbols
   - Optionally separates by phylum

6. **Missing ASV Analysis**:
   - Identifies ASVs missing from final taxonomy
   - Creates a list of missing ASVs
   - Generates a FASTA file with missing sequences
