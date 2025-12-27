# CoralNet Analysis Pipeline

This directory contains scripts for processing and analyzing CoralNet image annotation data from ARMS (Autonomous Reef Monitoring Structures) plates deployed during the CYCLE21 project.

## Overview

CoralNet annotations were performed on photographic images of ARMS settlement plates to quantify recruitment and percent cover of benthic organisms. This pipeline processes raw CoralNet output, integrates environmental metadata, and performs multivariate community analyses.

## Scripts

### Python Scripts

#### 1. `CoralNet_abundances.py`
**Purpose**: Primary data processing script that converts CoralNet point annotations to community matrices.

**Key Functions**:
- Parses image filenames to extract sample IDs, site localities, and plate identifiers
- Maps annotation codes to taxonomic groups (Phylum level)
- Generates phylum-level count and proportion matrices
- Integrates environmental metadata (turbidity, depth, temperature, salinity)
- Outputs R-ready data matrices for downstream analysis

**Inputs**:
- `ALL_ARMS_CoralNet_Annotations.csv` - Raw CoralNet annotation output
- `CYCLE_ARMS_metadata_FINAL.csv` - Environmental metadata

**Outputs**:
- `CoralNet_R_community_matrix_ALL.csv` - Phylum counts per plate
- `CoralNet_R_sample_metadata_ALL.csv` - Sample metadata with environmental variables

---

#### 2. `CoralNet_dataprep_postMETA.py`
**Purpose**: Post-processing script that enriches CoralNet annotations with metadata.

**Key Functions**:
- Adds taxonomic classifications (Code_deff and Phylum) to each annotation point
- Integrates site-level environmental data
- Maps site localities to environmental parameters
- Validates data completeness

**Inputs**:
- `ALL_ARMS_CoralNet_Annotations.csv` - Raw annotations
- `CoralNet_R_sample_metadata.csv` - Pre-processed metadata

**Outputs**:
- `CoralNet_AllPlates_Counts.csv` - Enriched annotation data with taxonomy and environmental variables

---

#### 3. `CN_proportion_check.py`
**Purpose**: Quality control script to validate proportion calculations.

**Key Functions**:
- Verifies that proportions sum to 1.0 for each sample
- Validates algae exclusion filters (when applied)
- Performs manual spot-checks on random samples
- Identifies problematic samples with proportion calculation errors

**Usage**:
```bash
python CN_proportion_check.py
```

**Validates**:
- Standard proportion files (all taxa)
- No-algae files (excludes Chlorophyta/Ochrophyta, retains Rhodophyta)
- Recruited-only files (excludes non-recruitment categories)

---

### R Markdown

#### `CoralNet_CYCLE21_All_Plates.Rmd`
**Purpose**: Comprehensive statistical analysis and visualization of CoralNet community data.

**Analyses Include**:
1. **Data Import & Preparation**
   - Combines replicate plates per ARMS unit
   - Calculates mean proportions and summed counts
   - Filters data to recruited organisms (excludes biofilm, sediment)

2. **Alpha Diversity**
   - Shannon diversity index
   - Simpson diversity
   - Species richness
   - Pielou's evenness

3. **Beta Diversity & Ordination**
   - PERMANOVA (turbidity, depth effects)
   - NMDS ordinations (Bray-Curtis dissimilarity)
   - Distance-based redundancy analysis (dbRDA)
   - SIMPER analysis (taxa driving community differences)

4. **Environmental Filtering**
   - Turbidity effects (binary and 3-category groupings)
   - Depth effects (shallow vs deep)
   - Taxa-environment associations
   - FDR correction for multiple testing

5. **Visualization**
   - Community composition barplots
   - NMDS ordination plots
   - Alpha diversity boxplots
   - dbRDA biplots
   - PCoA
   - Taxa abundance comparisons across environmental gradients

**Outputs**:
- Statistical test results (CSV files)
- Ordination plots (PNG/PDF)
- Manuscript-ready summary tables
- FDR-corrected results

---

## Taxonomic Classification

The pipeline uses a hierarchical classification system mapping CoralNet label codes to higher taxonomic groups:

**Major Groups Analyzed**:
- **Metazoans**: Porifera, Cnidaria, Bryozoa, Mollusca, Annelida, Arthropoda, Chordata (Ascidiacea)
- **Algae**: Rhodophyta (retained), Chlorophyta (optional exclusion), Ochrophyta (optional exclusion)
- **Other**: Biofilm, Sediment, No recruitment, Unavailable substrate

## Data Processing Workflow

```
Raw CoralNet Annotations
         ↓
CoralNet_abundances.py → Community matrices (counts & proportions)
         ↓
CoralNet_dataprep_postMETA.py → Enriched annotation data
         ↓
CN_proportion_check.py → Quality validation
         ↓
CoralNet_CYCLE21_All_Plates.Rmd → Statistical analysis & visualization
```

## Environmental Variables

The following environmental parameters are integrated from site-level metadata:
- **turbidity_std_rank**: Standardized turbidity ranking
- **turbidity_m**: Mean visibility estimate (m)(NTU)
- **depth**: Deployment depth (m)
- **temp_mean**: Mean temperature over deployment period (°C)
- **temp**: temperature during recovery
- **salinity**: Salinity (PSU)
- **latitude** / **longitude**: Geographic coordinates

## Site Localities

ARMS were deployed at 12 sites across multiple banks:
- **Diaphus Bank**: Diaphus_coral, Diaphus_background
- **Alderdice Bank**: Alderdice_coral, Alderdice_background, Alderdice_shallow
- **McGrail Bank**: McGrail_coral
- **Bright Bank**: Bright_coral, Bright_background, Bright_shallow
- **Stetson Bank**: Stetson_coral
- **East Flower Garden Bank**: EFGB_shallow, EFGB_deep

_coral and _background were later changed to _deep1 and _deep_2
## Requirements

**Python**:
- pandas
- numpy
- re, os

**R**:
- vegan (multivariate analysis)
- ggplot2, viridis, RColorBrewer (visualization)
- dplyr, tibble (data manipulation)
- gridExtra, cowplot (plot arrangement)
- ape (phylogenetic/distance methods)

## Notes

- CoralNet point annotations represent percent cover estimates (typically 100-200 points per image)
- Multiple plates per ARMS unit are combined (mean proportions, summed counts)
- 4th root transformation is applied to reduce influence of dominant taxa in ordinations
- FDR correction applied to all pairwise comparisons to control family-wise error rate

## Citation

If using these scripts, please cite:
Nicole C. Pittoors, Sarah M. Tweedt, Luke J. McCartin, Samuel A. Vohsen, Luisa Lopera, Sophia Mihalek, Jamie Lai, Kathleen Durkin, Lee Weigt, Marissa F. Nuttall, Annalisa Bracco, Christopher P. Meyer, Santiago Herrera. bioRxiv 2025.11.02.686126; doi: https://doi.org/10.1101/2025.11.02.686126

## Contact

Nicole Pittoors  
ncp220@lehigh.edu

---

*Last updated: December 2025*
