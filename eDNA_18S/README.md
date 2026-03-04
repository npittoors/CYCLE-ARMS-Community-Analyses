# eDNA_18S

R Markdown analysis pipeline for 18S environmental DNA (eDNA) metabarcoding data collected from mesophotic reef banks on the Texas–Louisiana continental shelf (Gulf of Mexico). This script is part of the broader analysis supporting:

> Pittoors et al. — *Environmental filtering shapes patch dynamics across isolated mesophotic reefs* 

---

## Overview

This pipeline processes 18S eDNA metabarcoding data from 33 water column samples collected across 12 site-depth combinations at six reef banks (Stetson, EFGB, Bright, McGrail, Alderdice, Diaphus). It performs quality control, rarefaction, alpha diversity analysis, community overlap assessment, beta diversity and multivariate ordination, and ANCOM-BC II differential abundance testing against environmental gradients (turbidity, depth, productivity).

**Input:** Pre-processed phyloseq object (`eDNA18S_phyloseq_20250915.rds`)  
830 ASVs × 33 samples × 20 sample variables, 7 taxonomic ranks

---

## Analysis Workflow

### 1. Data Loading & Quality Control
- Loads phyloseq object and inspects structure (ASV/sample counts, read distribution)
- Removes zero-read ASVs
- Maps `Site_locality` labels to harmonized `sitelocality` names consistent with ARMS datasets
- Appends mean primary productivity (`pp_mean`) values per site from ARMS metadata

### 2. Raw Data Visualization
- Bar plots of raw read counts per sample and per site
- Read distribution plots (ASVs and samples, log scale)

### 3. Rarefaction
- Generates rarefaction curves by site using `vegan::rrarefy`
- Selects rarefaction depth based on minimum sample read depth
- Produces rarefied phyloseq object for all downstream analyses

### 4. Alpha Diversity
Calculates per-sample diversity indices from rarefied data:
- Observed ASVs, Shannon (H'), Simpson (1-D), Pielou's evenness (J'), Margalef richness, Jaccard mean dissimilarity, Chao1, ACE

Visualizations:
- Shannon/Simpson box plots by site
- Observed ASVs and Pielou's evenness bar plots (mean ± SE) by site
- Site diversity rankings and hotspot/coldspot identification

### 5. Community Overlap & ASV Sharing
- Per-site ASV inventory (total, unique, core)
- Pairwise Jaccard similarity matrix and heatmap
- ASV sharing distribution (number of sites per ASV)
- Core ASVs (present in ≥75% of sites) and site-unique ASVs

### 6. Beta Diversity & Multivariate Analysis
- Bray-Curtis and Jaccard dissimilarity matrices
- NMDS/PCoA ordination colored by site
- PERMANOVA testing community composition against environmental variables (turbidity, depth, productivity) and site locality
- db-RDA (distance-based redundancy analysis) for constrained ordination
- Variance partitioning

### 7. Differential Abundance — ANCOM-BC II
Tests for taxon-level differential abundance across:
- **Environmental models:** turbidity (`Visibility_std_rank`), depth, and productivity (`pp_mean`) as continuous predictors
- **Site locality model:** `sitelocality` as a categorical predictor
- Analyses run at both Phylum and Family levels
- Prevalence cutoff: taxa present in ≥20% of samples (~7/33 samples)
- Results exported as CSVs; visualization of significant effects as bar plots ordered by environmental gradient

---

## Outputs

| File | Description |
|------|-------------|
| `eDNA_18S_RawReads_Sample_Site.pdf` | Raw read counts by sample, colored by site |
| `eDNA_18S_Reads_ASVs_Samples.pdf` | Read distribution across ASVs and samples |
| `eDNA_18S_Total_Reads_by_Site.pdf` | Total reads summed by site |
| `Alpha_Div_eDNA/eDNA_18S_alphaDiv_Locality_boxplot_rarefied.pdf` | Shannon/Simpson box plots |
| `Alpha_Div_eDNA/eDNA_18S_alphaDiv_Locality_barplot_rarefied_RvP.pdf` | ASV richness and Pielou's evenness bar plots |
| `Alpha_Div_eDNA/eDNA_site_similarity_heatmap.pdf` | Pairwise Jaccard similarity heatmap |
| `eDNA_Turbidity_Phyla_Effects.pdf` | ANCOM-BC II turbidity effects (phylum) |
| `eDNA_Depth_Phyla_Effects.pdf` | ANCOM-BC II depth effects (phylum) |
| `eDNA_Turbidity_Family_Effects_Top12.pdf` | ANCOM-BC II turbidity effects (top 12 families) |
| `eDNA_Depth_Family_Effects_Top12.pdf` | ANCOM-BC II depth effects (top 12 families) |
| `eDNA_Productivity_Phyla_Effects.pdf` | ANCOM-BC II productivity effects (phylum) |

---

## Dependencies

```r
# Data management
tidyr, plyr, dplyr, reshape2

# Diversity & community ecology
vegan, phyloseq

# Differential abundance
ANCOMBC  # ANCOM-BC II

# Mixed models & statistics
lme4, multcomp, emmeans, performance, coin, broom.mixed, car, MuMIn

# Visualization
ggplot2, RColorBrewer, gridExtra, scales, ggsci, patchwork, grid, kableExtra

# Output
openxlsx
```
