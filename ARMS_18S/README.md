# 18S rRNA ARMS Metabarcoding: Post-Bioinformatics Analysis Pipeline

**Author:** Nicole Pittoors, Lehigh Oceans Research Center, Lehigh University  
**Manuscript:** *"Environmental filtering shapes patch dynamics across isolated mesophotic reefs"*

---

## Overview

This repository contains the bioinformatics analysis pipeline for 18S rRNA metabarcoding data from Autonomous Reef Monitoring Structures (ARMS) deployed across mesophotic reef banks in the northwestern Gulf of Mexico. It covers all steps from the LULU-curated, taxonomically assigned ASV table through community diversity and differential abundance analyses. The pipeline is divided into two sequential R Markdown scripts:

1. **`FINAL_18S_PreTaxa.Rmd`** — Pre-taxonomic (ASV-level) biodiversity analyses: alpha diversity, beta diversity, ordination, PERMANOVA, dbRDA, and distance-based community assembly tests
2. **`PostTaxa_18S.Rmd`** — Post-taxonomic analyses: relative abundance visualization, differential abundance (ANCOM-BC2), and functional trait analysis

**Analysis scope:** 114 ARMS samples (36 ARMS units × 3 size fractions) across 12 sites on six reef banks.  
**Size fractions:** Sessile (encrusting), 500 µm motile, 100 µm motile.

---

## Data Availability

| Data | File | Location |
|---|---|---|
| Raw sequences (COI + 18S) | — | NCBI SRA, BioProject [PRJNA1159220](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1159220) |
| 18S ASV table | `18S_metazoo_ASVtab.csv` | `ARMS_metabarcoding/Data/` |
| COI ASV table | `COI_metazoan_ASVtab.csv` | `ARMS_metabarcoding/Data/` |
| 18S metadata | `metadata_18S_FINAL.csv` | `ARMS_metabarcoding/Data/` |
| COI metadata / environmental variables | `metadata_ARMS_env_COI.csv` | `ARMS_metabarcoding/Data/` |
| 18S taxonomy table | `taxa_metazoan_18S.tsv` | `ARMS_metabarcoding/Data/` |
| COI taxonomy table | `taxa_metazoan_COI_final.csv` | `ARMS_metabarcoding/Data/` |
| Functional trait mapping table | `trait_mapping_QC_v4.csv` | `ARMS_metabarcoding/Data/` |

---

## Complete Workflow

```
Step 1a: Sequence Processing  [QC_trim_DADA2.md]
            ↓
Step 1b: Taxonomic Assignment  [ARMS_Assign_Taxonomy_Pipeline]
            ↓
Step 2: Pre-Taxonomic Diversity Analysis  [FINAL_18S_PreTaxa.Rmd]
            ↓
Step 3: Taxonomic & Differential Abundance Analysis  [PostTaxa_18S.Rmd]
```

### Step 1a — Sequence Processing (DADA2 + LULU)

Raw paired-end Illumina amplicon sequences were processed using **DADA2** for quality filtering, denoising, ASV inference, and chimera removal, preceded by quality control with **FastQC/MultiQC** and primer trimming with **Cutadapt**.

> **The full pipeline for both COI and 18S markers is documented in [`QC_trim_DADA2.md`](./QC_trim_DADA2.md) in this repository.** Parameters for all steps are recorded there. Key parameters are also summarized below for reference.

#### Primer Trimming (Cutadapt)

Primers were removed using **Cutadapt** prior to DADA2 processing. Reads with no detectable forward primer were discarded (`--discard-untrimmed`).

**18S V4** 

| Parameter | Value |
|---|---|
| Forward primer (5′→3′) | `CCAGCASCYGCGGTAATTCC` (required) |
| Forward primer RC | `TCATYRATCAAGAACGAAAGT` (optional, internal) |
| Reverse primer (5′→3′) | `ACTTTCGTTCTTGATYRATGA` (required) |
| Reverse primer RC | `GGAATTACCGCRGSTGCTGG` (optional, internal) |
| Adapter syntax | Anchored linked adapters (`^PRIMER;required;...RC$;optional`) |
| `--overlap` / `-O` | 5 bp |
| `--minimum-length` | 5 bp |
| `--pair-filter` | `any` |
| `--action` | `trim` |
| `--discard-untrimmed` | Yes |
| Input | Paired-end `.fastq.gz`; processed in a `for` loop over `*_R1_001.fastq.gz` |

*Amplicon target:* ~536 bp fragment of the 18S V4 region  
*References:* Piredda et al. 2017; Tragin et al. 2018

**COI**

Primers used were the **Geller/Leray** COI primers targeting a ~313 bp fragment of mitochondrial cytochrome c oxidase subunit I, designed for marine invertebrates and metazoan metabarcoding.

| Parameter | Value |
|---|---|
| Forward primer (5′→3′) | `GGWACWGGWTGAACWGTWTAYCCYCC` (anchored) |
| Reverse primer (5′→3′) | `TAIACYTCIGGRTGICCRAARAAYCA` (anchored; contains inosine) |
| `-n` | 1 |
| `--overlap` / `-O` | 5 bp |
| `--minimum-length` | 5 bp |
| `--pair-filter` | `any` |
| `--action` | `trim` |
| `--discard-untrimmed` | Yes |
| Input | Paired-end `.fastq.gz`; processed in a `for` loop over `*_R1_001.fastq.gz` |

*Amplicon target:* ~313 bp fragment of COI (390 bp including adapters/primers)  
*References:* Geller et al. 2013 (*Mol Ecol Res* 13(5):851–861); Leray et al. 2013 (*Front Zool* 10(34):1–14)

#### Quality Filtering and Denoising (DADA2)

**COI** 

| Parameter | Value |
|---|---|
| `truncLen` | c(220, 190) — forward, reverse |
| `maxEE` | c(2, 3) — forward, reverse |
| `errorEstimationFunction` | `loessErrfun` (quality-score-ignoring) |
| `pool` | `"pseudo"` (pseudopooling; enables singleton detection) |
| `multithread` | `TRUE` |
| Chimera removal | `removeBimeraDenovo()` |

**18S** 

| Parameter | Value |
|---|---|
| `truncLen` | c(220, 180) — forward, reverse |
| `maxEE` | c(2, 3) — forward, reverse |
| `pool` | `"pseudo"` (pseudopooling; enables singleton detection) |
| `multithread` | `TRUE` |
| Chimera removal | `removeBimeraDenovo()` |

> **QC note:** Although the 18S V4 amplicon is nominally ~536 bp, truncation lengths of c(220, 180) produced a mean paired-read merge rate of 81.8% across 164 samples (range: 45.6–97.1%), confirming sufficient overlap for merging. The effective amplicon length is shorter than the nominal target across the metazoan assemblage sampled, consistent with known V4 length variability among marine invertebrate phyla. Two samples (both 500 µm fraction) fell below 50% merge rate; all other samples exceeded 65%. Overall read retention from input to non-chimeric reads was 65.2%.

#### LULU Curation
Following DADA2, ASV tables for both markers were curated using **LULU** to remove erroneous ASVs arising from co-amplification artifacts. The LULU curation script (`lulu_curation.R`) is available in this repository.

### Step 1b — Taxonomic Assignment

Taxonomic classification was performed using a custom hierarchical BLAST pipeline against curated marine reference databases, with taxonomy backfilling via WoRMS or GBIF and filtering to retain Metazoa and Rhodophyta only.

> 
> **[https://github.com/npittoors/CYCLE-ARMS-Community-Analyses/tree/main/ARMS_Assign_Taxonomy_Pipeline](https://github.com/npittoors/CYCLE-ARMS-Community-Analyses/tree/main/ARMS_Assign_Taxonomy_Pipeline)**

**Output of Steps 1a–1b:** `18S_metazoo_ASVtab.csv` — 7,089 ASVs × 114 ARMS samples (metazoan + Rhodophyta only)

---

## Required R Packages

### Pre-Taxonomic Script (`FINAL_18S_PreTaxa.Rmd`)

```r
# Data management
library(tidyr)
library(plyr)        # Load before dplyr
library(dplyr)
library(reshape2)

# Diversity & community analysis
library(vegan)
library(phyloseq)

# Statistical modeling
library(lme4)
library(multcomp)
library(emmeans)
library(car)
library(MuMIn)
library(performance)
library(coin)
library(broom.mixed)

# Visualization
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(ggsci)
library(patchwork)
library(grid)

# Optional
library(openxlsx)
```

### Post-Taxonomic Script (`PostTaxa_18S.Rmd`)

```r
# Core
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)

# Differential abundance
library(ANCOMBC)
library(conflicted)
library(foreach)
library(rngtools)
library(reshape2)

# Visualization
library(viridis)
library(gridExtra)
library(patchwork)
library(ggpattern)   # For combined 18S + COI plots
```

---

## Input Files

### Primary Inputs

| File | Description | Used In |
|---|---|---|
| `18S_metazoo_ASVtab.csv` | LULU-curated metazoan + Rhodophyta ASV count table (7,089 ASVs × 114 samples) | PreTaxa |
| `metadata_ARMS_env.csv` | Sample metadata with environmental variables (depth, turbidity, coordinates, HOBO data) | PreTaxa |
| `biotic_percs_18S_4th.Rdata` | ASV proportions (4th-root transformed, ordered west to east) | PostTaxa |
| `biotic_clean_18S.Rdata` | ASV raw count table (zero-read ASVs removed) | PostTaxa |
| `metadata_18S.Rdata` | Sample metadata including `turbidity_std_rank`, `sitelocality`, `fraction`, `ARMS_ID_corrected` | PostTaxa |
| `taxa_metazoan_18S.tsv` | Metazoan taxonomy table (rows = ASVs; columns = phylum through species) | PostTaxa |
| `missing_asvs.tsv` | ASV IDs with no taxonomy database hits (for unclassified read analysis) | PostTaxa |
| `18S_lulu_ASVtab.csv` | Complete post-LULU ASV count table (all samples including non-ARMS) | PostTaxa |
| `final_taxonomy_results.tsv` | Complete taxonomy assignments (metazoan + non-metazoan) | PostTaxa |
| `trait_mapping_QC_v4.csv` | Functional trait table mapping families to nutritional strategies | PostTaxa |
| `COI_phylum_turbidity_simple.csv` | Exported COI phylum turbidity results (for combined 18S + COI plots) | PostTaxa |
| `all_environmental_significant_results_18S.csv` | Exported from previous ANCOM-BC2 run (used in recovery chunk) | PostTaxa |
| `ancombc_raw_counts_physeq_18S.RData` | Saved phyloseq object from ANCOM-BC2 run (optional, for recovery) | PostTaxa |

### Recommended Intermediate R Objects (for Reproducibility)

| File | Description | Allows Starting From |
|---|---|---|
| `metazoan_ASVtab_18S.Rdata` | Clean ASV matrix | Part 2 (skip data wrangling) |
| `metadata_ARMS.Rdata` | Metadata synchronized with ASV table | Part 2 |
| `d_r17k.Rdata` | Rarefied phyloseq object (17,869 reads/sample) | Part 4 (alpha diversity) |
| `18S_all_data_subsets.RData` | All fraction-specific matrices and metadata | Part 8 (beta diversity) |
| `metadata_18S_FINAL.RData` | Final synchronized metadata for all analyses | Any part |

---

## Study Design

12 sites across six reef banks, ordered from high to low turbidity:

| Site Code | Site Locality |
|---|---|
| MCG | McGrail_deep |
| EFGshal | EFGB_shallow |
| ALDshal | Alderdice_shallow |
| BRIshal | Bright_shallow |
| DIAback | Diaphus_deep2 |
| EFGdeep | EFGB_deep |
| DIAcoral | Diaphus_deep1 |
| BRIback | Bright_deep2 |
| BRIcoral | Bright_deep1 |
| ALDcoral | Alderdice_deep1 |
| STE | Stetson_shallow |
| ALDback | Alderdice_deep2 |

---

## Part 1: Pre-Taxonomic Diversity Analysis (`FINAL_18S_PreTaxa.Rmd`)

> Analyses use ASV-level data without family/genus/species resolution, providing a taxonomically agnostic view of community structure.

### Data Filtering Applied Upstream
Retained: All Metazoa phyla, Rhodophyta  
Removed: Protists, fungi, non-rhodophyte algae, terrestrial/freshwater taxa, controls, eDNA water samples  
Result: 6,332 ASVs retained across 114 ARMS samples after zero-abundance removal

### Structure

**Part 1 — Data Import and Phyloseq Object Creation**  
Imports ASV table and metadata, synchronizes sample IDs, removes controls and eDNA columns, builds primary phyloseq object. Removes zero-abundance ASVs (7,089 → 6,332).

**Part 2 — Data Quality Assessment**  
Raw read counts, library size distributions, and rarefaction curves by site and fraction to assess sampling completeness.

*Output files:* `18S_RawReads_Sample_Site.pdf`, `18S_Reads_ASVs_Samples.pdf`, `18S_rarefaction_by_site.pdf`, `18S_rarefaction_by_fraction.pdf`

**Part 3 — Rarefaction**  
Rarefied to 17,869 reads (minimum library size). Seed = 42. Retains 5,852 / 7,089 ASVs (82.6%).

**Part 4 — Alpha Diversity Metrics**  
Observed richness, Shannon (H'), Simpson, Pielou's evenness (J'), Chao1, ACE, Margalef, and mean Jaccard dissimilarity calculated from rarefied data via `vegan`. Visualized as violin plots by site and fraction.

*Output files:* `18S_diversity_indices_rarefied.csv`, `18S_alphaDiv_Locality_violin_rarefied.pdf`, `18S_alphaDiv_Fraction_violin_rarefied.pdf`

**Part 5 — Linear Mixed-Effects Models (LMMs)**  
Three analyses accounting for paired fraction structure (random effect: ARMS unit):
- *Analysis 1:* Fraction effects — `Diversity ~ fraction + (1|ARMS_ID)`; sessile vs. pooled motile contrast
- *Analysis 2:* Environmental effects — `Diversity ~ depth + turbidity + fraction + (1|ARMS_ID)`; marginal R² reported
- *Analysis 3:* Site effects — `Diversity ~ fraction + sitelocality + (1|ARMS_ID)`; Tukey HSD post-hoc if significant

*Output files:* `18S_Analysis1_Fraction_Effects_4metrics.csv`, `18S_Analysis2_Environmental_Effects_4metrics.csv`, `18S_Analysis3_Site_Effects_4metrics.csv`, conditional site comparison CSVs

**Part 6 — Descriptive Statistics**  
Summary tables (mean ± SE) by fraction and site for all four diversity metrics. Includes site rankings and full sample-level data.

*Output files:* `TableS2A_Supplement_18S_Fraction_Stats.csv`, `TableS2_18S_Site_[metric].csv`, `TableS2_18S_Sample_Level_Data.csv`

**Part 7 — 4th Root Transformation for Beta Diversity**  
Converts raw counts to relative proportions then applies 4th-root transformation (`proportion^0.25`) to downweight numerically dominant ASVs. Creates fraction-specific subsets (all, sessile, motile, 100 µm, 500 µm) with raw, proportion, 4th-root, and presence/absence matrices.

*Output files:* `18S_all_data_subsets.RData`, `biotic_clean_18S.RData`, `biotic_percs_18S_4th.RData`, `metadata_18S_FINAL.RData/.csv`

**Part 8 — Beta Diversity Ordination (NMDS)**  
NMDS on Bray-Curtis dissimilarity (k=2, max 1,000 iterations, 20 random starts). Ordinations for all samples, by site, by fraction, sessile only, and motile only. Visualized with 95% ellipses and environmental vectors.

*Output files:* `18S_NMDS_AllSamples_BySite.pdf`, `18S_NMDS_AllSamples_ByFraction.pdf`, `18S_NMDS_Sessile.pdf`, `18S_NMDS_Motile.pdf`

**Part 9 — PERMANOVA**  
Bray-Curtis dissimilarity on 4th-root data, 999 permutations. Four model types:
- Environmental baseline: `~ depth + turbidity`
- Spatial: `~ latitude + longitude`
- Combined: `~ depth + turbidity + latitude + longitude`
- Site identity: `~ sitelocality`

Conditional tests partition pure environmental vs. pure spatial variance. Pre-tested with PERMDISP. Run separately for all samples, sessile, motile, 100 µm, 500 µm.

*Output files:* `18S_PERMANOVA_[subset].csv`

**Part 10 — dbRDA**  
Constrained ordination with hierarchical model series (individual predictors → combined → conditional). Restricted permutations within ARMS units. Reports sequential R² per predictor and unique contribution of turbidity beyond depth.

**Part 11 — Distance Matrix Analyses**  
Three complementary tests using Bray-Curtis (community), Euclidean on z-standardized env. variables (environmental distance), and least-cost geographic distance:
- Simple Mantel (9,999 permutations)
- Partial Mantel controlling for environment or geography
- Multiple Regression on Distance Matrices (MRM)

Run separately for all samples, sessile, and motile fractions.

*Output files:* `18S_Mantel_Results_nopp.csv`, `18S_Partial_Mantel_Results_nopp.csv`, `18S_MRM_Results_nopp.csv`

---

## Part 2: Post-Taxonomic Analysis (`PostTaxa_18S.Rmd`)

### Structure

**Section 1 — Setup and Data Import**  
Imports all input files, aligns taxonomy to ASV table, builds phyloseq object using 4th-root transformed proportions. Output directory: `./Taxa_Output_18S/`

**Section 2 — Relative Abundance Plots (Classified Reads)**  
Stacked bar plots of top 20 metazoan/rhodophyte phyla per sample (ordered by turbidity) and aggregated by site.

*Output files:* `18S_Top20_Phyla_RelativeAbundance_bySample.pdf`, `18S_Top20_Phyla_RelativeAbundance_BySite.pdf`

**Section 3 — Relative Abundance Plots (Including Unclassified Reads)**  
Adds ASVs with no taxonomy hits as an "Unclassified" category.

*Output files:* `18S_Top20_Phyla_RelativeAbundance_bySample_WITH_UNCLASSIFIED.pdf`, `18S_Top20_Phyla_RelativeAbundance_BySite_WITH_UNCLASSIFIED.pdf`

**Section 4 — ANCOM-BC2 Differential Abundance**  
Uses raw count data. Three sample subsets (all, motile, sessile) × two taxonomic levels (phylum, family) × three predictor sets:

| Analysis | Formula | Test Type |
|---|---|---|
| Fraction effects | `~ fraction` | Global + pairwise |
| Site effects | `~ sitelocality` | Global + pairwise |
| Environmental effects | `~ depth + turbidity_std_rank` | Per-variable |

Key parameters: `prv_cut = 0.10`, `p_adj_method = "BH"`, `pseudo_sens = TRUE`

> ⚠️ ANCOM-BC2 runs are computationally intensive (10–30+ min per call). A recovery chunk reloads results from previously exported CSVs and a saved `.RData` phyloseq object if the session crashes.

**Section 5 — Environmental Effects Visualization**  
Phylum relative abundance barplots ordered by turbidity (cividis palette) and depth (mako palette) for motile and sessile fractions. Log₂ fold-change (LFC) barplots for all significant phyla.

*Output files:* `Phylum_LFC_[subset]_Clean_18S.pdf/.png`, `[Fraction]_[Gradient]_Phyla_[palette].pdf`

**Section 6 — Combined 18S + COI Phylum LFC Plots**  
Overlays 18S turbidity effects (solid bars) with COI results (striped bars, from `COI_phylum_turbidity_simple.csv`).

*Output files:* `Combined_18S_COI_Phylum_[subset].pdf/.png`

**Section 7 — Functional Trait Analysis**  
Maps families to nutritional strategy via `trait_mapping_QC_v4.csv`. Tests whether feeding mechanism predicts turbidity sensitivity, with focus on Fine Dead-End Microfilterers (FDM; sponges — predicted to decrease with turbidity) vs. Active Dead-End Sieves (ADS; bivalves and barnacles — predicted to increase with turbidity).

*Output files:* `functional_group_summary.csv`, `suspension_mechanism_summary.csv`, `functional_analysis_combined.pdf`, and related plots

---

## Key Variables Reference

### Phyloseq Objects

| Object | Description |
|---|---|
| `cycle18Sdat` | Original non-rarefied phyloseq (6,332 ASVs, 114 samples) |
| `cycle18Sdat0` | Pre-filtering with zero-abundance ASVs (7,089 ASVs) |
| `d_r17k` | Rarefied phyloseq (5,852 ASVs, 17,869 reads/sample) |
| `cycle18Sdat_4th` | 4th-root transformed non-rarefied |
| `d_r17k_4th` | 4th-root transformed rarefied |

### Metadata Versions

| Version | Object | Key Additions |
|---|---|---|
| 1 | `metadata_ARMS_env.csv` | Raw input |
| 2 | `metadata_ARMS` | + `fileID_18S` for ASV table matching |
| 3 | `metadata_18S` | + abundance/richness columns |
| 4 | `meta_dr17k` | + 9 diversity metrics, `ARMS_ID_corrected` |
| 5 | `meta_consistent` | Alias of v4; used for all final analyses |




