# 18S rDNA Pre-Taxonomic Assignment Analysis Pipeline

## Overview
This R markdown pipeline performs comprehensive **pre-taxonomic** biodiversity analyses on 18S rDNA metabarcoding data from ARMS (Autonomous Reef Monitoring Structures) deployed across mesophotic Gulf of Mexico reef banks. Analyses are conducted at the ASV level before family/species-level taxonomic assignment, providing a taxonomically-agnostic view of community structure and diversity patterns.

**Analysis Scope:** 114 ARMS samples (36 ARMS units × 3 size fractions) across 12 sites and 3 size fractions (sessile, 500 µm motile, 100 µm motile).

**Important Note on Data Filtering:** The input ASV table (`18S_metazoo_ASVtab.csv`) contains **only marine benthic metazoans and Rhodophyta (red algae)**. Taxonomic classification was performed upstream (see taxonomic assignment pipeline) to confirm ASV identities before this analysis. All non-metazoan eukaryotes (protists, fungi, non-rhodophyte algae, etc.), controls, and eDNA samples were removed prior to this workflow. This ensures that "pre-taxa" analyses are still biologically meaningful, focusing on reef-associated multicellular organisms while avoiding taxonomic resolution biases at family/genus/species levels.


---

## Relationship to Taxonomic Assignment Pipeline

This "pre-taxa" analysis is **downstream** of initial taxonomic classification but **upstream** of fine-scale taxonomic analyses. Here's how it fits into the complete workflow:

### Complete Analysis Workflow:

**Step 1: Bioinformatics & Taxonomic Classification** (Upstream - Separate Pipeline)
- Raw sequence processing (DADA2)
- ASV inference and chimera removal
- LULU curation for co-amplification artifacts
- **Taxonomic classification to phylum level minimum**
  - Custom hierarchical pipeline: NMNH vouchers → MetaZooGene → NCBI
  - Forward-filling missing ranks using WoRMS API
  - Marine habitat verification
- **Filter to metazoans + Rhodophyta only**
- Remove controls, blanks, eDNA samples
- Output: `18S_metazoo_ASVtab.csv` (input for this pipeline)

**Step 2: Pre-Taxonomic Diversity Analysis** (This Pipeline)
- ASV-level alpha diversity (richness, Shannon, Simpson, evenness)
- ASV-level beta diversity (NMDS, PERMANOVA, dbRDA)
- Distance-based analyses (Mantel, MRM)
- Community structure without taxonomic labels

**Step 3: Taxonomic Analysis** (Downstream - Separate Pipeline)
- Phylum-level composition bar plots
- Family-level differential abundance (ANCOM-BC2)
- Taxon-specific environmental responses
- Functional group analyses

### Data Quality Assurance:

Taxonomic filtering ensures biological relevance while maintaining ASV-level resolution:

**What was filtered OUT (before this analysis):**
- Protists (Ciliophora, Cercozoa, etc.)
- Fungi (Ascomycota, Basidiomycota)
- Non-rhodophyte algae (Chlorophyta, Ochrophyta)
- Terrestrial/freshwater taxa
- Unclassified eukaryotes (no phylum match)
- Negative controls and extraction blanks
- eDNA water samples (analyzed separately)

**What was KEPT (for this analysis):**
- All Metazoa phyla: Porifera, Cnidaria, Platyhelminthes, Annelida, Mollusca, Arthropoda, Bryozoa, Echinodermata, Chordata, etc.
- Rhodophyta (red algae) - important calcifying reef organisms
- Both common and rare ASVs (no abundance filtering)
- All three size fractions (sessile, 500 µm, 100 µm)

**What was NOT required (intentionally):**
- Family-level classification
- Genus-level classification  
- Species-level classification
- Perfect taxonomic completeness


---

## Pipeline Structure

### **Part 1: Data Import and Phyloseq Object Creation**
- Load metazoan+Rhodophyta-filtered ASV table (LULU-curated)
  - **Upstream filtering already completed:** Taxonomic classification was performed to identify ASVs to at least phylum level
  - Only ASVs classified as **Metazoa** (all animal phyla) or **Rhodophyta** (red algae, important reef calcifiers) were retained
  - All other eukaryotes removed: protists, fungi, Chlorophyta, Ochrophyta, non-marine taxa
- Import ARMS metadata with environmental variables
- Synchronize file IDs between ASV table and metadata
- Remove controls and eDNA samples (columns 115-163)
- Create phyloseq object for integrated analyses
- Remove zero-abundance ASVs (7,089 → 6,332 ASVs retained)
  - Zero-abundance ASVs result from removing eDNA/control samples

**Key Outputs:**
- `biotic_18S` - Clean ASV matrix (6,332 metazoan/rhodophyte ASVs × 114 ARMS samples)
- `cycle18Sdat` - Primary phyloseq object
- `metadata_18S` - Enhanced metadata with diversity metrics

---

### **Part 2: Data Quality Assessment**
- **Raw read visualization** - Per-sample read counts colored by site
- **Library size distributions** - ASV and sample read summaries
- **Rarefaction curves** - Adequacy of sampling depth by site and fraction
  - Calculated at 100-read intervals up to minimum library size
  - Separate curves for 12 sites and 3 size fractions
  - Asymptote assessment for sampling completeness

**Outputs:**
- `18S_RawReads_Sample_Site.pdf`
- `18S_Reads_ASVs_Samples.pdf`
- `18S_rarefaction_by_site.pdf`
- `18S_rarefaction_by_fraction.pdf`

---

### **Part 3: Rarefaction to Standardize Library Size**
**Rationale:** 18S library sizes ranged from 17,869 to 55,729 reads. Rarefaction standardizes library size to enable fair alpha diversity comparisons.

**Rarefaction Parameters:**
- Target depth: 17,869 reads (minimum library size)
- Seed: 42 (for reproducibility)
- ASVs retained: 5,852 / 7,089 (82.6%)
- Removed ASVs: primarily rare singletons/doubletons

**Key Outputs:**
- `d_r17k` - Rarefied phyloseq object (17,869 reads/sample)
- Rarefaction report with ASV retention statistics

---

### **Part 4: Alpha Diversity Metrics**
Calculate diversity indices from rarefied data using `vegan` package:

**Primary Metrics:**
- **Observed richness** - Number of ASVs per sample
- **Shannon diversity (H')** - Entropy-based integrated diversity metric
- **Simpson diversity** - Probability-based integrated diversity metric  
- **Pielou's evenness (J')** - Shannon divided by maximum possible Shannon

**Additional Metrics:**
- Chao1 richness estimate
- ACE richness estimate
- Margalef diversity index
- Mean Jaccard dissimilarity per sample

**Visualization:**
- Violin plots by site (12 sites)
- Violin plots by fraction (sessile, 500 µm, 100 µm)
- Mean ± SE overlays with diamond symbols

**Outputs:**
- `18S_diversity_indices_rarefied.csv` - All metrics per sample
- `18S_alphaDiv_Locality_violin_rarefied.pdf`
- `18S_alphaDiv_Fraction_violin_rarefied.pdf`

---

### **Part 5: Linear Mixed-Effects Models (LMMs)**
Statistical framework accounting for paired sampling structure (3 fractions per ARMS unit).

**Analysis 1: Size Fraction Effects**
- **Model:** `Diversity ~ fraction + (1|ARMS_ID)`
- **Test:** Type II Wald χ² for overall fraction effect
- **Contrasts:** Planned comparison of sessile vs. pooled motile (100 µm + 500 µm)
- **Purpose:** Test hypothesis that sessile and motile communities differ fundamentally
- **Metrics:** Shannon, Simpson, richness, evenness

**Analysis 2: Environmental Gradient Effects**
- **Model:** `Diversity ~ depth + turbidity + fraction + (1|ARMS_ID)`
- **Predictors:**
  - Depth (54-85 m continuous)
  - Turbidity rank (standardized visibility estimates, 1=clear to 13=turbid)
- **Purpose:** Quantify how mesophotic depth and benthic nepheloid layer influence diversity
- **Output:** Marginal R² (variance explained by fixed effects)

**Analysis 3: Spatial Variation (Site Effects)**
- **Model:** `Diversity ~ fraction + sitelocality + (1|ARMS_ID)`
- **Test:** Type II Wald χ² for site effect
- **Post-hoc:** Conditional Tukey HSD pairwise comparisons if p < 0.05
- **Purpose:** Identify spatial heterogeneity across 12 sampling locations
- **Degrees of freedom:** Kenward-Roger approximation

**Model Diagnostics:**
- Random effect variance inspection
- Marginal R² via `MuMIn::r.squaredGLMM()`
- Residual plots (linearity, homoscedasticity, normality)

**Outputs:**
- `18S_Analysis1_Fraction_Effects_4metrics.csv`
- `18S_Analysis2_Environmental_Effects_4metrics.csv`
- `18S_Analysis3_Site_Effects_4metrics.csv`
- `18S_Site_Comparisons_[metric].csv` (conditional on significance)
- `18S_Significant_Site_Comparisons_[metric].csv`

---

### **Part 6: Descriptive Statistics**
Summary tables with mean ± SE format.

**Table S2A Components:**
1. **By Fraction:** Mean, SD, SE, Min, Max for all 4 metrics across 3 fractions
2. **By Site:** Mean, SD, SE for all 4 metrics across 12 sites
3. **Site Rankings:** Sites ranked by richness, evenness, and Shannon diversity
4. **Site × Fraction:** Richness and evenness breakdown for each combination
5. **Sample-Level Data:** Full dataset for transparency and verification

**Outputs:**
- `TableS2A_Supplement_18S_Fraction_Stats.csv`
- `TableS2_18S_Site_Richness.csv`
- `TableS2_18S_Site_Evenness.csv`
- `TableS2_18S_Site_Shannon.csv`
- `TableS2_18S_Site_Fraction_Detail.csv`
- `TableS2_18S_Sample_Level_Data.csv`

---

### **Part 7: 4th Root Transformation for Beta Diversity**
**Rationale:** Down-weight abundant ASVs to reveal community-wide patterns rather than being dominated by a few numerically abundant taxa. The 4th root transformation is mathematically equivalent to square-rooting twice and is more effective than Hellinger transformation for highly skewed metabarcoding data.

**Transformation Steps:**
1. Remove zero-abundance ASVs across all samples
2. Convert raw counts to relative abundance (proportions sum to 1)
3. Apply 4th root: `transformed_value = proportion^0.25`
4. Create transformed phyloseq objects

**Data Products:**
- `cycle18Sdat_4th` - Non-rarefied 4th root transformed phyloseq
- `d_r17k_4th` - Rarefied 4th root transformed phyloseq

**Fraction Subsets:**
Created for motile vs. sessile comparisons:
- All samples combined
- Sessile only (biofilm + encrusting organisms)
- Motile only (100 µm + 500 µm combined)
- 100 µm only
- 500 µm only

**For Each Subset:**
- Raw count matrix (ASVs × samples)
- Relative abundance matrix (proportions)
- 4th root transformed matrix (for ordination)
- Presence/absence matrix (for Jaccard)
- Metadata subset

**Outputs:**
- `18S_all_data_subsets.RData` - Complete workspace with all subsets
- `biotic_clean_18S.RData` - Clean full ASV matrix
- `biotic_percs_18S_4th.RData` - 4th root transformed matrix
- `metadata_18S_FINAL.RData` / `.csv` - Final metadata

---

### **Part 8: Beta Diversity Ordination (NMDS)**
Non-metric multidimensional scaling on Bray-Curtis dissimilarity matrices calculated from 4th root transformed data.

**NMDS Parameters:**
- Distance metric: Bray-Curtis
- Dimensions: k = 2
- Maximum iterations: 1,000
- Random starts: 20
- Convergence criterion: stress < 0.20

**Ordination Types:**
1. **All Samples** - Full dataset (n = 114)
2. **By Site** - 12 separate site ordinations
3. **By Fraction** - Sessile, 500 µm, 100 µm separate
4. **Sessile Only** - Biofilm and encrusting taxa (n = 38)
5. **Motile Only** - Free-living cryptofauna (n = 76)

**Visualization:**
- Color by: site, fraction, depth, turbidity
- Shape by: fraction
- 95% confidence ellipses by group
- Environmental vector fitting (depth, turbidity, salinity, temperature)

**Outputs:**
- `18S_NMDS_AllSamples_BySite.pdf`
- `18S_NMDS_AllSamples_ByFraction.pdf`
- `18S_NMDS_Sessile.pdf`
- `18S_NMDS_Motile.pdf`

---

### **Part 9: PERMANOVA (Community Composition)**
Permutational Multivariate Analysis of Variance testing effects of environmental and spatial variables on community composition using Bray-Curtis dissimilarity from 4th root transformed data.

**Model Types:**

**Model 1: Environmental Baseline**
- Formula: `Community ~ depth + turbidity`
- Purpose: Quantify environmental filtering effects
- Permutations: 999

**Model 2: Spatial Variables**
- Formula: `Community ~ latitude + longitude`
- Purpose: Test for geographic distance effects
- Interpretation: Dispersal limitation if significant

**Model 3: Combined Model**
- Formula: `Community ~ depth + turbidity + latitude + longitude`
- Purpose: Partition variance between environment and space

**Model 4: Site Identity**
- Formula: `Community ~ sitelocality`
- Purpose: Test for site-specific effects beyond measured variables
- Interpretation: Unmeasured local processes or microhabitat heterogeneity

**Conditional Tests (Variance Partitioning):**
- `Environment | Geography` - Pure environmental effects controlling for spatial autocorrelation
- `Geography | Environment` - Pure spatial effects (dispersal limitation) controlling for environment
- `Site | All` - Local heterogeneity beyond all regional predictors

**Separate Analyses For:**
- All samples combined
- Sessile fraction only
- Motile fraction only (100 µm + 500 µm combined)
- 100 µm fraction only
- 500 µm fraction only

**Pre-test Diagnostics:**
- PERMDISP (test for homogeneity of multivariate dispersions)
- If PERMDISP significant → use Type III SS or interpret with caution

**Key Outputs:**
- `18S_PERMANOVA_AllSamples.csv`
- `18S_PERMANOVA_Sessile.csv`
- `18S_PERMANOVA_Motile.csv`
- `18S_PERMANOVA_100um.csv`
- `18S_PERMANOVA_500um.csv`
- R² values (variance explained)
- F-statistics and p-values

---

### **Part 10: Distance-Based Redundancy Analysis (dbRDA)**
Constrained ordination partitioning community variation among environmental and spatial predictors.

**Hierarchical Model Series:**
1. Individual environmental drivers (depth, turbidity)
2. Combined environmental model (depth + turbidity)
3. Spatial model (latitude + longitude)
4. Full model (environment + space)

**Conditional Analyses:**
- `Geography | Environment` - Spatial effects after removing environmental variance
- `Environment | Geography` - Environmental effects after removing spatial autocorrelation
- `Turbidity | Depth` - Unique contribution of nepheloid layer beyond depth-associated gradients

**Permutation Strategy:**
- 999 permutations
- Restricted permutations for paired ARMS fractions (permute within ARMS unit)
- Term-by-term ANOVA for sequential variance partitioning

**Outputs:**
- dbRDA ordination plots with environmental vectors
- Variance partitioning diagrams
- Sequential R² contributions per predictor

---

### **Part 11: Distance Matrix Analyses**
Complementary hypothesis tests using pairwise dissimilarity matrices to test environmental filtering vs. dispersal limitation.

**Distance Matrices:**
- **Community dissimilarity:** Bray-Curtis on 4th root data
- **Environmental distance:** Euclidean distance on z-standardized depth + turbidity
- **Geographic distance:** Least-cost distance through navigable water (km)

**Statistical Tests:**

**1. Simple Mantel Tests (Bivariate Correlations)**
- Test: Correlation between community dissimilarity and environmental distance
- Test: Correlation between community dissimilarity and geographic distance
- Interpretation: Which distance type better predicts community turnover?
- Permutations: 9,999
- Statistic: Pearson's r

**2. Partial Mantel Tests (Independent Effects)**
- `Environment | Geography` - Environmental effects controlling for geographic distance
- `Geography | Environment` - Geographic effects controlling for environmental distance
- Interpretation: Evidence for environmental filtering vs. dispersal limitation
- Permutations: 9,999

**3. Multiple Regression on Distance Matrices (MRM)**
- Model: `Community_Dissimilarity ~ Environmental_Distance + Geographic_Distance`
- Provides: Standardized coefficients (β), p-values, overall R²
- Advantage: Simultaneously evaluates both predictors
- Disadvantage: Assumes linear relationships

**Separate Analyses For:**
- All samples
- Sessile only
- Motile only

**Key Outputs:**
- `18S_Mantel_Results_nopp.csv`
- `18S_Partial_Mantel_Results_nopp.csv`
- `18S_MRM_Results_nopp.csv`
- Scatter plots: Community dissimilarity vs. environmental/geographic distance

---

## Software Requirements

### R Version
- R ≥ 4.0.0 recommended

### Required R Packages

**Data Management:**
```r
tidyr          # Data reshaping
plyr           # Data manipulation (load before dplyr!)
dplyr          # Data manipulation
reshape2       # Data restructuring
```

**Diversity & Community Analysis:**
```r
vegan          # Diversity indices, ordination, PERMANOVA
phyloseq       # Microbiome/metabarcoding data handling
```

**Statistical Modeling:**
```r
lme4           # Linear mixed-effects models
multcomp       # Multiple comparisons
emmeans        # Estimated marginal means and contrasts
car            # Type II/III ANOVA tests
MuMIn          # Model selection and R² calculations
performance    # Model diagnostics
coin           # Permutation tests
broom.mixed    # Tidy model outputs
```

**Visualization:**
```r
ggplot2        # Publication-quality plots
RColorBrewer   # Color palettes
gridExtra      # Multi-panel plots
scales         # Scale functions for axes
ggsci          # Scientific journal color palettes
patchwork      # Combining ggplot objects
grid           # Low-level plotting
```

**Other:**
```r
openxlsx       # Excel file export (optional)
```

```

---

## Required Input Files for Reproducibility

1. **`18S_metazoo_ASVtab.csv`** (Primary input)
   - Marine benthic metazoan + Rhodophyta ASV abundance matrix
   - Rows: ASVs (7,089 ASVs after upstream filtering)
   - Columns: Samples (114 ARMS samples)
   - Format: CSV with row names (ASV IDs)
   - **Pre-filtered for biological relevance:**
     - Taxonomic classification performed upstream (see taxonomic assignment pipeline)
     - Only ASVs classified as **Metazoa** (all animal phyla) retained
     - Rhodophyta (red algae) included as key reef calcifiers
     - Non-metazoan eukaryotes removed: protists, fungi, Chlorophyta, Ochrophyta, other algae
     - Non-marine taxa removed based on WoRMS habitat classification
     - Controls and eDNA samples removed (original n = 163 samples → 114 ARMS samples)
   - **Still "pre-taxa" because:** Analyses use ASV-level data without family/genus/species assignments, avoiding lower-rank taxonomic biases

2. **`metadata_ARMS_env.csv`** (Primary input)
   - Sample metadata with environmental variables
   - Rows: 114 samples (matches ASV table columns)
   - Key columns:
     - Sample identifiers: `fileID0_18S`, `fileID_18S`
     - Site info: `sitelocality` (12 sites), `arms` (ARMS unit ID)
     - Fraction: `fraction` (sessile, 500, 100)
     - Environmental: `depth`, `turbidity_std_rank`, `latitude`, `longitude`
     - HOBO data: `temp`, `salinity`, `DO`, `light`
   - Format: CSV with row names (Sample IDs)

3. **`FINAL_18S_PreTaxa.Rmd`** (Analysis script)
   - Complete R markdown with all analyses
   - **Note:** Contains personal file paths that need sanitization (see below)

### Intermediate R Objects (Recommended for Reproducibility):

4. **`metazoan_ASVtab_18S.Rdata`**
   - Clean ASV matrix as R object
   - Allows users to skip data wrangling steps

5. **`metadata_ARMS.Rdata`**
   - Metadata as R object with fileID corrections
   - Synchronized with ASV table

6. **`d_r17k.Rdata`**
   - Rarefied phyloseq object (17,869 reads)
   - Users can start from "Part 4: Alpha Diversity"

7. **`18S_all_data_subsets.RData`**
   - All fraction-specific matrices and metadata
   - Users can start from "Part 8: Beta Diversity"

---

## Quick Start Guide

### Option 1: Start from CSV Files (Full Pipeline)
```r
# 1. Set working directory to repository root
setwd("path/to/FINAL_18S_Pipeline")

# 2. Install required packages (if needed)
install.packages(c("tidyr", "dplyr", "vegan", "phyloseq", "ggplot2", 
                   "lme4", "emmeans", "car", "MuMIn"))

# 3. Run from beginning
# Open FINAL_18S_PreTaxa.Rmd
# Execute chunks sequentially starting from "Load ASV table..."
```

### Option 2: Start from R Objects (Skip Data Wrangling)
```r
# 1. Navigate to "START HERE" section (line ~134)
# 2. Load pre-made R objects:
load("./metazoan_ASVtab_18S.Rdata")
load("./metadata_ARMS.Rdata")

# 3. Continue with phyloseq object creation
```

### Option 3: Start from Rarefied Data (Alpha/Beta Diversity Only)
```r
# 1. Navigate to "Part 4: Alpha Diversity Metrics"
# 2. Load rarefied phyloseq:
load("./d_r17k.Rdata")
load("./meta_dr17k.Rdata")

# 3. Continue with diversity calculations
```

---

## Key Variables Reference

### Phyloseq Objects:
- `cycle18Sdat` - Original non-rarefied phyloseq (6,332 ASVs, 114 samples)
- `cycle18Sdat0` - Pre-filtering phyloseq with zero-abundance ASVs (7,089 ASVs)
- `d_r17k` - Rarefied phyloseq (5,852 ASVs, 17,869 reads/sample)
- `cycle18Sdat_4th` - 4th root transformed non-rarefied phyloseq
- `d_r17k_4th` - 4th root transformed rarefied phyloseq

### ASV Matrices:
- `biotic_18S` - Raw count matrix (ASVs × samples)
- `biotic_clean_18S` - Zero-abundance ASVs removed
- `biotic_percs_18S_4th` - 4th root transformed proportions
- `biotic_sessile`, `biotic_motile` - Fraction-specific matrices

### Metadata:
- `metadata_ARMS` - Original imported metadata
- `metadata_18S` - Enhanced with diversity metrics
- `meta_dr17k` - Metadata for rarefied data (matches `d_r17k`)
- `meta_consistent` - Final synchronized metadata for all analyses
### Metadata Evolution

| Version | File/Object | Source | Added Columns | Total Cols | Saved As | Use Case |
|---------|------------|--------|---------------|------------|----------|----------|
| 1 | `metadata_ARMS_env.csv` | External | N/A | ~44 | CSV input | **Upload to GitHub** |
| 2 | `metadata_ARMS` | Import + sync | +1 fileID_18S | ~45 | metadata_ARMS.Rdata | ASV table matching |
| 3 | `metadata_18S` | + diversity | +2 abundance/richness | ~47 | basic_18S_pretaxa_inputs.Rdata | Phyloseq creation |
| 4 | `meta_dr17k` | + rarefied diversity | +9 diversity metrics, +2 ARMS_ID | ~56 | meta_dr17k.Rdata | Alpha diversity LMMs |
| 5 | `meta_consistent` | Alias of #4 | None (copy) | ~56 | metadata_18S_FINAL.RData/CSV | **Final all analyses** |

### Sample Vectors:
- `samples18S_100` - Sample IDs for 100 µm fraction (n = 38)
- `samples18S_500` - Sample IDs for 500 µm fraction (n = 38)
- `samples18S_sessile` - Sample IDs for sessile fraction (n = 38)
- `samples18S_motile` - Combined motile sample IDs (n = 76)

---

## Citation


*[Manuscript citation once published]*

**Key References:**
- Rarefaction: McMurdie PJ, Holmes S (2014) Waste not, want not: why rarefying microbiome data is inadmissible. PLoS Comput Biol 10(4):e1003531
- Phyloseq: McMurdie PJ, Holmes S (2013) phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE 8(4):e61217
- Vegan: Oksanen J et al. (2020) vegan: Community Ecology Package. R package version 2.5-7
- PERMANOVA: Anderson MJ (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecology 26:32-46
- 4th root transformation: Clarke KR, Gorley RN (2015) PRIMER v7: User Manual/Tutorial. PRIMER-E Ltd, Plymouth

---

## Troubleshooting

### Common Issues:

**1. "cannot open file" errors**
- Ensure all input files are in the correct directories
- Check that relative paths are used (not absolute paths)
- Verify file names match exactly (case-sensitive)

**2. Sample order mismatches**
- Check that `all(colnames(ASV_table) == rownames(metadata))` returns `TRUE`
- Use `metadata <- metadata[colnames(ASV_table), ]` to reorder

**3. Rarefaction removes too many ASVs**
- Check library size distribution: `summary(sample_sums(phyloseq_object))`
- Consider less conservative rarefaction depth if appropriate

**4. NMDS fails to converge (high stress)**
- Increase `trymax` to 2000
- Check for outlier samples: `plot(sample_sums(phyloseq_object))`
- Consider removing samples with very low reads

**5. PERMANOVA p-values all = 0.001**
- Increase permutations to 9999: `permutations = 9999`
- Not a problem if effect is genuinely strong

**6. LMM singular fit warnings**
- Random effect variance near zero (expected for some metrics)
- Not problematic - fixed effects are still valid
- Report warning in manuscript methods

