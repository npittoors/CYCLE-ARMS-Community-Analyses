# Motile 2mm Organism Analysis Pipeline

This directory contains scripts for processing and analyzing motile invertebrate fauna collected from ARMS (Autonomous Reef Monitoring Structures) during the CYCLE21 project. Organisms retained on a 2mm sieve were sorted, identified morphologically and/or molecularly, and quantified.

## Overview

Motile organisms (cryptofauna) were extracted from ARMS units, size-fractionated using a 2mm sieve, preserved, and identified to the lowest taxonomic level possible. This pipeline processes taxonomic count data, integrates environmental metadata, and performs multivariate community analyses to examine spatial patterns and environmental drivers of motile invertebrate assemblages.

## Scripts

### Python Script

#### `motile_organism_processing.py`
**Purpose**: Primary data processing script that converts taxonomic occurrence records to community matrices at multiple taxonomic levels.

**Key Functions**:
- Parses ARMS event IDs to extract sample identifiers and site localities
- Implements taxonomic forward-filling using symbols to denote lowest available classification:
  - `!` = Class level
  - `#` = Order level
  - `^` = Family level
  - `*` = Genus level
  - `**` = Species level
- Generates count and proportion matrices at multiple taxonomic resolutions
- Creates order-level aggregated matrices for broader taxonomic comparisons
- Generates embedded R script for NMDS ordination and visualization

**Inputs**:
- `CYCLE_ARMS_2mm_counts.csv` - Raw taxonomic occurrence data with individual counts

**Outputs**:
- `Motile_count_matrix.csv` - Species-level count matrix (samples × taxa)
- `Motile_proportions_with_site.csv` - Species-level relative abundance
- `Motile_order_counts_with_site.csv` - Order-level count matrix
- `Motile_order_proportions_with_site.csv` - Order-level relative abundance
- `Motile_metadata.csv` - Sample metadata for R analyses
- `Motile_processed.csv` - Full processed dataset with taxonomic forward-filling

**Taxonomic Forward-Filling Logic**:
When specimens cannot be identified to species level, the script fills lower taxonomic ranks with the lowest available classification plus an appropriate symbol. For example:
- A specimen identified only to Family Caprellidae would have:
  - `family`: "Caprellidae"
  - `genus`: "Caprellidae^"
  - `scientificName`: "Caprellidae^"

This preserves taxonomic information while allowing aggregation at different levels.

---

### R Markdown

#### `Motile_2mm_CYCLE21.Rmd`
**Purpose**: Comprehensive statistical analysis and visualization of motile organism community data.

**Analyses Include**:

1. **Multi-Level Count Table Generation**
   - Phylum-level counts (broadest resolution)
   - Order-level counts (primary analysis level)
   - Family-level counts (intermediate resolution)
   - Species-level counts (finest resolution with morphospecies)
   - Quality checks to ensure count totals are preserved across taxonomic levels

2. **Taxonomic Inventory**
   - Exports unique family list with full taxonomic hierarchy
   - Removes incomplete classifications and non-metazoan entries
   - Creates clean reference taxonomy for manuscript

3. **Alpha Diversity Metrics**
   - Shannon diversity index (entropy-based)
   - Simpson diversity (dominance-based)
   - Species richness (observed taxa)
   - Pielou's evenness (equitability)
   - Calculated at multiple taxonomic levels (order and species)

4. **Beta Diversity & Ordination**
   - PERMANOVA testing environmental effects (turbidity, depth, site)
   - NMDS ordinations using Bray-Curtis dissimilarity
   - Distance-based redundancy analysis (dbRDA) with environmental constraints
   - Convex hull plots for site and environmental groupings
   - 4th root transformation to balance rare and abundant taxa

5. **Community Composition**
   - Taxonomic bar plots showing relative abundance by site
   - Order-level and family-level compositional summaries
   - Top taxa identification and visualization
   - Comparison across environmental gradients

6. **Environmental Filtering (SIMPER Analysis)**
   - **Turbidity effects**: Binary (high/low) and 3-category groupings
   - **Depth effects**: Binary (shallow/deep) and 3-category groupings
   - Taxa driving community differences between environmental groups
   - FDR correction for multiple pairwise comparisons
   - Effect size calculations and fold-change metrics

7. **Comparative Analyses**
   - Direct comparison framework with CoralNet sessile communities
   - Evaluation of motile vs. sessile environmental responses
   - Assessment of behavioral vs. environmental determinism
   - Manuscript-ready summary tables

**Outputs**:
- Count matrices at 4 taxonomic levels (CSV)
- Statistical test results with FDR correction (CSV)
- NMDS ordination plots (PNG/PDF)
- Alpha diversity comparisons (boxplots, PNG/PDF)
- Taxonomic composition plots (barplots, PNG/PDF)
- SIMPER results for turbidity and depth (CSV)
- Manuscript-ready summary tables (CSV)

---

## Data Structure

### Input Data Format
The raw data (`CYCLE_ARMS_2mm_counts.csv`) contains the following key columns:
- **eventID**: ARMS deployment identifier (e.g., `CYCLE_2021_ARMS_05_DIAcoral`)
- **individualCount**: Number of individuals per taxonomic unit
- **kingdom**, **phylum**, **class**, **order**, **family**, **genus**, **scientificName**: Taxonomic hierarchy

### Taxonomic Resolution
Organisms are identified to the lowest possible taxonomic level:
- Molecular identification (COI barcoding) when possible
- Morphological identification to family, genus, or species
- Morphospecies designations when formal identification unavailable
- Symbol-based forward-filling preserves partial taxonomic information

---

## Data Processing Workflow

```
Raw Occurrence Data (individualCount per taxon)
         ↓
motile_organism_processing.py
         ├─> Taxonomic forward-filling
         ├─> Sample ID extraction
         ├─> Count matrix generation (species & order level)
         └─> Proportion calculations
         ↓
Motile_2mm_CYCLE21.Rmd
         ├─> Multi-level count tables (phylum, order, family, species)
         ├─> Alpha diversity calculations
         ├─> PERMANOVA & NMDS ordinations
         ├─> SIMPER environmental filtering analysis
         └─> Visualization & statistical outputs
```

---

## Environmental Variables

Environmental data integrated with community matrices:
- **turbidity_std_rank**: Standardized turbidity rank (continuous)
- **turbidity_m**: Mean estimated visibility (m)
- **depth**: Deployment depth (m)
- **temp_mean**: Mean temperature during deployment (°C)
- **temp**: In situ temperature at collection
- **salinity**: Salinity (PSU)
- **latitude** / **longitude**: Geographic coordinates

**Environmental Groupings**:
- **Turbidity Binary**: High vs. Low (median split)
- **Turbidity 3-Category**: Low, Medium, High (tertiles)
- **Depth Binary**: Shallow (<30m) vs. Deep (≥30m)
- **Depth 3-Category**: Shallow, Intermediate, Deep (tertiles)

---

## Site Localities

ARMS deployment sites across Gulf of Mexico banks:
- **Diaphus Bank**: Diaphus_coral, Diaphus_background
- **Alderdice Bank**: Alderdice_coral, Alderdice_background, Alderdice_shallow
- **McGrail Bank**: McGrail_coral
- **Bright Bank**: Bright_coral, Bright_background, Bright_shallow
- **Stetson Bank**: Stetson_coral
- **East Flower Garden Bank**: EFGB_shallow, EFGB_deep
updated naming: 
coral = deep1 
background = deep2
---

## Taxonomic Coverage

**Major Phyla Recovered**:
- **Arthropoda**: Amphipods, isopods, decapods, copepods, ostracods, etc.
- **Mollusca**: Gastropods, bivalves, polyplacophorans
- **Annelida**: Polychaetes, oligochaetes
- **Echinodermata**: Ophiuroids, echinoids, holothurians
- **Chordata**: Ascidians (tunicates)
- **Cnidaria**: Hydroids, corals
- **Bryozoa**: Colonial bryozoans
- **Porifera**: Small sponges and fragments
- **Others**: Nemerteans, platyhelminthes, minor phyla

**Typical Taxonomic Resolution**:
- Decapoda: Often to genus or species
- Amphipoda: Family to genus level
- Polychaeta: Family level (some to genus)
- Gastropoda: Family to genus level
- Ophiuroidea: Often to genus

---

## Statistical Approach

### Transformation
- **4th root transformation** applied to count data before ordination
- Reduces influence of highly abundant taxa
- Allows detection of patterns driven by rarer organisms
- Enables direct comparison with CoralNet sessile analyses

### Multiple Testing Correction
- **FDR (False Discovery Rate)** correction applied to all pairwise SIMPER comparisons
- Controls for family-wise error rate across multiple taxa and comparisons
- Thresholds: p < 0.05 (significant), p < 0.1 (marginally significant)

### Dissimilarity Metric
- **Bray-Curtis dissimilarity** for all ordinations
- Appropriate for count/abundance data
- Unaffected by double absences

---

## Requirements

**Python**:
- pandas
- numpy
- re, os

**R**:
- vegan (multivariate ecology)
- ggplot2, viridis, RColorBrewer (visualization)
- dplyr, tibble (data manipulation)
- gridExtra, cowplot (multi-panel plots)
- ape (phylogenetic/distance methods)
- reshape2 (data reshaping for plots)

---

## Key Analytical Comparisons

### Sessile vs. Motile Communities
This analysis is designed for direct comparison with CoralNet sessile community results:

| Aspect | CoralNet (Sessile) | Motile 2mm |
|--------|-------------------|------------|
| **Data type** | % cover from point counts | Individual counts |
| **Taxonomic level** | Phylum | Order (primary) |
| **Transformation** | 4th root | 4th root |
| **SIMPER analysis** | Environmental filtering | Environmental filtering |
| **Key question** | Which sessile taxa respond to turbidity/depth? | Do motile taxa show similar environmental responses? |


---

## Notes

- Individual counts represent total abundance per ARMS unit (all compartments combined)
- Taxonomic identifications combine molecular (COI) and morphological approaches
- Forward-filling symbols denote lowest available classification level
- Multiple taxonomic levels analyzed to balance resolution and data completeness
- SIMPER results identify taxa driving environmental differences
- 4th root transformation balances contributions of rare and common taxa
- Direct statistical comparison with sessile communities enables mechanistic inference

---

## Interpretation Guidelines

**Strong Turbidity Signal** → Environmental filtering during larval settlement or post-recruitment mortality

**Weak Turbidity Signal** → Behavioral dispersal ability of mobile organisms buffers against environmental gradients

**Depth Effects** → May reflect water mass characteristics, food availability, or habitat structure rather than depth per se

**Sessile-Motile Concordance** → Suggests environmental filtering at settlement phase (larvae)

**Sessile-Motile Discordance** → Suggests post-settlement processes (behavior, predation, competition) modify patterns

---

## Citation

If using these scripts, please cite:
Nicole C. Pittoors, Sarah M. Tweedt, Luke J. McCartin, Samuel A. Vohsen, Luisa Lopera, Sophia Mihalek, Jamie Lai, Kathleen Durkin, Lee Weigt, Marissa F. Nuttall, Annalisa Bracco, Christopher P. Meyer, Santiago Herrera. bioRxiv 2025.11.02.686126; doi: https://doi.org/10.1101/2025.11.02.686126

## Contact

Nicole Pittoors  
ncp220@lehigh.edu

---

*Last updated: December 2025*