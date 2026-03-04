# CYCLE-ARMS Community Analyses

**Authors:** Nicole C. Pittoors, Sarah M. Tweedt, Luke J. McCartin, Samuel A. Vohsen, Luisa Lopera, Sophia Mihalek, Jamie Lai, Kathleen Durkin, Lee Weigt, Marissa F. Nuttall, Annalisa Bracco, Christopher P. Meyer, Santiago Herrera  
**Contact:** Nicole Pittoors — ncp220@lehigh.edu  
**Corresponding author:** Santiago Herrera — santiago.herrera@alum.mit.edu

**Manuscript:** *"Environmental filtering shapes patch dynamics across isolated mesophotic reefs"*  
**Preprint:** bioRxiv 2025.11.02.686126 — https://doi.org/10.1101/2025.11.02.686126

---

## Overview

This repository contains all code and associated data files for community analyses of cryptobenthic invertebrate assemblages across six mesophotic reef banks on the Texas-Louisiana continental shelf, Gulf of Mexico. Samples were collected using Autonomous Reef Monitoring Structures (ARMS) during the CYCLE 2021 expedition.

The study integrates six complementary biodiversity datasets to test whether environmental filtering or dispersal limitation governs community assembly across geographically isolated mesophotic reef patches:

| Dataset | Method | Primary question |
|---|---|---|
| COI metabarcoding (ARMS) | DNA metabarcoding, Folmer/Leray primers | Community composition, beta diversity, dispersal vs. filtering |
| 18S metabarcoding (ARMS) | DNA metabarcoding, V4_18SNext primers | Community composition, beta diversity, dispersal vs. filtering |
| 18S eDNA (water column) | DNA metabarcoding of filtered seawater, V4_18SNext primers | Water column community composition, environmental filtering |
| CoralNet image annotations | Point-count annotation of ARMS settlement plates | Sessile macrofaunal cover and recruitment |
| 2mm motile macrofauna | Morphological + molecular identification of sieved fauna | Motile invertebrate community structure |
| Environmental / seascape model | Satellite remote sensing + HOBO logger data | Environmental distance matrices for Mantel/MRM analyses |

---

## Repository Structure

```
CYCLE-ARMS-Community-Analyses/
├── ARMS_Metabarcoding/          # COI and 18S metabarcoding pipeline (DADA2 → diversity analyses)
├── ARMS_Assign_Taxonomy_Pipeline/  # Custom BLAST-based taxonomy assignment pipeline
├── Environmental_Model/         # Environmental and geographic distance matrix construction
├── eDNA_18S/                    # 18S eDNA water column metabarcoding analyses
├── CoralNet_ARMS/               # CoralNet image annotation processing and community analyses
├── 2mm_motile_vouchers/         # 2mm sieved motile fauna processing and community analyses
└── README.md                    # This file
```

---

## Directory Descriptions

### [`ARMS_Metabarcoding/`](./ARMS_Metabarcoding)

Full bioinformatics pipeline for COI and 18S rRNA metabarcoding data. Covers quality control (FastQC/MultiQC), primer trimming (Cutadapt), ASV inference (DADA2), LULU curation, and all downstream diversity analyses, including alpha diversity (LMMs), beta diversity (NMDS, PERMANOVA, dbRDA), distance-based community assembly tests (Mantel, partial Mantel, MRM), differential abundance (ANCOM-BC2), and functional trait analyses.

**Scripts:** `QC_trim_DADA2.md`, `FINAL_18S_PreTaxa.Rmd`, `PostTaxa_18S.Rmd`, `lulu_curation.R`  
**See:** [`ARMS_Metabarcoding/README.md`](./ARMS_Metabarcoding/README.md)

---

### [`ARMS_Assign_Taxonomy_Pipeline/`](./ARMS_Assign_Taxonomy_Pipeline)

Custom hierarchical BLAST pipeline for taxonomic assignment of COI and 18S ASVs against curated marine reference databases (NMNH vouchers → MetaZooGene → NCBI). Includes taxonomy backfilling via WoRMS or GBIF, marine habitat verification, and filtering to retain Metazoa and Rhodophyta.

**Scripts:** `taxonomy_master.sh`, `MZG_BLAST.py`, `nmnh_taxonomy_processor.py`, `best_taxonomy_selector.py`, `final_process_taxonomy.py`  
**See:** [`ARMS_Assign_Taxonomy_Pipeline/README.md`](./ARMS_Assign_Taxonomy_Pipeline/README.md)

---

### [`Environmental_Model/`](./Environmental_Model)

Constructs environmental and geographic distance matrices used in downstream Mantel and MRM analyses. Extracts satellite-derived oceanographic variables (primary productivity, chlorophyll-a, turbidity proxies, SPM) from Copernicus Marine Service (CMEMS) NetCDF files, builds scaled Euclidean environmental distance matrices, and calculates bathymetry-constrained least-cost geographic distances using MARMAP.

**Scripts:** `seascape_env_model.Rmd`  
**See:** [`Environmental_Model/README.md`](./Environmental_Model/README.md)

---

### [`eDNA_18S/`](./eDNA_18S)

Analysis pipeline for 18S rRNA metabarcoding of water column environmental DNA (eDNA) collected at 12 site-depth combinations across six reef banks (33 samples, 830 ASVs). Covers quality control, metadata harmonization (site label mapping, environmental variable integration), rarefaction, alpha diversity (Shannon, Simpson, Pielou's evenness, Chao1, ACE), community overlap and ASV sharing analyses, beta diversity (NMDS, PERMANOVA, dbRDA), and ANCOM-BC II differential abundance testing against turbidity, depth, and primary productivity gradients at phylum and family levels. Site color scheme and ordering are harmonized with ARMS metabarcoding analyses for cross-dataset comparisons.

**Scripts:** `eDNA18S_pub.Rmd`  
**See:** [`eDNA_18S/README.md`](./eDNA_18S/README.md)

---

### [`CoralNet_ARMS/`](./CoralNet_ARMS)

Processes CoralNet point-count annotations from ARMS settlement plate images to quantify percent cover and recruitment of benthic taxa at the phylum level. Integrates environmental metadata and performs multivariate community analyses (PERMANOVA, NMDS, dbRDA, SIMPER) to assess environmental filtering of sessile macrofaunal communities.

**Scripts:** `CoralNet_abundances.py`, `CoralNet_dataprep_postMETA.py`, `CN_proportion_check.py`, `CoralNet_CYCLE21_All_Plates.Rmd`  
**See:** [`CoralNet_ARMS/README.md`](./CoralNet_ARMS/README.md)

---

### [`2mm_motile_vouchers/`](./2mm_motile_vouchers)

Processes taxonomic count data for motile macrofauna (≥2mm) sieved from ARMS units. Organisms were identified morphologically and/or by COI barcoding. Pipeline converts occurrence records to community matrices at multiple taxonomic levels, calculates alpha and beta diversity, and performs SIMPER analyses to assess turbidity and depth effects — enabling direct comparison with sessile CoralNet communities to distinguish environmental from behavioral drivers of community structure.

**Scripts:** `motile_organism_processing.py`, `Motile_2mm_CYCLE21.Rmd`  
**See:** [`2mm_motile_vouchers/README.md`](./2mm_motile_vouchers/README.md)

---

## Data Availability

| Data | File | Location |
|---|---|---|
| Raw sequences (COI + 18S) | — | NCBI SRA, BioProject [PRJNA1159220](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1159220) |
| 18S ASV table | `18S_metazoo_ASVtab.csv` | `ARMS_Metabarcoding/Data/` |
| COI ASV table | `COI_metazoan_ASVtab.csv` | `ARMS_Metabarcoding/Data/` |
| 18S taxonomy table | `taxa_metazoan_18S.tsv` | `ARMS_Metabarcoding/Data/` |
| COI taxonomy table | `taxa_metazoan_COI_final.csv` | `ARMS_Metabarcoding/Data/` |
| 18S metadata | `metadata_18S_FINAL.csv` | `ARMS_Metabarcoding/Data/` |
| COI metadata / environmental variables | `metadata_ARMS_env_COI.csv` | `ARMS_Metabarcoding/Data/` |
| Functional trait mapping table | `trait_mapping_QC_v4.csv` | `ARMS_Metabarcoding/Data/` |
| eDNA 18S phyloseq object | `eDNA18S_phyloseq_20250915.rds` | `eDNA_18S/` |
| eDNA 18S cleaned data + site mapping | `eDNA_18S_cleaned_data_with_mapping.Rdata` | `eDNA_18S/` |
| CoralNet annotations | `ALL_ARMS_CoralNet_Annotations.csv` | `CoralNet_ARMS/` |
| 2mm motile counts | `CYCLE_ARMS_2mm_counts.csv` | `2mm_motile_vouchers/` |

---

## Complete Analysis Workflow

```
Raw sequences (NCBI SRA: PRJNA1159220)
        ↓
[ARMS_Metabarcoding] QC_trim_DADA2.md
FastQC/MultiQC → Cutadapt → DADA2 → LULU
        ↓
[ARMS_Assign_Taxonomy_Pipeline] taxonomy_master.sh
BLAST (NMNH → MetaZooGene → NCBI) → WoRMS/GBIF backfill
        ↓
[Environmental_Model] seascape_env_model.Rmd
CMEMS satellite data → Environmental distance matrices
MARMAP bathymetry → Geographic distance matrices
        ↓
[ARMS_Metabarcoding] FINAL_18S_PreTaxa.Rmd
Alpha diversity (LMMs) → Beta diversity (NMDS, PERMANOVA, dbRDA)
Mantel / partial Mantel / MRM → Environmental filtering vs. dispersal
        ↓
[ARMS_Metabarcoding] PostTaxa_18S.Rmd
Relative abundance → ANCOM-BC2 differential abundance
Functional trait analysis
        ↓
[eDNA_18S] eDNA18S_pub.Rmd                          [CoralNet_ARMS] + [2mm_motile_vouchers]
Water column 18S eDNA (33 samples, 830 ASVs)         Sessile cover (CoralNet) + Motile fauna (2mm)
Alpha/beta diversity → ANCOM-BC2 env. filtering      PERMANOVA / NMDS / SIMPER
        ↓                                                        ↓
                    Cross-dataset comparison of environmental filtering signals
```

---

## Citation

If using any scripts or data from this repository, please cite:

> Pittoors NC, Tweedt SM, McCartin LJ, Vohsen SA, Lopera L, Mihalek S, Lai J, Durkin K, Weigt L, Nuttall MF, Bracco A, Meyer CP, Herrera SP. *Environmental filtering shapes patch dynamics across isolated mesophotic reefs.* bioRxiv 2025.11.02.686126. https://doi.org/10.1101/2025.11.02.686126

