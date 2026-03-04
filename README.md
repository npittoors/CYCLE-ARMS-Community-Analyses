# Seascape Environmental Distance Matrix Pipeline
**Author:** Nicole Pittoors, Lehigh Oceans Research Center, Lehigh University  
**Manuscript:** *"Environmental filtering shapes patch dynamics across isolated mesophotic reefs"*  
**Script:** `seascape_env_model.Rmd`

---

## Overview

This R Markdown script constructs the environmental and geographic distance matrices used in downstream Mantel tests, partial Mantel tests, and multiple regression on distance matrices (MRM) to test environmental filtering versus dispersal limitation across mesophotic reef banks in the northwestern Gulf of Mexico.

The pipeline has three major components:

1. **Remote sensing data extraction** — Downloads and processes satellite-derived oceanographic variables from Copernicus Marine Service (CMEMS) NetCDF files and extracts site-level mean values
2. **Environmental distance matrices** — Builds scaled Euclidean environmental distance matrices for all samples, and separately for motile and sessile fractions
3. **Geographic distance matrices** — Calculates least-cost (bathymetry-constrained) geographic distances between sites using MARMAP

Output distance matrices feed directly into the community assembly analyses in `FINAL_18S_PreTaxa.Rmd`.

---

## Required R Packages

```r
library(ncdf4)           # NetCDF file handling
library(lubridate)       # Date/time conversion for NetCDF time variables
library(RColorBrewer)    # Color palettes
library(lattice)         # Visualization
library(dplyr)           # Data manipulation
library(oce)             # Oceanographic data processing
library(raster)          # Raster data handling and extraction
library(seacarb)         # Seawater carbonate chemistry (if needed)
library(rnaturalearth)   # Country/coastline polygons for maps
library(rnaturalearthdata)
library(marmap)          # Bathymetric data and least-cost distance calculation
```

---

## Input Files

| File | Description |
|---|---|
| `metadata_ARMS_env.csv` | ARMS sample metadata with site coordinates (`latitude`, `longitude`), `sitelocality`, `depth`, `turbidity_std_rank`, `fraction` |
| `EnvMatrix_ind_corrected.csv` | Environmental variable matrix at the individual ARMS level (rows = SH sample IDs); corrected for SH ID–environment value mismatches |
| `coord_sites_indv.csv` | Geographic coordinates for individual ARMS samples (rows = SH IDs) |
| `cmems_obs-oc_glo_bgc-optics_my_l4-multi-4km_P1M_*.nc` | CMEMS optical variables NetCDF (BBP, CDM) |
| `cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_*.nc` | CMEMS chlorophyll-a (CHL) NetCDF |
| `cmems_obs-oc_glo_bgc-pp_my_l4-multi-4km_P1M_*.nc` | CMEMS primary productivity (PP) NetCDF |
| `cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M_*.nc` | CMEMS transparency variables NetCDF (KD490, SPM, ZSD) |

> **Note on CMEMS files:** NetCDF files were downloaded from the [Copernicus Marine Service](https://marine.copernicus.eu/). The filenames above include timestamps from the original download; filenames may differ if re-downloaded. All variables are monthly composites at 4 km resolution.

---

## Environmental Variables Extracted

| Variable | Source | Description |
|---|---|---|
| `depth` | Metadata | Site depth (m) |
| `turbidity_std_rank` | Metadata | Standardized turbidity rank (benthic nepheloid layer proxy) |
| `PP_mean` | CMEMS PP | Mean primary productivity (mg C m⁻² d⁻¹) at site coordinates |
| `PP_sd` | CMEMS PP | Standard deviation of primary productivity over time |
| `CHL_mean` | CMEMS Plankton | Mean chlorophyll-a concentration (mg m⁻³) |
| `CHL_sd` | CMEMS Plankton | Standard deviation of chlorophyll-a over time |
| `BBP` | CMEMS Optics | Particulate backscattering coefficient at 443 nm |
| `CDM` | CMEMS Optics | Coloured dissolved and detrital organic matter |
| `KD490` | CMEMS Transp | Diffuse attenuation coefficient at 490 nm |
| `SPM` | CMEMS Transp | Suspended particulate matter (g m⁻³) |
| `ZSD` | CMEMS Transp | Secchi disk depth (m) |

---

## Script Structure

### Section 1 — Coordinate Setup

Derives site-level mean coordinates from `metadata_ARMS` for the 12 sampling sites (`coord_sites`) and extracts individual ARMS-level coordinates (`coord_ind`).

### Section 2 — Remote Sensing Data Processing

For each CMEMS NetCDF variable (PP, CHL, BBP/CDM optics, KD490/SPM/ZSD transparency):

1. Opens the NetCDF file with `ncdf4`
2. Converts the time variable from seconds-since-epoch to readable dates using `lubridate`
3. Converts the NetCDF to a `RasterBrick` object with correct spatial extent and WGS84 projection
4. Exports as GeoTIFF for archival
5. Loops through `coord_sites` to extract time-series data at the nearest raster cell to each site
6. Calculates mean and standard deviation across the time dimension
7. Appends results to `metadata_ARMS`

*Output files:* `PP_data.tif`, `CHL_data.tif`, and equivalent GeoTIFFs for each variable

### Section 3 — Environmental Matrix Construction

Assembles the full environmental matrix (`EnvMatrix_reduced`) from selected variables across the 12 sites and applies z-score standardization (`scale(center = T, scale = T)`). Includes:

- Boxplots before and after standardization for QC
- Heatmap visualization of the full environmental matrix and a subset of key variables (`turbidity_std_rank`, `depth`, `PP_mean`)
- PCA (`prcomp()`) to explore multicollinearity and variance structure; PC loadings and variance explained plotted for PC1 and PC2

*Output files:* `all_variables_heatmap.pdf`, `selected_variables_heatmap.pdf`

### Section 4 — Environmental Distance Matrices

Loads the individual-level environmental matrix (`EnvMatrix_ind_corrected.csv`; rows = SH IDs, corrected August 2025) and computes Euclidean distance matrices on z-score standardized variables.

**Three distance matrices calculated:**

| Object | Samples | Description |
|---|---|---|
| `env_dist` / `env_dist_df` | All 114 ARMS | Full environmental distance matrix |
| `env_dist_sess` | 36 sessile samples | Subset for sessile fraction Mantel tests |
| `env_dist_mot` | 78 motile samples | Subset for motile fraction Mantel tests |

Individual distance matrices for key variables also computed separately (`depth_dist`, `turb_dist`, `PP_dist`) for univariate comparisons.

**Fraction subsets** are defined by hardcoded SH ID vectors:
- `sessile_rows` — 36 SH IDs (sessile fraction)
- `motile_rows` — 78 SH IDs (100 µm + 500 µm motile fractions)

*Output files:* `env_dist_matrix.csv`, `env_dist_input.Rdata`

### Section 5 — Geographic Distance Matrices (MARMAP)

Calculates least-cost geographic distances constrained by bathymetry using `marmap`, ensuring distances follow navigable water paths rather than straight-line overwater distances.

**Steps:**
1. Downloads NOAA bathymetric data for the study region (lon: −95 to −90, lat: 27 to 30; resolution: 0.2 arc-min) via `getNOAA.bathy()`
2. Creates a transition matrix (`trans.mat()`) with `min.depth = -0.1` to constrain paths to water
3. Calculates least-cost distance matrices for sites and individual ARMS using `lc.dist()`

**Six distance matrices calculated:**

| Object | Samples | Output file |
|---|---|---|
| `leastDist.km` | 12 sites | `leastDist_km.csv` |
| `leastDist.km_ind` | All 114 ARMS | `leastDist_km_ind.csv` |
| `leastDist.km_mot` | 78 motile ARMS | `leastDist_km_motile.csv` |
| `leastDist.km_sess` | 36 sessile ARMS | `leastDist_km_sessile.csv` |

> **Note:** Minimum inter-site distances within the same bank approach zero (sites positioned very close together within a bank), which is expected given the study design.

---

## Output Files

| File | Description | Used In |
|---|---|---|
| `env_dist_matrix.csv` | Full 114×114 environmental distance matrix | `FINAL_18S_PreTaxa.Rmd` Mantel/MRM |
| `env_dist_input.Rdata` | R objects: `sEnvMatrix_ind`, `env_dist_df`, `env_dist_mot`, `env_dist_sess` | `FINAL_18S_PreTaxa.Rmd` |
| `leastDist_km.csv` | 12×12 least-cost geographic distance matrix (sites) | Visualization |
| `leastDist_km_ind.csv` | 114×114 least-cost geographic distance matrix (all ARMS) | `FINAL_18S_PreTaxa.Rmd` Mantel/MRM |
| `leastDist_km_motile.csv` | 78×78 geographic distance matrix (motile) | `FINAL_18S_PreTaxa.Rmd` |
| `leastDist_km_sessile.csv` | 36×36 geographic distance matrix (sessile) | `FINAL_18S_PreTaxa.Rmd` |
| `leastDist.km_ind.Rdata` | R dist object, all ARMS geographic distances | `FINAL_18S_PreTaxa.Rmd` |
| `leastDist.km_mot.Rdata` | R dist object, motile geographic distances | `FINAL_18S_PreTaxa.Rmd` |
| `leastDist.km_sess.Rdata` | R dist object, sessile geographic distances | `FINAL_18S_PreTaxa.Rmd` |
| `PP_data.tif`, `CHL_data.tif`, etc. | GeoTIFF rasters for each remote sensing variable | Archival |
| `all_variables_heatmap.pdf` | Heatmap of full environmental matrix | QC / visualization |
| `selected_variables_heatmap.pdf` | Heatmap of key variables (turbidity, depth, PP) | QC / visualization |

---

## Design Notes

- **Z-score standardization is required** before computing Euclidean distances because environmental variables span different units and scales (e.g., depth in meters vs. turbidity rank 1–13 vs. PP in mg C m⁻² d⁻¹). Failing to standardize would allow variables with larger absolute values to dominate the distance metric.
- **Least-cost rather than straight-line geographic distances** are used to respect the bathymetric structure of the Texas-Louisiana shelf, where shallow muddy substrate between banks may constrain larval dispersal pathways.
- **`EnvMatrix_ind_corrected.csv` should be used** (not the uncorrected version) — SH ID to environmental value mismatches were corrected in August 2025.
- **Fraction-specific matrices are necessary** because Mantel and MRM tests must use distance matrices with matching dimensions to the Bray-Curtis community dissimilarity matrices computed in `FINAL_18S_PreTaxa.Rmd`.

---

## Usage Notes

Set the working directory at the top of the script to the folder containing your environmental data and CMEMS NetCDF files:

```r
setwd("path/to/env_data")
```

CMEMS NetCDF files are large and the raster extraction loop can be slow for many sites. If re-running, load previously exported GeoTIFFs directly with `raster::brick("PP_data.tif")` rather than re-processing the NetCDF files.

The `marmap::getNOAA.bathy()` call requires an internet connection to download bathymetric data from NOAA. The downloaded `bathydata` object and transition matrix `t` can be saved and reloaded to avoid repeated downloads:

```r
save(bathydata, t, file = "bathy_transition.Rdata")
load("bathy_transition.Rdata")
```

---

**references:**
- MARMAP: Pante E & Simon-Bouhet B (2013) marmap: A Package for Importing, Plotting and Analyzing Bathymetric and Topographic Data in R. *PLoS ONE* 8(9):e73051
- CMEMS: Copernicus Marine Service, https://marine.copernicus.eu/
