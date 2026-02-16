# Evaluation of Sentinel-2 Vegetation Indices for Estimating Leaf Area Index in Cassava Plots

This repository contains the data and code accompanying the manuscript:

## Overview

This study systematically evaluates 13 Sentinel-2-derived vegetation indices (VIs) for estimating ground-measured Leaf Area Index (LAI) in cassava (*Manihot esculenta* Crantz) plots across different growth stages. Ground-LAI was measured monthly using a SunScan Canopy Analyzer from January to June 2022 (2–7 months after planting) in 47 cassava plots in Nakhon Ratchasima Province, Thailand.

## Repository Structure

```
├── DATA/                                    # Input datasets
├── GEE-script.js                            # Google Earth Engine script for Sentinel-2
│                                            #   image acquisition, preprocessing, and
│                                            #   vegetation index computation
├── Zonal-statistic.R                        # Extraction of plot-level VI values using
│                                            #   zonal statistics (median aggregation)
├── VIs-time-series.R                        # Temporal analysis and visualization of
│                                            #   VIs and ground-LAI across growth stages
├── Linear-mixed-effects-model.R             # Linear mixed-effects models evaluating VI
│                                            #   performance across all growth stages
├── Simple-linear-regression-model.R         # Stage-specific linear regression models
│                                            #   (per MAP) for each VI
├── Bootstraping-Mixed-Effects-Models.R      # Bootstrap validation comparing model
│                                            #   performance between with-tree and
│                                            #   without-tree plots
└── README.md
```
## Vegetation Indices

The 13 Sentinel-2-derived vegetation indices evaluated:

| Index | Name | Bands Used |
|-------|------|------------|
| NDVI | Normalized Difference Vegetation Index | NIR, Red |
| SAVI | Soil-Adjusted Vegetation Index | NIR, Red |
| EVI | Enhanced Vegetation Index | NIR, Red, Blue |
| BNDVI | Blue Normalized Difference Vegetation Index | NIR, Blue |
| CIG | Chlorophyll Index Green | NIR, Green |
| DVI | Difference Vegetation Index | NIR, Red |
| GNDVI | Green Normalized Difference Vegetation Index | NIR, Green |
| GRVI | Green Ratio Vegetation Index | NIR, Green |
| NDWI | Normalized Difference Water Index | Green, NIR |
| RVI | Ratio Vegetation Index | Red Edge 2, Red |
| SeLI | Sentinel-2 LAI Green Index | NIR2, Red Edge 1 |
| TCARI | Transformed Chlorophyll Absorption in Reflectance Index | Green, Red, Red Edge 1 |
| VIG | Vegetation Index Green | Green, Red |

## Study Area

- **Location**: Nong Bua Sa-Art Sub-district, Bua Yai District, Nakhon Ratchasima Province, Thailand (102.262°E, 15.575°N)
- **Plots**: 47 cassava plots (~0.71 ha total area)
- **Growing season**: 2021/2022 late rainy season (planting: Dec 2021 – Jan 2022)
- **LAI measurement period**: January – June 2022 (2–7 MAP)

## How to Use

1. **Run the GEE script** (`GEE-script.js`) in the [Google Earth Engine Code Editor](https://code.earthengine.google.com/) to acquire and process Sentinel-2 imagery and compute vegetation indices. Export the results as GeoTIFF files.
2. **Extract plot-level VI values** by running `Zonal-statistic.R` with the exported GeoTIFFs and plot boundary shapefiles.
3. **Run the analysis scripts** in R:
   - `VIs-time-series.R` for temporal pattern visualization
   - `Linear-mixed-effects-model.R` for overall VI performance evaluation
   - `Simple-linear-regression-model.R` for stage-specific analysis
   - `Bootstraping-Mixed-Effects-Models.R` for spatial heterogeneity assessment

## Contact

if you have inquiries about the code and data, please contact;
- Kanokporn Promnikorn — kanokporn.promn@ku.th
- Ekaphan Kraichak (corresponding author) — ekaphan.k@ku.th

Department of Botany, Faculty of Science, Kasetsart University, Bangkok, Thailand
