# Species Distribution Modeling (SDM) Pipeline

This repository contains the R pipeline to run a complete **Species Distribution Modeling (SDM)** workflow using multiple algorithms, spatial thinning, multicollinearity analysis, and ensemble strategies (*Ensemble Modeling* and *Ensemble of Small Models - ESM*).

The script is structured by default using ***Nasua nasua*** as an example species, but it can easily be adapted for any target species.

---

## Directory Structure

To ensure the script runs smoothly, make sure your local directory matches the following structure:

```text
.
├── 1_OCCURENCES/
│   ├── species_occ.shp      # Species occurrence points (Shapefile)
│   └── study_area.shp       # Study area boundary (Shapefile)
├── ENV_VARIABLES/
│   └── *.tif                # Environmental variable raster files
├── 4_RESULTS/               # Output directory for generated tables and figures
├── scripts/
│   └── sdm_script.R         # Main R script
└── README.md                # Project documentation
```
