<!-- README.md is generated from README.Rmd. Please edit that file -->

# cafloodr <img src="https://img.shields.io/badge/R-CAFlood-blue.svg" align="right" />

<!-- badges: start -->
<!-- You can add badges like CRAN, GitHub Actions, etc., if desired -->
<!-- badges: end -->

## Overview

**`cafloodr`** is an R package that provides functions to preprocess Digital Elevation Models (DEMs) and run the CADDIES/CAFLOOD 2D hydrodynamic flood model.  
It is designed for use on **Windows operating systems** and integrates R-based tools with Python and compiled CAFLOOD simulations.  

Core features include:

- Automated preparation of rainfall input, DEM carving, and flow direction.
- Generation of control files and outlet configurations.
- Execution of the CAFLOOD model directly from R.
- Support for Python model integration (e.g., Forward-Mole 1D).

---

## Installation

To install the development version of `cafloodr` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("marcosrbenso/cafloodr")
```

## Install CADDIES/CAFLOOD

To use **`cafloodr`**, you need the standalone **CAFlood** binary (Windows):

### How to Download

1. Go to the **CAFlood Software** page:  
   [https://cafloodpro.com/caflood-software/](https://cafloodpro.com/caflood-software/) :contentReference[oaicite:1]{index=1}

2. Click **"DOWNLOAD CAFLOOD FREE"** under the free version.

3. Complete the form (name, email, organization, etc.).  
   After submission, you'll receive a download link for the Windows executable.:contentReference[oaicite:2]{index=2}

---

### What You’ll Get

- A Windows `.zip` or `.exe` package containing:
  - `caflood.exe` (the core hydrodynamic simulation engine)
  - Example input files and documentation (e.g., `CADDIES-manual-caflood-110.pdf`):contentReference[oaicite:3]{index=3}

---

## Run CADDIES/CAFLOOD in R

After downloading:

```bash
# Unzip and copy caflood.exe to a folder, e.g.:
C:/Models/CAFlood/caflood.exe

## Run CAflood in R

```r

library(cafloodr)

# Define model parameters
param <- list(
  alpha_fraction = 0.1,
  roughness_global = 0.035,
  slope_tolerance = 0.5,
  boundary_ele = 580,
  Raster_WD_Tolerance = 0.01,
  Upstream_Reduction = 0.5
)

# Define paths
caflood_path <- "C:/Models/CAFlood/caflood.exe"
Input <- "C:/Projects/FloodSim/Input/"
Output <- "C:/Projects/FloodSim/Output/"
event_name <- "storm_2025_01"
dem_name <- "dem_sp"
path2 <- "DEMs/"
outlet_path <- "C:/Projects/FloodSim/Shapefiles/outlets.shp"
rain_path <- "C:/Projects/FloodSim/Rain/rain_event.csv"
stream_network <- "C:/Projects/FloodSim/Streams/stream"
snap_dist <- 30  # max distance (m) to snap outlets to stream

# Run the CAFlood simulation
cafloodr_A01(
  param = param,
  caflood_path = caflood_path,
  Input = Input,
  Output = Output,
  event_name = event_name,
  dem_name = dem_name,
  path2 = path2,
  outlet_path = outlet_path,
  rain_path = rain_path,
  stream_network = stream_network,
  snap_dist = snap_dist
)





```
