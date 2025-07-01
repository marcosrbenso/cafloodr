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
  - `caflood.exe` (the core hydrodynamic simulation engine) that must be saved 
  - Example input files and documentation

---

## Input data

### DEM

The **Digital Elevation Model (DEM)** provides the terrain surface for hydrological and hydrodynamic calculations.

**Requirements:**

- **Format**: GeoTIFF (`.tif`) or ASCII Grid (`.asc`)
- **Coordinate Reference System**: Projected in **metric units** (e.g., UTM)
- **Units**: Elevation values in **meters**
- **Resolution**: High-resolution (1–10 meters) recommended for local-scale or urban studies

Use the `terra::rast()` or `raster::raster()` functions to load and inspect your DEM in R:

```r
library(terra)
dem <- rast("path/to/dem.tif")
plot(dem)
```
### Outlets
Outlet points define where the model extracts simulated hydrographs (e.g., water depth and velocity over time).

Format Requirements

- **File format**: ESRI Shapefile (`.shp`) with accompanying `.dbf`, `.shx`, and `.prj` files.
- **Geometry type**: `POINT`
- **Coordinate Reference System**: Must match the DEM's CRS (typically a UTM projection in meters).
- **Attribute field**: Must include a numeric column named `ID` to uniquely identify each outlet.

Example

| ID   | geometry (X, Y)       |
|------|------------------------|
| 11   | POINT (450300, 7456300)|
| 157  | POINT (451000, 7455000)|
| 1001 | POINT (452000, 7454200)|


```r
library(sf)

outlets <- st_read("data/outlets.shp")
head(outlets)

# Check if 'ID' exists
if (!"ID" %in% names(outlets)) {
  stop("Outlet shapefile must contain an 'ID' column.")
}

# Visualize
plot(outlets["ID"], main = "Outlet Points")
```

### Precipitation file

Rainfall drives the flood simulation and must be provided as a CSV file with regular time steps.

Format Requirements
- File format: `.csv`
- Columns:
  - `Row`
  - `Start`: Start time of rainfall step (YYYY-MM-DD HH:MM:SS)
  - `End`: End time of rainfall step
  - `Events`: Name of the rainfall event
  - `Prec`: Rainfall depth/intensity in millimeters (mm/10min)


| Row | Start               | End                 | Events     | Prec               |
|-----|---------------------|---------------------|------------|--------------------|
| 1   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |
| 2   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |
| 3   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |
| 4   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |
| 5   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0.00749712655111925|
| 6   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |
| 7   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |
| 8   | 2017-04-06 20:00:00 | 2017-04-07 18:00:00 | Event_1_SP | 0                  |

```r
rain <- read.csv("data/rain_event.csv")

# Check structure
str(rain)

# Optional: Impute missing values
library(imputeTS)
rain$prec <- na_kalman(rain$prec)

# Visualize time series
plot(rain$start, rain$prec, type = "l", col = "blue",
     xlab = "Time", ylab = "Rainfall (mm)", main = "Rainfall Time Series")
```

## Run CADDIES/CAFLOOD in R

After downloading:

```bash
# Unzip and copy caflood.exe to a folder, e.g.:
C:/Models/CAFlood/caflood.exe

## Run CAflood in R
```

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



```
