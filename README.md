<!-- README.md is generated from README.Rmd. Please edit that file -->

# Outline

[cafloodr](#cafloodr) Introduction to CAFLOOD in R 

[vignettes](#vignettes) Full workflow on how to run CAFLOOD in R

# cafloodr: Caddies/Caflood for R

<img src="https://img.shields.io/badge/R-CAFlood-blue.svg" align="right" />

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
   [https://cafloodpro.com/caflood-software/](https://cafloodpro.com/caflood-software/)

2. Click **"DOWNLOAD CAFLOOD FREE"** under the free version.

3. Complete the form (name, email, organization, etc.).  
   After submission, you'll receive a download link for the Windows executable.

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


# Vignettes

## Import libraries

```
library(imputeTS)
library(zoo)
library(gstat)
library(sf)
library(sp)
library(cafloodr)

```

## Interpolate precipitation data

We start by preparing spatially distributed rainfall input via **Inverse Distance Weighting (IDW)**. This example processes data from multiple rainfall stations, interpolates over a catchment shapefile (aricanduva), and writes a precipitation file in CSV format compatible with cafloodr.

```


dataset <- read.csv("Data/rainfall_2021.csv")
dataset

stations <- (read.csv("Data/CoordenadasEstacoes (coordinate).csv"))
coordinates(stations) <- ~X+Y

aricanduva <- st_read("Data/aricanduva.shp")

# IDW

events <-
  data.frame(
    id = c("Event_1_SP"),
    start = as.POSIXct(c("2021-01-25 20:00"),tz = "UTC"),
    end =  as.POSIXct(c("2021-02-19 18:00"),tz = "UTC")
  )

rain_events <- c()

for(t in 1:nrow(events)){
  dates <- seq(from = (events$start[t]),
               to = (events$end[t]),by = "10 min")


  results <- rep(NA,length(dates))
  for(i in 1:length(dates)){
    tryCatch({
      sub_sel <- dataset |>
        subset(DATA == dates[i]) |>
        merge(stations,
              by.x = "Posto",
              by.y = "ID") |>
        dplyr::select(X,Y,Rain3) |>
        na.omit() |>
        st_as_sf(coords = c("X",'Y'),crs="WGS84") |>
        st_transform(crs = st_crs(aricanduva))


      rainSP <- as(sub_sel, 'Spatial')
      neighbors <- length(rainSP) - 1
      beta <- 2

      idw_rain <- gstat::gstat(
        formula = Rain3 ~ 1, # intercept-only model
        data = rainSP,
        nmax = neighbors,
        set = list(idp = beta)
      )

      predict(
        idw_rain, aricanduva
      )$var1.pred -> results[i]

    },error = function(e){
      NA
    })

  }

  rain_events[[t]] <- unlist(results)


  write.csv(
    data.frame(start = events$start[t],
               end = events$end[t],
               events = events$id[t],
               prec = unlist(results)),
    paste0(events$id[t],"prec.csv")
  )
}


```

## Setting up the working paths

Define paths to the CAFLOOD executable, input folder, and output directory.

```



caflood_path <- "C:/Path/CAFLOOD"  # Path to caflood.exe
Input <- "Data/"                   # Folder containing model inputs
Output <- "E:/flood_model_results/"  # Folder to store results


```

## Running CAFLOOD


### Set the parameters
Set up the key simulation parameters using a named list. These parameters control the CAFLOOD simulation engine.

```{r, fig.show='hold'}

param <- list(
  alpha_fraction = 0.1,
  roughness_global = 0.009,       # Manning’s n (lower = faster flow)
  slope_tolerance = 0.528,
  boundary_ele = -9,              # Default boundary elevation
  Raster_WD_Tolerance = 0.1,
  Upstream_Reduction = 0.5
)


```

### Run CAFLOOD

Call the `cafloodr_A01()` wrapper to run a full simulation, including:

- Preprocessing DEM and rainfall
- Creating control files
- Snapping outlet points
- Executing the flood simulation

```

cafloodr_A01(
  param = param,
  caflood_path = caflood_path,
  Input = Input,
  Output = Output,
  event_name = "Event_test",
  dem_name = "ari_lidar_30m",  # name of DEM (without extension)
  path2 = "",
  outlet_path = "estacoes.shp",
  rain_path = "Event_1_SP_prec.csv",
  stream_network = "aricanduva_streams2",  # name (without extension)
  snap_dist = 30
)

```
