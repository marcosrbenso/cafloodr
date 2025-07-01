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
