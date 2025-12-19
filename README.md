# FATGBM: Fractal Activity Time Geometric Brownian Motion

This repository contains the R implementation of the **FATGBM model**,  and option pricing methods presented in the research paper:

> **[Insert Paper Title Here]** > *Authors: [Insert Author Names]* > *Journal/Conference: [Insert Journal Name, Year]* > [Link to paper if available (e.g., arXiv or DOI)]

## Overview

The FATGBM model generalizes standard Geometric Brownian Motion by incorporating **Fractal Activity Time** to better capture the leptokurtic nature and volatility clustering observed in financial markets. 

This repository provides:
1.  **Stochastic Process Simulation:** Algorithms to simulate the **superposition of Ornstein-Uhlenbeck (supOU)** process, which is used to simulate the fractal activity time 
2.  **FATGBM Model:** Core functions to simulate the Fractal Activity Time Geometric Brownian Motion process.
3.  **Option Pricing:** Numerical methods for pricing barrier options under FATGBM.
4.  **Reproduction Scripts:** Code to reproduce the figures and tables presented in the paper.

## Repository Structure

### Core Models
* `FATGBM.R`: Main script containing the FATGBM model definition and option pricing logic.
* `RGamma_supOU_s.R`: Implementation of the Random Gamma process / Superposed Ornstein-Uhlenbeck (supOU) process used for activity time modeling.

### Figures (Reproduction)
The following scripts generate the figures found in the manuscript:
* `fig123.R`: Generates Figures 1, 2, and 3 (Trajectories and densities).
* `fig5.R` - `fig7.R`: Generates Figures 5 through 7.
* `fig8-9.R`: Generates Figures 8 and 9.
* `fig10.R` & `fig11.R`: Generates volatility surfaces and calibration plots.

### Tables (Pricing & Analysis)
* `table2.R`: Code to reproduce Table 2 (Parameter estimation/comparison).
* `table3pricing.R`: Code to reproduce Table 3 (Option pricing results).

## Prerequisites

To run these scripts, you need **R** installed.
You may also need the following R packages (install via `install.packages("package_name")`):

* `stats`
* `MASS`
* *[List any other libraries you used, e.g., 'numDeriv', 'ggplot2', etc.]*

## Usage

### 1. Simulating the Model
To use the core functions, source the main files at the top of your script:

```r
source("FATGBM.R")
source("RGamma_supOU_s.R")

# Example function call (replace with actual example from your code)
# price <- FATGBM_pricing(S0=100, K=100, T=1, ...)
