# FATGBM: Fractal Activity Time Geometric Brownian Motion

This repository contains the R implementation of the **FATGBM model**,  and option pricing methods presented in the research paper:

> **[A risky asset Student model with long-range dependence through fractal activity time]** > *Authors: [Nikolai N. Leonenko, Andrey Pepelyshev and Bader Saidan]* > 

## Overview

The FATGBM is a stochastic model that generalizes the classical geometric Brownian motion (GBM) by accommodating a number of properties that are not captured by Brownian motion but which are observed in the data, such as dependence in absolute and squared returns (but not returns themselves) and a marginal distribution that is heavier tailed and higher peaked than Gaussian,

This repository provides:
1.  **Stochastic Process Simulation:** Algorithms to simulate the **superposition of Ornstein-Uhlenbeck (supOU)** process with inverse gamma marginals. This algorithms is used to simulate the fractal activity time. 
2.  **FATGBM Model:** Core functions to simulate the fractal activity time geometric Brownian motion process.
3.  **Option Pricing:** Numerical methods for pricing barrier options under FATGBM.
4.  **Reproduction Scripts:** Code to reproduce the figures and tables presented in the paper.

## Repository Structure

### Core Models
* `RGamma_supOU_s.R`: Implementation of supOU process with inverse gamma marginals used for activity time modeling.
* * `FATGBM.R`: Implementation of FATGBM model and option pricing.

### Figures (Reproduction)
The following scripts generate the figures found in the manuscript:
* `fig123.R`: Generates Figures 1, 2, and 3 .
* `fig5.R` - `fig7.R`: Generates Figures 5 through 7.
* `fig8-9.R`: Generates Figures 8 and 9.
* `fig10.R` & `fig11.R`: Generates Figure 10 and 11.

### Tables (Pricing & Analysis)
* `table2.R`: Code to reproduce Table 2 (Parameter estimation /comparison).
* `table3pricing.R`: Code to reproduce Table 3 (Option pricing comparison).

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
