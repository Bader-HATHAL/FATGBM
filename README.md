# FATGBM: Fractal Activity Time Geometric Brownian Motion

This repository contains the R implementation of the **FATGBM model** and the option pricing methods presented in the following research papers and PhD thesis:

> 1. **A risky asset Student model with long-range dependence through fractal activity time**
>    *Authors: Nikolai N. Leonenko, Andrey Pepelyshev, and Bader Saidan*
> 
> 2. **Exotic options in fractal activity time models with the Student distribution of log-returns**
>    *Authors: Nikolai N. Leonenko, Andrey Pepelyshev, and Bader Saidan*
> 
> 3. **Fractal Activity Time Risky Asset Models with Dependence and Heavy-Tailed Distributions** (PhD Thesis)
>    *Author: Bader Saidan*

## Overview

The **FATGBM** is a stochastic model that generalizes the classical Geometric Brownian Motion (GBM) by accommodating a number of properties that are not captured by standard Brownian motion but are frequently observed in empirical financial data. These include dependence in absolute and squared returns (but not the returns themselves) and a marginal distribution that is heavier-tailed and higher-peaked than Gaussian.

This repository provides:

1. **Stochastic Process Simulation:** Algorithms to simulate both the **Ornstein-Uhlenbeck (OU)** type process and the **superposition of Ornstein-Uhlenbeck (supOU)** process with inverse gamma marginals. This algorithm is used to simulate the fractal activity time using two different constructions: specifically, short-range dependence (OU-type) and long-range dependence (supOU).
2. **FATGBM Model:** Core functions to simulate the Fractal Activity Time Geometric Brownian Motion process.
3. **Option Pricing:** Numerical pricing methods implemented for the following contract types:
    * **Standard European Options:** Pricing for both **Call** and **Put** options (under FATGBM and GBM).
    * **Barrier Options:** Pricing for **Up-and-Out Call** and **Down-and-Out Put** options (under FATGBM and GBM).
4. **Reproduction Scripts:** Code to reproduce the figures and tables presented in the associated research papers and thesis.

## Repository Structure

### Core Models
* `RGamma_supOU_s.R`: Implementation of the supOU process with inverse gamma marginals used for activity time modeling.
* `FATGBM.R`: Implementation of the FATGBM model under both Long-Range Dependence (LRD) and Short-Range Dependence (SRD) frameworks. The script simulates the fractional activity time, $T_Y$, fitting it to an ex-Gaussian distribution for LRD and a Skew-Normal distribution for SRD. It also provides exact and Monte Carlo pricing functions for Standard European (Call/Put) and Barrier (Up-and-Out, Down-and-Out) options.
* `RGamma_OU.R`: Implementation of the OU-type process with inverse gamma marginals used for activity time modeling.
  
### Figures (Reproduction)
The following scripts generate the figures found in the first paper:
* `fig123.R`: Generates Figures 1, 2, and 3.
* `fig5.R` - `fig7.R`: Generates Figures 5 through 7.
* `fig8-9.R`: Generates Figures 8 and 9.
* `fig10.R` & `fig11.R`: Generates Figures 10 and 11.

### Tables (Pricing & Analysis)
The following scripts generate the tables found in the first paper:
* `table2.R`: Code to reproduce Table 2 (Parameter estimation/comparison).
* `table3pricing.R`: Code to reproduce Table 3 (Option pricing comparison).

## Prerequisites

To run these scripts, you must have **R** installed. 
You also require the following R packages, which can be installed via `install.packages("package_name")`:

* `gamlss`
* `gamlss.dist`
* `invgamma`
* `fitdistrplus`
* `ghyp`
* `moments`
* `zoo`

## Usage

### 1. Simulating the Model
To use the core functions, ensure you source the dependent scripts first. The FATGBM models depend on the respective activity time models:

```r
source("RGamma_supOU_s.R")
source("FATGBM.R")
source("RGamma_OU.R")

