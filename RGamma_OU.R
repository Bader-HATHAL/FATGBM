### =========================================================================
### Simulation of the OU process with inverse gamma (Reciprocal Gamma) marginals 
### =========================================================================

### =========================================================================
### Part 1: Helper Functions for OU Generation 
### =========================================================================

#' Calculate the g-function based on Bessel functions.
#' This function is a core component of the Lévy measure density for the GIG family.
#' @param x A numeric vector of positive integration values.
#' @param a The index parameter (nu) for the Bessel functions.
#' @return A numeric vector representing the evaluated g-function.
gFunction <- function(x, a) {
  sx <- sqrt(x)
  bJ <- besselJ(sx, a)
  bY <- besselY(sx, a)
  2 / (x * pi^2 * (bJ^2 + bY^2))
}

### =========================================================================
### Part 2: Random Variate Generation for DCRGamma
### =========================================================================

#' Sample from the Deconvolution of Reciprocal Gamma (DCRGamma) distribution.
#' This generates the innovations (Lévy increments) for the Reciprocal Gamma O-U process.
#' @param zeta The discounting parameter, typically exp(-lambda * dt / 2).
#' @param a The shape parameter (nu) of the Inverse Gamma marginal.
#' @param b The scale parameter (alpha) of the Inverse Gamma marginal.
#' @param n The number of internal exponential variables for the series representation.
#' @param N The number of integration points for the Riemann sum approximation.
#' @return A single numeric random variate from the DCRGamma distribution.
sampleFromDCRGamma <- function(zeta, a, b, n = 50, N = 600) {
  E <- rexp(n)
  S <- cumsum(E)
  
  # Calculate the intensity mapping for the jumps
  J0 <- 4 * (1 - zeta)^2 * b / (pi * S^2)
  W <- runif(n)
  
  Integ <- numeric(n)
  uu <- rev((1:N) / (N + 1))
  mloguu <- -log(uu)
  
  # Numerical integration to evaluate survival probabilities
  for(i in seq_len(n)) {
    z <- J0[i] / (4 * b)
    xx <- mloguu / z
    dd <- c(xx[1], diff(xx))
    
    integrand <- (exp(-J0[i] * xx / (4 * b)) - exp(-J0[i] * xx / (4 * b * zeta^2))) * gFunction(xx, a)
    Integ[i] <- sum(dd * integrand)
  }
  
  # Rejection step
  LHS <- sqrt(pi * J0) * Integ / (2 * (1 - zeta) * sqrt(b))
  J0[LHS < W] <- 0
  
  return(sum(J0))
}

### =========================================================================
### Part 3: Core Inverse Gamma OU Simulation
### =========================================================================

#' Simulate an Inverse Gamma Ornstein-Uhlenbeck process.
#' @param a The shape parameter (nu) of the marginal distribution.
#' @param b The scale parameter (alpha) of the marginal distribution.
#' @param lambda The mean reversion rate (decay parameter).
#' @param T The total time for the simulation.
#' @param dt The time increment (Delta).
#' @param n_burn The number of burn-in steps to reach stationarity.
#' @return A list containing the time vector (t) and the simulated process values (X).
rReciprocalGammaOU <- function(a, b, lambda, T, dt, n_burn = 10) {
  t <- seq(0, T, dt)
  n_steps <- length(t)
  X <- numeric(n_steps) # Pre-allocate vector for memory efficiency
  X0 <- 0
  
  # Pre-compute exponential decay terms
  zeta <- exp(-lambda * dt / 2)
  decay <- exp(-lambda * dt)
  
  # Warm-up / Burn-in period to ensure stationarity
  for(i in seq_len(n_burn)) {
    X0 <- X0 * decay + sampleFromDCRGamma(zeta, a, b)
  }
  
  X[1] <- X0
  
  # Path generation
  for(i in 2:n_steps) {
    X[i] <- X[i - 1] * decay + sampleFromDCRGamma(zeta, a, b)
  }
  
  return(list(t = t, X = X))
}

### =========================================================================
### Part 4: Process Initialization and Visualization
### =========================================================================

# Generate the process
ou_process <- rReciprocalGammaOU(a = 3, b = 2, lambda = 10, T = 11, dt = 0.01)

# Plot the trajectory
plot(ou_process$t, ou_process$X, 
     type = "l", 
     col = "darkblue", 
     xlab = "Time (t)", 
     ylab = "X(t)", 
     main = "Simulation of Reciprocal Gamma O-U Process")
### =========================================================================
### Part 5: Precomputation and Random Variate Generation
### =========================================================================
# (Note: gFunction is assumed to be loaded from Part 1)

#' Precompute the Inverse CDF for the Inverse Gamma Innovations
#' This generates a static lookup table to allow O(1) random sampling.
#' @param nu The shape parameter.
#' @param rho The scale parameter.
#' @param lambdaEps The intensity parameter for the Poisson jumps.
#' @param xn Number of integration points for the CDF grid.
#' @param un Number of quantiles for the final inverse CDF array.
#' @return A numeric vector representing the interpolated inverse CDF.
buildInvCdfLookup <- function(nu, rho, lambdaEps, xn = 400, un = 400) {
  # Define the upper bound using the base inverse gamma quantile
  xm <- qinvgamma(0.999, nu, rho)
  
  # Construct a non-linear grid to resolve the heavy left tail
  xv <- (1:xn) * (xm / xn)
  xv <- c(xv[1] * seq(0.09, 0.98, length.out = 50)^2, xv)
  xn_total <- length(xv)
  
  exp_gv <- function(xi, x) {
    exp(-1 / (4 * rho) * x * xi) * gFunction(xi, nu) / 2
  }
  
  cdfF <- numeric(xn_total)
  for(i in seq_len(xn_total)) {
    # Use tryCatch internally if needed to handle integration instability
    cdfF[i] <- 1 - (1 / lambdaEps) * integrate(exp_gv, 0, Inf, x = xv[i])$value
  }
  
  # Create the inverse mapping
  uv <- ((1:un) - 0.5) / un
  invCdf <- approx(cdfF, xv, uv)$y
  invCdf[is.na(invCdf)] <- 0
  
  return(invCdf)
}

#' Sample from the precomputed Inverse CDF lookup table
#' @param n The number of random variates to generate.
#' @param invCdfArray The precomputed numeric vector from buildInvCdfLookup.
#' @return A numeric vector of length n containing the random variates.
sampleFromLookup <- function(n, invCdfArray) {
  # O(1) array lookup using uniform random indexing
  max_index <- length(invCdfArray)
  indices <- ceiling(max_index * runif(n))
  return(invCdfArray[indices])
}


### =========================================================================
### Part 6: Core Simulation (Fast Shot-Noise Method)
### =========================================================================

#' Fast Simulation of the Inverse Gamma OU Process
#' Utilizes a shot-noise representation and vectorized path updates.
#' @param lambdaEps Jump intensity parameter.
#' @param A Mean reversion/decay parameter.
#' @param dt Time step size.
#' @param Y Total simulation time (formerly T).
#' @param T_min Warm-up / burn-in period to ensure stationarity.
#' @param invCdfArray The precomputed lookup table for jump sizes.
#' @return A list containing the time vector (t) and process values (X).
rReciprocalGammaOU_Fast <- function(lambdaEps, A, dt, Y, T_min, invCdfArray) {  
  t <- seq(0, Y, dt)
  n <- length(t) 
  X_sup <- numeric(n) # Pre-allocate vector
  
  # Simulate the first arrival time
  arrival <- T_min + rexp(1, rate = A * lambdaEps) 
  
  while(arrival < t[n]) {
    Z <- sampleFromLookup(1, invCdfArray) 
    ind <- floor(arrival / dt) + 2
    
    # Vectorized path updating
    if(arrival < 0) {
      X_sup <- X_sup + Z * exp(A * (arrival - t)) 
    } else if (ind <= n) {
      X_sup[ind:n] <- X_sup[ind:n] + Z * exp(A * (arrival - t[ind:n]))
    }
    
    # Simulate next arrival time
    arrival <- arrival + rexp(1, rate = A * lambdaEps)
  }
  
  return(list(t = t, X = X_sup))
}