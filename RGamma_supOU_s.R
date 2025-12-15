### =========================================================================
### Simulation of the supOU process with inverse gamma marginals  
### =========================================================================

### =========================================================================
### Part 1: Helper Functions for supOU Generation
### =========================================================================

#' Calculate the g-function based on Bessel functions.
#' This is a core component of the Lévy measure density.
#' @param t A numeric value.
#' @param a The shape parameter alpha.
#' @return The value of the g-function.
gFunction <- function(t, a) {
  st <- sqrt(t)
  bJ <- besselJ(st, a)
  bY <- besselY(st, a)
  2 / (t * pi^2 * (bJ^2 + bY^2))
}

#' Create an integrand for the m-function.
#' This is an exponential transform of the g-function.
#' @param t The integration variable.
#' @param x A numeric value.
#' @param a The shape parameter alpha.
#' @param b The scale parameter beta.
#' @return The value of the integrand.
expGv <- function(t, x, a, b) {
  exp(-(x * t) / (4 * b)) * gFunction(t, a) / (2 * x)
}

#' Calculate the m-function via numerical integration.
#' @return The numeric value of the integral.
mFunction <- function(x, a, b) {
  # Use tryCatch to handle potential integration errors gracefully
  tryCatch({
    integrate(expGv, 0, Inf, x = x, a = a, b = b)$value
  }, error = function(e) {
    warning("Integration failed for x = ", x, ". Returning NA.")
    return(NA)
  })
}

#' Calculate the distribution function Fe.
#' @return A numeric vector of probabilities.
feFunction <- function(x, a, b, epsilon, thetaEps) {
  # Vectorized calculation of the Fe distribution function
  sapply(x, function(val) {
    if (val > epsilon) {
      1 - val / thetaEps * mFunction(val, a, b)
    } else {
      0
    }
  })
}

#' Construct the inverse CDF of the Fe distribution.
#' This is used for efficient random sampling.
#' @return A list containing the quantile function (Q) and thetaEps.
getInvCdf <- function(a, b, epsilon, n_points = 600, n_quantiles = 500) {
  thetaEps <- epsilon * mFunction(epsilon, a, b)

  # Find a reasonable upper bound for the distribution
  xMax <- 1
  while (feFunction(xMax, a, b, epsilon, thetaEps) < 0.99999) {
    xMax <- xMax * 1.5
  }
  # Create a grid and compute the CDF values
  # Use cube-root spacing to resolve heavy tails
  xVector <- seq(epsilon^(1/3), xMax^(1/3), length.out=n_points)^3
  # Vectorized the call to feFunction
  feVector <- feFunction(xVector, a, b, epsilon, thetaEps)
  # Create the inverse mapping using linear approximation
  uVector <- ((1:n_quantiles) - 0.5) / n_quantiles
  invFeVector <- approx(feVector, xVector, uVector, rule = 2)$y
  return(list(Q = invFeVector, thetaEps = thetaEps))
}


### =========================================================================
### Part 2: Random Variate Generation
### =========================================================================

#' Sample a random variate from the pre-computed inverse Fe CDF.
#' @param invFeVector A numeric vector representing the quantile function.
#' @return A single random number.
sampleFromFe <- function(invFeVector) {
  # Sample with replacement from the quantile vector
  invFeVector[sample.int(length(invFeVector), 1)]
}


#' Set the piMeasure list to ensure the correct 'eta'.
setPiMeasure <- function(type, prm, T_WarmUp=1, B=1, beta=1) {
  switch(type,
         "delta" = return(list(type=type, prm=prm, 
                               T_WarmUp=T_WarmUp,eta=1/prm)),
         "gamma" = return(list(type=type, prm=prm, 
                               T_WarmUp=T_WarmUp,eta=1)),
         "gamma2" = return(list(type=type,prm=prm, beta=beta,
                                T_WarmUp=T_WarmUp,eta=beta/prm)),
         "beta" = return(list(type=type, prm=prm, B=B, 
                              T_WarmUp=T_WarmUp,eta=(1+prm)/prm/B)),
         stop("Unknown piMeasure type specified.")
  )
}
  
#' Sample a random variate from the mixing measure pi.
#' Supports delta, gamma, and beta distributions.
#' @param piMeasure A list describing the mixing measure.
#' @return A single random number.
sampleFromPi <- function(piMeasure) {
  switch(piMeasure$type,
    "delta" = return(piMeasure$prm),
    "gamma"  = return(max(1e-4,rgamma(1, shape = piMeasure$prm + 1, 
                                      rate = piMeasure$prm))),
    "gamma2" = return(max(1e-4,rgamma(1, shape = piMeasure$prm + 1, 
                                      rate = piMeasure$beta))),
    "beta" = {
      d <- 1 / (1 + piMeasure$prm)
      u <- runif(1, 1e-4, 1)
      return((u^d) * piMeasure$B)
    },
    stop("Unknown piMeasure type specified.")
  )
}


### =========================================================================
### Part 3: Core supOU Simulation and Analysis
### =========================================================================

#' Simulate a superposition of Ornstein-Uhlenbeck (supOU) processes.
#' @param invFe A list from getInverseCdf.
#' @param piMeasure A list describing the mixing measure.
#' @param dt Time step.
#' @param T Total simulation time.
#' @return A list containing the time vector (t) and process values (Y).
r_supOU <- function(invFe, piMeasure, dt, T) {
  intensityAdj <- invFe$thetaEps / piMeasure$eta
  t <- seq(0, T, dt)
  n <- length(t)
  Y <- numeric(n) # Pre-allocate vector for efficiency
  # Start with a warm-up period to ensure stationarity
  jumpTime <- -piMeasure$T_WarmUp + rexp(1, rate = intensityAdj)
  while (jumpTime < t[n]) {
    R_k <- sampleFromPi(piMeasure)   # Damping factor
    H_k <- sampleFromFe(invFe$Q)     # Jump height
    if (jumpTime < 0) {
      # Contribution from jumps during the warm-up period
      Y <- Y + H_k * exp(R_k * (jumpTime - t))
    } else {
      # Find the index where the jump occurs and update the rest of the path
      startIndex <- floor(jumpTime / dt) + 2
      if (startIndex <= n) {
        Y[startIndex:n] <- Y[startIndex:n] + 
          H_k * exp(R_k * (jumpTime - t[startIndex:n]))
      }
    }
    # Generate the next jump time
    jumpTime <- jumpTime + rexp(1, rate = intensityAdj)
  }
  return(list(t = t, Y = Y))
}

#' Estimate the decay rate 'alpha' from ACF via log-log regression.
#' @param acfObject An object returned by the acf() function.
#' @param plot Logical, whether to plot the log-log regression.
#' @return The estimated alpha value.
estimateAlpha <- function(acfObject, plot = TRUE) {
  # Exclude the first two lags (lag 0 and 1)
  acfValues <- acfObject$acf[-c(1, 2)]
  lags <- acfObject$lag[-c(1, 2)]
  # Truncate where ACF becomes too noisy (e.g., drops below 0.05)
  k <- which(acfValues < 0.02)
  if (length(k) > 0) {
    acfValues <- acfValues[1:k[1]]
    lags <- lags[1:k[1]]
  }
  # Log-transform for linear regression
  logLags <- log(lags)
  logAcf <- log(abs(acfValues))
  # Fit linear model: log(ACF) ~ log(Lag)
  model <- lm(logAcf ~ logLags)
  if (plot) { # Optional plot
    plot(logLags, logAcf,
      ylim = c(-5, -0.2),
      xlab = "log(Lag)", ylab = "log(ACF)",
      main = "Log-Log Plot of ACF for Alpha Estimation"
    )
    abline(model, col = "blue", lwd = 2)
    grid()
  }
  # The slope of the log-log plot is -alpha
  hatAlpha <- -coef(model)[2]
  return(hatAlpha)
}


### =========================================================================
### Part 4: "Smart" Simulation Wrapper 
### =========================================================================

#' Validate the piMeasure list to ensure it has all required components.
validatePiMeasure <- function(piMeasure) {
  required_names <- c("type", "prm", "eta", "T_WarmUp")
  if (!all(required_names %in% names(piMeasure))) {
    stop("piMeasure must be a list with elements: type, prm, eta, T_WarmUp")
  }
  if (piMeasure$type == "beta" && is.null(piMeasure$B)) {
    stop("piMeasure of type 'beta' must also include element 'B'")
  }
  if (piMeasure$type == "gamma2" && is.null(piMeasure$beta)) {
    stop("piMeasure of type 'gamma2' must also include element 'beta'")
  }
}

#' "Smart" supOU simulation that finds the best realization matching a target alpha.
#' This function over-simulates and selects the sub-series whose estimated alpha
#' is closest to the theoretical parameter.
#' @param invFe A list from getInvCdf.
#' @param piMeasure A list describing the mixing measure with target alpha (`prm`).
#' @param dt Time step.
#' @param T Total simulation time for the final output.
#' @return An object with the supOU process.
r_supOU_smart <- function(invFe, piMeasure, dt, T) {
  validatePiMeasure(piMeasure)
  # --- Over-simulation ---
  # Simulate a process 3x longer to find a good matching window
  longSim <- r_supOU(invFe, piMeasure, dt, T * 3)
  # --- Find the best matching window ---
  n <- round(T / dt) + 1 # Desired length of final series
  num_windows <- 20
  best_j <- 0
  best_diff <- Inf
  cat("Searching for best matching window...\n")
  for (j in 0:num_windows) {
    # Define the start index for the sliding window
    start_index <- 1 + floor(j / num_windows * (length(longSim$Y) - n))
    window_indices <- start_index:(start_index + n - 1)
    # Estimate alpha for the current window
    current_acf <- acf(longSim$Y[window_indices], lag.max = 250, plot = FALSE)
    est_alpha <- estimateAlpha(current_acf, plot = FALSE)
    current_diff <- abs(est_alpha - piMeasure$prm)
    cat(sprintf("Attempt %2d: Target Alpha=%.3f, Estimated Alpha=%.3f\n", j, piMeasure$prm, est_alpha))
    if (current_diff < best_diff) {
      best_diff <- current_diff
      best_j <- j
    }
  }
  cat("Best window found at iteration:", best_j, "with difference:", round(best_diff, 4), "\n")
  # Extract the best window
  start_index <- 1 + floor(best_j / num_windows * (length(longSim$Y) - n))
  window_indices <- start_index:(start_index + n - 1)
  result <- list(
    t = seq(0, T, dt),
    Y = longSim$Y[window_indices],
    params = list(a = a, b = b, piMeasure = piMeasure)
  )
  return(result)
}

