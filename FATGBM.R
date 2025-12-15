# Load required libraries
# NOTE: 'gamlss' and 'gamlss.dist' are required for histDist and dexGAUS
library(gamlss)
library(gamlss.dist) 
# NOTE: the package that provides 'r_supOU'
source("RGamma_supOU_s.R") 


##################################################################
# Section 1: FATGBM Simulation
##################################################################

#' Simulate a FATGBM Path
#'
#' @param S0 Initial stock price.
#' @param r Risk-free interest rate.
#' @param sigma Volatility.
#' @param Y Total calendar time (in years).
#' @param dt Time step (e.g., 1/252).
#' @param InvFeVector Parameters for the subordinator's Levy measure.
#' @param piMeasure Parameters for the subordinator's Levy measure.
#'
#' @return A list containing the simulated path `St` and the `ou` process.
#' @export
simulate_FATGBM <- function(S0, r, sigma, Y, dt, InvFeVector, piMeasure) {
  ou <- r_supOU(InvFeVector, piMeasure, dt = dt, T = Y) 
  ou$Y[1] <- 0
  Tt <- cumsum(dt * ou$Y)
  dW <- rnorm(length(Tt), mean = 0, sd = sigma * sqrt(dt * ou$Y))
  W <- cumsum(dW)
  St <- S0 * exp(r * ou$t - (sigma^2 / 2) * Tt + W)
  return(list(St = St, ou = ou))
}

#' Get the Distribution of T_Y
#'
#' Simulates the fractal activity time T_Y and fits an ex-Gaussian distribution.
#'
#' @param Y Time to expiry (in years).
#' @param InvFeVector Parameters for r_supOU.
#' @param piMeasure Parameters for r_supOU.
#' @param m Number of simulations.
#' @param DaysInYear Trading days in a year.
#' @param seedn Random seed for reproducibility.
#'
#' @return A list containing simulated T_Y values, the histogram object, and
#'         the estimated ex-Gaussian parameters.
#' @export
getTY <- function(Y, InvFeVector, piMeasure, m = 500, 
                  DaysInYear = 252, seedn = 127) {
  set.seed(seedn)
  dt <- 1 / DaysInYear
  ty <- rep(0, m)
  for (j in (1:m)) {
    ou <- r_supOU(InvFeVector, piMeasure, dt = dt, T = Y)
    ou$Y[1] <- 0
    Tt <- cumsum(dt * ou$Y) # Fractal activity time
    ty[j] <- Tt[length(Tt)] # T_Y
  }
  # Fit ex-Gaussian distribution to T_Y and plot the fit
  hh <- histDist(ty, "exGAUS", nbins = 70)
  est <- c(round(hh$mu, 2), round(hh$sigma, 3), round(hh$nu, 2))
  s1 <- paste0('a=', InvFeVector$a, ' alpha=', piMeasure$prm)
  s2 <- paste('exGAUS(', est[1], ',', est[2], ',', est[3], ')', sep = "")
  print(paste(s1, s2))
  mtext(s1, side = 3, line = -1)
  xx <- seq(0, 2, by = 0.01)
  ee <- dexGAUS(xx, hh$mu, hh$sigma, hh$nu)
  lines(xx, ee, type = "l", col = "blue", lwd = 3)
  return(list(TY = ty, h = hh, prm = est))
}

#' Monte Carlo Simulation for Barrier Options using FATGBM Model
#'
#' Simulates barrier options (Up-and-Out Call and Put) prices using the FATGBM model
#' via Monte Carlo.
#'
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param B Barrier price (Up-and-Out).
#' @param r Risk-free interest rate.
#' @param sigma Volatility.
#' @param Y Time to expiry (in years).
#' @param dt Time step (e.g., 1/252 for daily steps).
#' @param InvFeVector Parameters for the subordinator's Lévy measure (for simulate_FATGBM).
#' @param piMeasure Parameters for the subordinator's Lévy measure (for simulate_FATGBM).
#' @param m Number of Monte Carlo simulations.
#' @param conf_level Confidence level for the interval (e.g., 0.95 for 95% CI).
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{CALL.m}: Estimated price of the Barrier Call option.
#'     \item \code{CALL.SE}: Standard error of the Barrier Call option price.
#'     \item \code{CALL.CI_lower}: Lower bound of the confidence interval for Call.
#'     \item \code{CALL.CI_upper}: Upper bound of the confidence interval for Call.
#'     \item \code{PUT.m}: Estimated price of the Barrier Put option.
#'     \item \code{PUT.SE}: Standard error of the Barrier Put option price.
#'     \item \code{PUT.CI_lower}: Lower bound of the confidence interval for Put.
#'     \item \code{PUT.CI_upper}: Upper bound of the confidence interval for Put.
#'   }
#' @export
simulate_BarrierOption <- function(S0, K, B, r, sigma, Y, dt, 
                                   InvFeVector, piMeasure, 
                                   m = 10000,           # Increased default simulations for better accuracy
                                   conf_level = 0.95,   # Default 95% confidence interval
                                   seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }
  # --- Initial Checks ---
  # For an Up-and-Out barrier, if S0 is already above the barrier, the option is worthless
  if (S0 >= B) {
    message("Warning: Initial stock price S0 (", S0, ") is already at or above barrier B (", B, ").")
    message("Returning 0 for barrier option prices as the option is immediately knocked out.")
    return(list(
      CALL.m = 0, CALL.SE = 0, CALL.CI_lower = 0, CALL.CI_upper = 0,
      PUT.m = 0, PUT.SE = 0, PUT.CI_lower = 0, PUT.CI_upper = 0
    ))
  }
  # --- Pre-calculations ---
  # Discount factor
  discount_factor <- exp(-r * Y)
  # Quantile for confidence interval (e.g., 1.96 for 95% CI)
  z_alpha <- qnorm(1 - (1 - conf_level) / 2)
  # Initialize payoff vectors
  call_payoffs <- numeric(m)
  put_payoffs <- numeric(m)
  # --- Monte Carlo Loop ---
  message("Starting Monte Carlo simulation for Barrier Options (Up-and-Out)...")
  pb <- txtProgressBar(min = 0, max = m, initial = 0, style = 3) # Style 3 for dynamic bar
  for (i in 1:m) {
    fatgbm <- simulate_FATGBM(S0, r, sigma, Y, dt, InvFeVector, piMeasure)
    spath <- fatgbm$St
    # Barrier condition: Check if max price *with a slight adjustment* exceeds barrier B
    # NOTE: The factor 0.583 * volsdt is a specific approximation for continuous barrier monitoring.
    # For discrete barrier monitoring, you would typically use `max(spath) >= B` or `any(spath >= B)`.
    volsdt <- sigma * sqrt(dt)
    barrier_breached <- (max(spath) * exp(0.583 * volsdt) >= B)
    if (!barrier_breached) {
      final_price <- spath[length(spath)]
      call_payoffs[i] <- pmax(final_price - K, 0)
      put_payoffs[i] <- pmax(K - final_price, 0)
    }
    # If barrier_breached is TRUE, payoffs remain 0 (from initialization)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # --- Calculate Statistics ---
  # Call Option
  mean_call_payoff <- mean(call_payoffs)
  sd_call_payoff <- sd(call_payoffs)
  call_price <- discount_factor * mean_call_payoff
  call_se <- discount_factor * sd_call_payoff / sqrt(m) # CORRECTED SE calculation
  call_ci_lower <- call_price - z_alpha * call_se
  call_ci_upper <- call_price + z_alpha * call_se
  # Put Option
  mean_put_payoff <- mean(put_payoffs)
  sd_put_payoff <- sd(put_payoffs)
  put_price <- discount_factor * mean_put_payoff
  put_se <- discount_factor * sd_put_payoff / sqrt(m) # CORRECTED SE calculation
  put_ci_lower <- put_price - z_alpha * put_se
  put_ci_upper <- put_price + z_alpha * put_se
  # --- Output Messages ---
  message("\n--- Monte Carlo Simulation Results (FATGBM Barrier Options) ---")
  message("Parameters: K=", K, ", B=", B, ", S0=", S0, ", Y=", Y, ", r=", r, ", sigma=", sigma)
  message("            InvFeVector$a=", InvFeVector$a, ", piMeasure$prm=", piMeasure$prm, ", m=", m)
  message(sprintf("Barrier Call Price = %.2f (SE = %.2f, %g%% CI: [%.2f, %.2f])", 
                  call_price, call_se, conf_level * 100, call_ci_lower, call_ci_upper))
  message(sprintf("Barrier Put Price  = %.2f (SE = %.2f, %g%% CI: [%.2f, %.2f])", 
                  put_price, put_se, conf_level * 100, put_ci_lower, put_ci_upper))
  # --- Return Results ---
  return(list(
    CALL.m = call_price, 
    CALL.SE = call_se, 
    CALL.CI_lower = call_ci_lower, 
    CALL.CI_upper = call_ci_upper,
    PUT.m = put_price, 
    PUT.SE = put_se, 
    PUT.CI_lower = put_ci_lower, 
    PUT.CI_upper = put_ci_upper
  ))
}

#####################################################################
# Section 2: Call Option Pricing
#####################################################################

# Internal helper function for FATUpOutCall
# This is the integrand for the (K <= B) case
FATUpOutCall_subf <- function(t, Y, S0, K, B, sigma, r) {
  b <- 1 / sigma * (log(B / S0))
  k <- 1 / sigma * (log(K / S0))
  rys <- r * Y / sigma
  ybr <- 2 * Y * b * r / (t * sigma)
  
  d1 <- function(t, sigma, rys, k) (rys + sigma * t / 2 - k) / (sqrt(t))
  d2 <- function(t, sigma, rys, b) (rys + sigma * t / 2 - b) / (sqrt(t))
  d3 <- function(t, sigma, rys, k, b) (rys + sigma * t / 2 + 2 * b - k) / (sqrt(t))
  d4 <- function(t, sigma, rys, b) (rys + sigma * t / 2 + b) / (sqrt(t))
  d5 <- function(t, sigma, rys, k) (rys - sigma * t / 2 - k) / (sqrt(t))
  d6 <- function(t, sigma, rys, b) (rys - sigma * t / 2 - b) / (sqrt(t))
  d7 <- function(t, sigma, rys, k, b) (rys - sigma * t / 2 + 2 * b - k) / (sqrt(t))
  d8 <- function(t, sigma, rys, b) (rys - sigma * t / 2 + b) / (sqrt(t))
  
  term1 <- S0 * exp(r * Y) *
    (pnorm(d1(t, sigma, rys, k)) - pnorm(d2(t, sigma, rys, b)))
  term2 <- S0 * exp(r * Y) * exp(b * sigma + ybr) *
    (pnorm(d3(t, sigma, rys, k, b)) - pnorm(d4(t, sigma, rys, b)))
  term3 <- K *
    (pnorm(d5(t, sigma, rys, k)) - pnorm(d6(t, sigma, rys, b)))
  term4 <- K * exp(-b * sigma + ybr) *
    (pnorm(d7(t, sigma, rys, k, b)) - pnorm(d8(t, sigma, rys, b)))
  term2[is.na(term2)] <- 0
  term4[is.na(term4)] <- 0
  price <- term1 - term2 - term3 + term4
  return(price)
}

#' Price a FATGBM Barrier Call Option (Up-and-Out)
#'
#' @param Y Time to expiry (in years).
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param B Barrier price.
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param dTYprm Parameters of the ex-Gaussian density for T_Y.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
FATUpOutCall <- function(Y, S0, K, B, sigma, r, dTYprm, echo = FALSE) {
  # Standard check: if barrier is already breached, option is worthless.
  if (B < S0) return(0) 
  if (K > B) {
    #warning("FATUpOutCall formula is only implemented for K <= B. Returning 0.")
    return(0)
  }
  integrandc <- function(t, dTYprm) {
    density <- dexGAUS(t, dTYprm[1], dTYprm[2], dTYprm[3])
    cbbc_price <- FATUpOutCall_subf(t, Y, S0, K, B, sigma, r)
    density * cbbc_price
  }
  res <- integrate(integrandc, lower = 0, upper = Inf, dTYprm = dTYprm)
  price <- res$value
  if (echo) message("FATGBM Barrier Call Option Price: ", round(price, 4))
  return(price)
}

#' Price a GBM Barrier Call Option (Up-and-Out)
#'
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param B Barrier price.
#' @param Y Time to expiry (in years).
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
GBMUpOutCall <- function(S0, K, B, Y, sigma, r, echo = FALSE) {
  if (B < S0) return(0)
  dp <- function(Y, r, sigma, s) {
    (log(s) + (r + 0.5 * sigma^2) * Y) / (sigma * sqrt(Y))
  }
  dm <- function(Y, r, v, s) {
    (log(s) + (r - 0.5 * sigma^2) * Y) / (sigma * sqrt(Y))
  }
  term1 <- S0 * (pnorm(dp(Y, r, sigma, S0 / K)) - pnorm(dp(Y, r, sigma, S0 / B)))
  term2 <- S0 * (B / S0)^(1 + 2 * r / sigma^2) *
    (pnorm(dp(Y, r, sigma, B^2 / (K * S0))) - pnorm(dp(Y, r, sigma, B / S0)))
  M34 <- K * exp(-r * Y)
  term3 <- M34 * (pnorm(dm(Y, r, sigma, S0 / K)) - pnorm(dm(Y, r, sigma, S0 / B)))
  term4 <- M34 * (S0 / B)^(1 - 2 * r / sigma^2) *
    (pnorm(dm(Y, r, sigma, B^2 / (K * S0))) - pnorm(dm(Y, r, sigma, B / S0)))
  price <- term1 - term2 - (term3 - term4)
  if (echo) message("   GBM Barrier Call Option Price: ", round(price, 4))
  return(price)
}

#' Price a GBM Standard European Call Option
#'
#' @param Y Time to expiry (in years).
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
GBMEurCall <- function(Y, S0, K, sigma, r, echo = FALSE) {
  cd1 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y + (sigma^2 / 2) * t) / (sigma * sqrt(t))
  cd2 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y - (sigma^2 / 2) * t) / (sigma * sqrt(t))
  BS <- function(t, Y, S0, K, sigma, r) {
    S0 * pnorm(cd1(t, Y, S0, K, sigma, r)) -
      K * exp(-r * Y) * pnorm(cd2(t, Y, S0, K, sigma, r))
  }
  price <- BS(Y, Y, S0, K, sigma, r)
  if (echo) message("   GBM Call Option Price: ", round(price, 4))
  return(price)
}

#' Price a FATGBM Standard European Call Option
#'
#' @param Y Time to expiry (in years).
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param dTYprm Parameters of the ex-Gaussian density for T_Y.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
FATEurCall <- function(Y, S0, K, sigma, r, dTYprm, echo = FALSE) {
  cd1 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y + (sigma^2 / 2) * t) / (sigma * sqrt(t))
  cd2 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y - (sigma^2 / 2) * t) / (sigma * sqrt(t))
  BS <- function(t, Y, S0, K, sigma, r) {
    S0 * pnorm(cd1(t, Y, S0, K, sigma, r)) -
      K * exp(-r * Y) * pnorm(cd2(t, Y, S0, K, sigma, r))
  }
  integrandc <- function(t, dTYprm) {
    density <- dexGAUS(t, dTYprm[1], dTYprm[2], dTYprm[3])
    density * BS(t, Y, S0, K, sigma, r)
  }
  res <- integrate(integrandc, lower = 0, upper = Inf, dTYprm = dTYprm)
  price <- res$value
  if (echo) message("FATGBM Call Option Price: ", round(price, 4))
  return(price)
}

#####################################################################
# Section 3: Put Option Pricing
#####################################################################

# Internal helper: FATGBM Barrier Put (Up-and-Out) (K <= B)
FATUpOutPut_subf_BaK <- function(t, Y, S0, K, B, sigma, r) {
  b <- 1 / sigma * (log(B / S0))
  k <- 1 / sigma * (log(K / S0))
  rys <- r * Y / sigma
  ybr <- 2 * Y * b * r / (t * sigma)
  d1 <- function(t, sigma, rys, k) (rys - sigma * t / 2 - k) / (sqrt(t))
  d2 <- function(t, sigma, rys, k, b) (rys - sigma * t / 2 + 2 * b - k) / (sqrt(t))
  d3 <- function(t, sigma, rys, k) (rys + sigma * t / 2 - k) / (sqrt(t))
  d4 <- function(t, sigma, rys, k, b) (rys + sigma * t / 2 + 2 * b - k) / (sqrt(t))
  term1 <- K * pnorm(-d1(t, sigma, rys, k))
  term2 <- K * exp(-b * sigma + ybr) * pnorm(-d2(t, sigma, rys, k, b))
  term3 <- S0 * exp(r * Y) * pnorm(-d3(t, sigma, rys, k))
  term4 <- S0 * exp(r * Y) * exp(b * sigma + ybr) * pnorm(-d4(t, sigma, rys, k, b))
  term2[is.na(term2)] <- 0
  term4[is.na(term4)] <- 0
  price <- term1 - term2 - term3 + term4
  return(price)
}

# Internal helper: FATGBM Barrier Put (Up-and-Out) (K >= B)
FATUpOutPut_subf_KaB <- function(t, Y, S0, K, B, sigma, r) {
  b <- 1 / sigma * (log(B / S0))
  rys <- r * Y / sigma
  ybr <- 2 * Y * b * r / (t * sigma)
  d1 <- function(t, sigma, rys, b) (rys - sigma * t / 2 - b) / (sqrt(t))
  d2 <- function(t, sigma, rys, b) (rys - sigma * t / 2 + b) / (sqrt(t))
  d3 <- function(t, sigma, rys, b) (rys + sigma * t / 2 - b) / (sqrt(t))
  d4 <- function(t, sigma, rys, b) (rys + sigma * t / 2 + b) / (sqrt(t))
  term1 <- K * pnorm(-d1(t, sigma, rys, b))
  term2 <- K * exp(-b * sigma + ybr) * pnorm(-d2(t, sigma, rys, b))
  term3 <- S0 * exp(r * Y) * pnorm(-d3(t, sigma, rys, b))
  term4 <- S0 * exp(r * Y) * exp(b * sigma + ybr) * pnorm(-d4(t, sigma, rys, b))
  term2[is.na(term2)] <- 0
  term4[is.na(term4)] <- 0
  price <- term1 - term2 - term3 + term4
  return(price)
}

#' Price a FATGBM Barrier Put Option (Up-and-Out)
#'
#' @param Y Time to expiry (in years).
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param B Barrier price.
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param dTYprm Parameters of the ex-Gaussian density for T_Y.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
FATUpOutPut <- function(Y, S0, K, B, sigma, r, dTYprm, echo = FALSE) {
  # Standard check: if barrier is already breached, option is worthless.
  if (B < S0) return(0)
  integrandc <- function(t, dTYprm) {
    density <- dexGAUS(t, dTYprm[1], dTYprm[2], dTYprm[3])
    if (K <= B) {
      cbbc_price <- FATUpOutPut_subf_BaK(t, Y, S0, K, B, sigma, r)
    } else {
      cbbc_price <- FATUpOutPut_subf_KaB(t, Y, S0, K, B, sigma, r)
    }
    density * cbbc_price
  }
  res <- integrate(integrandc, lower = 0, upper = Inf, dTYprm = dTYprm)
  price <- res$value
  if (echo) message("FATGBM Barrier Put Option Price: ", round(price, 4))
  return(price)
}

#' Price a FATGBM Standard European Put Option
#'
#' @param Y Time to expiry (in years).
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param dTYprm Parameters of the ex-Gaussian density for T_Y.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
FATEurPut <- function(Y, S0, K, sigma, r, dTYprm, echo = FALSE) {
  cd1 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y + (sigma^2 / 2) * t) / (sigma * sqrt(t))
  cd2 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y - (sigma^2 / 2) * t) / (sigma * sqrt(t))
  BS <- function(t, Y, S0, K, sigma, r) {
    K * exp(-r * Y) * pnorm(-cd2(t, Y, S0, K, sigma, r)) -
      S0 * pnorm(-cd1(t, Y, S0, K, sigma, r))
  }
  integrandc <- function(t, dTYprm) {
    density <- dexGAUS(t, dTYprm[1], dTYprm[2], dTYprm[3])
    density * BS(t, Y, S0, K, sigma, r)
  }
  res <- integrate(integrandc, lower = 0, upper = Inf, dTYprm = dTYprm)
  price <- res$value
  if (echo) message("FATGBM Put Option Price: ", round(price, 4))
  return(price)
}

#' Price a GBM Barrier Put Option (Up-and-Out)
#'
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param B Barrier price.
#' @param Y Time to expiry (in years).
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
GBMUpOutPut <- function(S0, K, B, Y, sigma, r, echo = FALSE) {
  if (B < S0) return(0)
  dp <- function(Y, r, sigma, s) {
    (log(s) + (r + 0.5 * sigma^2) * Y) / (sigma * sqrt(Y))
  }
  dm <- function(Y, r, sigma, s) {
    (log(s) + (r - 0.5 * sigma^2) * Y) / (sigma * sqrt(Y))
  }
  M34 <- K * exp(-r * Y)
  if (K <= B) {
    term1 <- M34 * pnorm(-dm(Y, r, sigma, S0 / K))
    term2 <- M34 * (S0 / B)^(1 - 2 * r / sigma^2) *
      pnorm(-dm(Y, r, sigma, B^2 / (K * S0)))
    term3 <- S0 * pnorm(-dp(Y, r, sigma, S0 / K))
    term4 <- S0 * (B / S0)^(1 + 2 * r / sigma^2) *
      pnorm(-dp(Y, r, sigma, B^2 / (K * S0)))
  } else {
    term1 <- M34 * pnorm(-dm(Y, r, sigma, S0 / B))
    term2 <- M34 * (S0 / B)^(1 - 2 * r / sigma^2) *
      pnorm(-dm(Y, r, sigma, B / S0))
    term3 <- S0 * pnorm(-dp(Y, r, sigma, S0 / B))
    term4 <- S0 * (B / S0)^(1 + 2 * r / sigma^2) *
      pnorm(-dp(Y, r, sigma, B / S0))
  }
  price <- term1 - term2 - (term3 - term4)
  if (echo) message("   GBM Barrier Put Option Price: ", round(price, 4))
  return(price)
}

#' Price a GBM Standard European Put Option
#'
#' @param Y Time to expiry (in years).
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param sigma Volatility.
#' @param r Risk-free interest rate.
#' @param echo Logical, whether to print the price.
#'
#' @return Option price.
#' @export
GBMEurPut <- function(Y, S0, K, sigma, r, echo = FALSE) {
  cd1 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y + (sigma^2 / 2) * t) / (sigma * sqrt(t))
  cd2 <- function(t, Y, S0, K, sigma, r)
    (log(S0 / K) + r * Y - (sigma^2 / 2) * t) / (sigma * sqrt(t))
  BS <- function(t, Y, S0, K, sigma, r) {
    K * exp(-r * Y) * pnorm(-cd2(t, Y, S0, K, sigma, r)) -
      S0 * pnorm(-cd1(t, Y, S0, K, sigma, r))
  }
  price <- BS(Y, Y, S0, K, sigma, r)
  if (echo) message("   GBM Put Option Price: ", round(price, 4))
  return(price)
}

