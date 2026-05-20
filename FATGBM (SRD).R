# Load required libraries
# NOTE: 'gamlss' and 'gamlss.dist' are required for histDist and SN
library(gamlss)
library(gamlss.dist) 
# NOTE: the package that provides 'Rg.OU'
source("RGamma_OU.R") 

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
#' @param a The shape parameter  of the Inverse Gamma marginal.
#' @param b The scale parameter of the Inverse Gamma marginal.
#' @param lambdaEps Jump intensity parameter.
#'
#' @return A list containing the simulated path `St` and the `ou` process.
#' @export
simulate_FATGBM <- function(S0, r, sigma, Y, n, a, b, lambdaEps, Delt){
  dt = Y / n  
  t = seq(0, Y, by = dt) 
  ou=Rg.OU(lambdaEps=lambdaEps,A=A,dt=dt,Y=Y,T_min=-1)
  Tt = cumsum(dt * ou$X)
  dW = rnorm(n, mean = 0, sd = sigma * sqrt(dt * ou$X))
  W = c(0,cumsum(dW))
  St = S0 * exp(r * ou*t - (sigma^2 / 2) * Tt + W)
  return(list(St = St, ou = u))
}
#' Get the Distribution of T_Y
#'
#' Simulates the fractal activity time T_Y and fits a skewed normal distribution.
#'
#' @param Y Time to expiry (in years).
#' @param InvFeVector Parameters for r_supOU.
#' @param piMeasure Parameters for r_supOU.
#' @param m Number of simulations.
#' @param DaysInYear Trading days in a year.
#' @param seedn Random seed for reproducibility.
#'
#' @return A list containing simulated T_Y values, the histogram object, and
#'         the estimated skewed normal parameters.
#' @export
#' 
getTY <- function(Y, lambdaEps, lambda, m = 500, DaysInYear = 252, seedn = 127){
set.seed(seedn)
dt <- 1 / DaysInYear
ty <- rep(0, m)
for (j in (1:m)) {
  ou=Rg.OU(lambdaEps=lambdaEps, lambda=lambda, dt=dt, Y=Y, T_min=-1)
  ou$Y[1] <- 0
  Tt <- cumsum(dt * ou$X) # Fractal activity time
  ty[j] <- Tt[length(Tt)] # T_Y
}
# Fit skewed normal distribution to T_Y and plot the fit
hh <- histDist(ty, "SN2", nbins = 70) 
est <- c(round(hh$mu, 2), round(hh$sigma, 3), round(hh$nu, 2))
s1 <- paste0('nu=',nu, 'lambda',lambda)
s2 <- paste('SN(', est[1], ',', est[2], ',', est[3], ')', sep = "")
print(paste(s1, s2))
mtext(s1, side = 3, line = -1)
xx <- seq(0, 2, by = 0.01)
ee <- dSN2(xx, hh$mu, hh$sigma, hh$nu)
lines(xx, ee, type = "l", col = "blue", lwd = 3)
return(list(TY = ty, h = hh, prm = est))
}

#' Monte Carlo Simulation for Barrier Options using FATGBM Model
#'
#' Simulates barrier options (Up-and-Out Call and Put) (Down-and-Out Call and Put) prices using the FATGBM model
#' via Monte Carlo.
#'
#' @param S0 Initial stock price.
#' @param K Strike price.
#' @param B Barrier price (Up-and-Out and Down-and-Out Call and Put).
#' @param r Risk-free interest rate.
#' @param sigma Volatility.
#' @param Y Time to expiry (in years).
#' @param dt Time step (e.g., 1/252 for daily steps).
#' @param lambdaEps Jump intensity parameter.
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



