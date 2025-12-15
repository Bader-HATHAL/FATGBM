# Figure 8 and 9 in the paper
# Load required libraries 
library(invgamma)
library(fitdistrplus)
library(ghyp)
library(moments)
library(zoo)
source("RGamma_supOU_s.R") 
source("FATGBM.R")
  
ToPDF <- FALSE
ToPDF <- TRUE

# setup of parameters 
DaysInYear <- 252
S0 <- 100       
Y <- 1          
r <- 0.04       
sigma <- 0.3    
a <- 1.5
b <- a - 1
epsilon <- 0.001
dt <- 1 / DaysInYear
InvFeVector <- getInvCdf(a, b, epsilon)

# Set the measure for supOU simulation
beta <- 1/40 
alpha <- 0.6
piMeasure1 <- setPiMeasure(type='gamma2', prm=alpha, 
                           beta=beta, T_WarmUp=2)

if(ToPDF) pdf('fig8.pdf',width=6,height=2.5)
# Plot setup

# Simulate Multiple Realizations
m_paths <- 10  # number of realizations
seedn <- 138
set.seed(seedn)
paths_matrix <- matrix(NA, nrow = Y/dt+1, ncol = m_paths)
for (i in 1:m_paths) {
  result <- simulate_FATGBM(S0, r, sigma, Y, dt, InvFeVector, piMeasure1)
  paths_matrix[, i] <- result$St
}
t <- result$ou$t

par(bty="n",mar=c(2.1, 2.1, 1.1, 1.1))
matplot(t, paths_matrix, type = "l", lty = 1, col = rainbow(m_paths),
        xlab = "Time ", ylab = expression(S[t]), lwd=2)
box(bty = "n")
grid()

Y <- 20
seedn <- 135 
set.seed(seedn)
paths_matrix <- matrix(NA, nrow = Y/dt+1, ncol = m_paths)
for (i in 1:m_paths) {
  result <- simulate_FATGBM(S0, r, sigma, Y, dt, InvFeVector, piMeasure1)
  paths_matrix[, i] <- result$St
}
t <- result$ou$t

par(bty="n",mar=c(2.1, 2.1, 1.1, 1.1))
matplot(t, paths_matrix, type = "l", lty = 1, col = rainbow(m_paths),
        xlab = "Time ", ylab = expression(S[t]), lwd=2)
box(bty = "n")
grid()

LSSest <- matrix(0, nrow=m_paths, ncol=3)
GHest <- matrix(0, nrow=m_paths, ncol=5)
par(mar=c(2,2,1,1))
for (i in 1:m_paths) {
  log_returns <- diff(log(paths_matrix[, i]))
  ed <- density(log_returns, adjust = 1.25)
  if(i==1)
    plot(ed, log='y', xlim=c(-0.05,0.05), lwd=2,
         ylim=c(0.1,1.4*max(ed$y)), main="",
         col=rainbow(m_paths)[1])  
  else
    lines(ed, col = rainbow(m_paths)[i], lwd=2)
  print(paste("iteration",i))
  fitted_t <- fit.tuv(data = log_returns,symmetric = TRUE,silent=TRUE)
  est <- coef(fitted_t)
  LSSest[i,] <- c(est$nu,est$sigma,est$mu)
  fitted_ghyp <- fit.ghypuv(data = log_returns, mu = mean(log_returns), 
                            sigma = sd(log_returns), silent=TRUE)
  est <- coef(fitted_ghyp)
  GHest[i,] <- c(est$lambda,est$alpha.bar,est$gamma,est$sigma,est$mu)
}
log_returns <- diff(log(paths_matrix[, 5]))
# Fit normal distribution
fit_norm <- fitdistr(log_returns, "normal")
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
norm_density <- dnorm(x_vals, mean = fit_norm$estimate["mean"], 
                      sd = fit_norm$estimate["sd"])
lines(x_vals, norm_density, col = "black", lwd = 3)
mtext(paste0("Log-densities of Log-Returns"), line=0, side=3)


par(mar=c(2,2,1,1))
for (i in 1:m_paths) {
  log_returns <- diff(log(paths_matrix[, i]))
  acf_sq <- acf(log_returns^2, plot = FALSE, lag.max = 250)
  if(i==1)
    plot(acf_sq$lag, acf_sq$acf, lwd=2, type='l', xlim=c(0,40),
        main="", col=rainbow(m_paths)[1])  
  else
    lines(acf_sq$lag, acf_sq$acf, col=rainbow(m_paths)[i], lwd=2)
  hata=estimateAlpha(acf_sq)
  print(paste0("iteration=",i," hatAlpha=",round(hata,3),
               " ann-sigma=",round(sqrt(252)*sd(log_returns),3)))
}
n <- length(log_returns)
conf_limit <- qnorm((1 + 0.95)/2) / sqrt(n)
# Add confidence bands
abline(h = c(-conf_limit, conf_limit), col = "darkgray", lty = 2)
mtext("ACFs", side=3, line=0)
if(ToPDF) dev.off()

print(colMeans(LSSest))

print(colMeans(GHest))


if(ToPDF) pdf('fig8A.pdf',width=6,height=2.5)
par(mar=c(2,2,1,1))
window_size <- 126
JJ=which(round(t)==t)
for (i in 1:m_paths) {
  print(paste0(i," sd=",round(sd(log_returns),3),
               " kurtosis=",round(kurtosis(log_returns),3)))
  log_returns <- diff(log(paths_matrix[, i]))
  local_sd_zoo <- rollapply(
    log_returns, 
    width = window_size, 
    FUN = sd, 
    align = "right", 
    fill = NA
  )
  #if(i%%2==1)
  plot(local_sd_zoo,type="l",col=rainbow(m_paths)[i],ylim=c(0,0.04),
       xaxt = "n", xlab="")
  #else
  #  lines(local_sd_zoo,col=rainbow(m_paths)[i])
  axis(1, at = JJ, labels=t[JJ], las=2)
  
}
if(ToPDF) dev.off()
