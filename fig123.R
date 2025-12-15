# Load required libraries
library(MASS)
library(ghyp)
library(fitdistrplus)
library(lubridate)
library(moments)
library(zoo)
source("RGamma_supOU_s.R")

oname=c("SPY","MSFT")[1]
# Read the data
if(oname=="MSFT") 
  df <- read.csv("Microsoft.csv", stringsAsFactors = FALSE)
if(oname=="SPY") 
  df <- read.csv("SPY 1.csv", stringsAsFactors = FALSE)

# Convert Date and sortdf$DateS=df$Date
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
df <- df[order(df$Date), ]

# Compute log-returns
log_returns <- diff(log(df$Price))

sample_skewness = skewness(log_returns)
sample_kurtosis = kurtosis(log_returns)
hat_daily_sigma  =  sd(log_returns)

# Fit t-univariate distribution
fitted_t <- fit.tuv(data = log_returns, symmetric = TRUE, silent=TRUE)
print(fitted_t)

# Fit GHyperbolic distribution
fitted_ghyp <- fit.ghypuv(data = log_returns, mu = mean(log_returns), 
                          sigma = sd(log_returns), silent=TRUE)
print(fitted_ghyp)

# Fit normal distribution
fit_norm <- fitdistr(log_returns, "normal")

# Fit Student's t-distribution (another parametrization)
fit_t <- fitdistr(log_returns, densfun = "t", 
                  start = list(df = 5, m = mean(log_returns), 
                               s = sd(log_returns)))

# Create density values
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
norm_density <- dnorm(x_vals, mean = fit_norm$estimate["mean"], 
                      sd = fit_norm$estimate["sd"])
t_density <- dt((x_vals - fit_t$estimate["m"]) / fit_t$estimate["s"], 
                df = fit_t$estimate["df"]) / fit_t$estimate["s"]

pdf(paste0("fig1",oname,".pdf"),width=6,height=2.5)
par(mar=c(2,2,1,1))
# Plot 1A
hist(log_returns, 
     breaks = seq(1.1*min(log_returns),1.1*max(log_returns),0.0025), 
     freq = FALSE, col = "lightgray", main="",
     xlab = "Log-Return", ylab = "Density",
     xlim=c(-0.05,0.05))
mtext(paste0("Densities of Log-Returns for ",oname),line=0,side=3)
lines(x_vals, norm_density, col = "red", lwd = 2)
lines(x_vals, t_density, col = "blue", lwd = 2, lty = 1)
legend("topright", legend = c("Normal Fit", "Student's t Fit"),
       col = c("red", "blue"), lwd = 2, lty = c(1, 1))
legend("topleft", bty="n", legend = c(
  paste0("kurtosis=",round(sample_kurtosis,2)), 
  paste0("sd=",round(hat_daily_sigma,4))))

# Plot 1B
ed=density(log_returns, adjust = 0.25)
plot(ed, type='l', lwd=2, log='y', bty="n",
     xlim=c(-0.05,0.05), ylim=c(0.1,1.4*max(ed$y)), main="")
mtext(paste0("Log-densities of Log-Returns for ",oname),line=0,side=3)
lines(x_vals, norm_density, col = "red", lwd = 2)
lines(x_vals, t_density, col = "blue", lwd = 2, lty = 1)
legend("topright", legend = c("Normal Fit", "Student's t Fit"),
       col = c("red", "blue"), lwd = 2, lty = c(1, 1))

# Compute absolute and squared log-returns
abs_log_returns <- abs(log_returns)
squared_log_returns <- log_returns^2
# Compute ACFs
acf_log <- acf(log_returns, plot = FALSE, lag.max = 40)
acf_abs <- acf(abs_log_returns, plot = FALSE, lag.max = 40)
acf_sq <- acf(squared_log_returns, plot = FALSE, lag.max = 40)
# Plot 2
plot(acf_log$lag, acf_log$acf, type = "l", col = "black", lwd = 2,
     xlab = "", ylab = "", main = "", bty="n",
     ylim = range(c(acf_log$acf, acf_abs$acf, acf_sq$acf)))
grid()
lines(acf_abs$lag, acf_abs$acf, col = "blue", lwd = 2, lty = 1)
lines(acf_sq$lag, acf_sq$acf, col = "red", lwd = 2, lty = 1)
# Confidence interval (approximate for white noise)
n <- length(log_returns)
conf_limit <- qnorm((1 + 0.95)/2) / sqrt(n)
# Add confidence bands
abline(h = c(-conf_limit, conf_limit), col = "darkgray", lty = 2)
legend("topright", 
       legend = c("Log-returns", "Absolute log-returns", 
                  "Squared log-returns"),
       col = c("black", "blue", "red"), lwd = 2, lty = c(1, 1, 1))
mtext(paste0("ACFs for ",oname),side=3,line=0)


# Plot 3
hat_alpha <- estimateAlpha(acf_sq)
mtext(paste0("log-log ACF for ",oname),side=3,line=0)
mtext(paste0("fit: alpha=",round(hat_alpha,3)),side=3,line=-1.5)
hat_H  <-  1 / (1 + hat_alpha)
print(paste("alpha,H=",hat_alpha,hat_H))

# plot the time series of Price
par(mar=c(3,1,1,1))
plot(df$Price, type='l', xaxt = "n", xlab="")
m <- format(df$Date, "%m")
d <- format(df$Date, "%d")
y <- format(df$Date, "%Y")
J2 <- which( m == "01" & d == "03")
JJ <- J2[1]+(0:19)*252
axis(1, at = JJ, labels=y[JJ], las=2)
mtext(oname,line=0,side=3)

window_size <- 126
local_sd_zoo <- rollapply(
  log_returns, 
  width = window_size, 
  FUN = sd, 
  align = "right", 
  fill = NA
)
plot(local_sd_zoo, type="l", ylim=c(0,max(local_sd_zoo,na.rm=T)),
     xaxt = "n", xlab="")
#abline(h=sd(log_returns),col="blue")
axis(1, at = JJ, labels=y[JJ], las=2)
mtext(paste0("rolling volatility for ",oname),line=0,side=3)

dev.off()
