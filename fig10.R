############################################################
## S&P 500 (SPY) – Actual market price vs FATGBM & GBM
############################################################

# ----- Load required libraries and custom supOU / FATGBM functions -----
source("RGamma_supOU_s.R")
source("FATGBM.R")

# ----- Read market data and compute log-returns -----
# NOTE: Consider passing "SPY 1.csv" as an argument to a function for robustness.
df <- read.csv("SPY 1.csv", stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
df <- df[order(df$Date), ]
log_returns <- diff(log(df$Price))

# ----- Define/Calculate Core Parameters -----
hatDailySigma  <-  sd(log_returns)
sigma          <-  hatDailySigma * sqrt(252) # Annualized sigma
DaysInYear     <- 252
dt             <- 1 / DaysInYear
Y              <- 0.92 # Time to maturity
S0             <- 636.94 # Initial Stock Price
r              <- 0.045 # Risk-free rate 
strike_prices  <- c(540, 545, 550, 555, 560, 565, 570, 575, 580, 
                    585, 590, 595, 600, 605, 610, 615, 620, 625)
market_prices  <- c(124.61, 120.37, 116.17, 111.99, 107.85, 103.51, 
                    99.64, 95.81, 91.81, 87.69, 83.84, 79.92, 76.10,
                    72.28, 68.45, 64.78, 61.12, 57.53)
a <- 1.5
b <- a-1           
epsilon <- 0.001
InvFeVector <- getInvCdf(a, b, epsilon)
beta <- 1/40 
alpha <- 0.55
piMeasure <- setPiMeasure(type='gamma2', prm=alpha, 
                          beta=beta, T_WarmUp=2)
#dTYprm = getTY(Y=Y, InvFeVector, piMeasure)$prm
dTYprm <- c(0.23,0.050,0.46) # a=1.5

# ----- Calculate Option Prices -----
# Set echo=FALSE for final plots
fatgbm_prices <- sapply(strike_prices, function(K) {
  FATEurCall(Y=Y, S0=S0, K=K, sigma=sigma, r=r, dTYprm=dTYprm, echo=FALSE)
})
gbm_prices    <- sapply(strike_prices, function(K) {
  GBMEurCall(Y=Y, S0=S0, K=K, sigma=sigma, r=r, echo = FALSE)
})

# ----- Prepare Tidy Data for Plotting -----
plot_data <- data.frame(
  Strike = strike_prices,
  Market = market_prices,
  FATGBM = fatgbm_prices,
  GBM    = gbm_prices
)

# ----- Generate Plot (Base R - Corrected) -----
pdf("fig10.pdf", width=6, height=3)
par(mar=c(3, 3.2, 1, 1)) # Adjusted margins for better label fit

# Dynamically set y-limits based on all price data
y_range <- range(c(plot_data$Market, plot_data$FATGBM, plot_data$GBM), na.rm = TRUE)

plot(plot_data$Strike, plot_data$Market, 
     type = "p", 
     pch = 19, 
     col = "green",
     cex = 0.8,
     # Use dynamic limits and include labels in the main call
     xlim = range(plot_data$Strike, na.rm = TRUE),        
     ylim = y_range,         
     xlab = "", ylab = "", main = "",
     las = 1, # Horizontal y-axis labels
     bty = "l" # L-shaped plot frame
) 
mtext("Strike Price (K)",side=1,line=2)
mtext("SPY option Price",side=2,line=2.3)

# Add model points
points(plot_data$Strike, plot_data$FATGBM, pch = 3, col = "black", cex = 0.8)
points(plot_data$Strike, plot_data$GBM , pch = 4, col = "black", cex = 0.8)

# Add Legend
legend("topright", 
       legend = c("Market Price", "FATGBM Model", "GBM Model"),
       col = c("green", "black", "black"),
       pch = c(19, 3, 4),
       cex = 0.8,
       bty = "n")
grid()
dev.off() # Close the PDF device to save the file

# ----- Compute RMSE to evaluate pricing errors for each model -----
rmse_fat = sqrt(mean((fatgbm_prices - market_prices)^2))
rmse_gbm = sqrt(mean((gbm_prices  - market_prices)^2))
cat("RMSE (FATGBM):", round(rmse_fat, 4), "\n")
cat("RMSE (GBM):", round(rmse_gbm, 4), "\n")
