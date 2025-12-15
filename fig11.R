###############MSFT) – Actual market price vs FATGBM & GBM
############################################################

# ----- Load required libraries and custom supOU / FATGBM functions -----
source("RGamma_supOU_s.R")
source("FATGBM.R")

# ----- Read market data and compute log-returns -----
# NOTE: Consider passing "SPY 1.csv" as an argument to a function for robustness.
df <- read.csv("Microsoft.csv", stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")
df <- df[order(df$Date), ]
log_returns <- diff(log(df$Price))

# ----- Define/Calculate Core Parameters -----
hatDailySigma  <-  sd(log_returns)
sigma          <-  hatDailySigma * sqrt(252) # Annualized sigma
DaysInYear     <- 252
dt             <- 1 / DaysInYear
Y              <- 1.08 # Time to maturity
S0             <- 524.11 # Initial Stock Price
r              <- 0.045 # Risk-free rate 
strike_prices  <- c(350, 360, 370,380, 390, 400, 410, 420, 430, 
                    440, 450, 460, 470, 480, 490, 500)
market_prices  <- c(193.25, 184.075, 175, 166.95,157.85, 149.67, 
                     141.37, 133.37,125.52, 118 ,110.45, 103.3, 
                     96.275, 89.57, 82.975, 76.9)
a <- 2.0
b <- a-1           
epsilon <- 0.001
InvFeVector <- getInvCdf(a, b, epsilon)
beta <- 1/40 
alpha <- 0.45
piMeasure <- setPiMeasure(type='gamma2', prm=alpha, 
                          beta=beta, T_WarmUp=2)
#dTYprm = getTY(Y=Y, InvFeVector, piMeasure)$prm
dTYprm <- c(0.39,0.094,0.47) # a=2

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
pdf("fig11.pdf", width=6, height=3)
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
mtext("MSFT option Price",side=2,line=2.3)

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
