# Figure 7 in the paper
# Load required libraries
source("RGamma_supOU_s.R") 

ToPDF <- FALSE
ToPDF <- TRUE

# setup of parameters 
DaysInYear <- 252
a <- 1.5
b <- a - 1
epsilon <- 0.001
InvFeVector <- getInvCdf(a, b, epsilon)

# Set the measure for supOU simulation
beta <- 1/40 
alpha <- 0.6
piMeasure1=setPiMeasure(type='gamma2', prm=alpha, beta=beta, T_WarmUp=2)

n_paths <- 10  # Number of realizations to simulate and plot

# Function: Simulate and plot realizations of T_t over [0, Y]

plot_Tt_realizations <- function(Y, title_suffix) {
  dt <- 1 / DaysInYear
  t <- seq(0, Y, by = dt)
  N <- length(t)
  Tt_matrix <- matrix(0, nrow = N, ncol = n_paths)
  for (j in 1:n_paths) {
    ou <- r_supOU(InvFeVector, piMeasure1, dt = dt, T = Y)
    Tt_matrix[, j] <- cumsum(dt * ou$Y)
  }
  plot(t, Tt_matrix[, 1], type = "l", col = 1, lwd = 2,
       ylim = range(Tt_matrix), bty="n",
       xlab = "Time ", ylab = "Fractal activity time")
  grid()
  for (j in 2:n_paths) {
    lines(t, Tt_matrix[, j], col = j, lwd = 2)
  }
}

if(ToPDF) pdf('fig6.pdf',width=5,height=4)
par(mar=c(2.1, 2.1, 1.1, 1.1))
seedn <- 154
set.seed(seedn)
plot_Tt_realizations(Y = 1, title_suffix = "(One Year)")
plot_Tt_realizations(Y = 10, title_suffix = "(Ten Years)")
if(ToPDF) dev.off()



