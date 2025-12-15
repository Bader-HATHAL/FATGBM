# Load required libraries
source("FATGBM.R")

# Set the measure for supOU simulation
beta <- 1/40 
alpha <- 0.5
piMeasure5 <- setPiMeasure(type='gamma2', prm=alpha, 
                           beta=beta, T_WarmUp=2)
piMeasure6 <- setPiMeasure(type='gamma2', prm=0.6, 
                           beta=beta, T_WarmUp=2)
piMeasure7 <- setPiMeasure(type='gamma2', prm=0.7, 
                           beta=beta, T_WarmUp=2)

epsilon <- 0.001
for(nu in c(1.5,2,2.5,3)*2)
{
  a <- nu/2
  b <- a - 1
  InvFeVector <- getInvCdf(a, b, epsilon)
  res <- getTY(Y=1, InvFeVector, piMeasure5, m=2000)
  res <- getTY(Y=1, InvFeVector, piMeasure6, m=2000)
  res <- getTY(Y=1, InvFeVector, piMeasure7, m=2000)
}
