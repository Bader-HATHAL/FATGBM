# Load required libraries
source("FATGBM.R")

# Parameter setup
DaysInYear=252
S0 =100    
K = 110
B = 145 
Y = 1 
r = 0.04  
sigma = 0.3 
dt <- 1 / DaysInYear

# Set the measure for supOU simulation
beta <- 1/40 
alpha <- 0.5 # 0.5
piMeasure <- setPiMeasure(type='gamma2', prm=alpha, 
                           beta=beta, T_WarmUp=2)
a <- 1.5 # 1.5
b <- a - 1
epsilon <- 0.001
InvFeVector <- getInvCdf(a, b, epsilon)
#dTYprm <- getTY(Y=1, InvFeVector, piMeasure)$prm
dTYprm <- c(0.26,0.061,0.45)

for(K in c(105,110,115))
{
  set.seed(139)
  rs=simulate_BarrierOption(S0, K, B, r, sigma, Y, dt, InvFeVector, piMeasure, m=2000) 
  FATUpOutCall(Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
  FATUpOutPut(Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
}

B=110
for(K in c(115,120,125))
  {
  set.seed(139)
  rs=simulate_BarrierOption(S0, K, B, r, sigma, Y, dt, InvFeVector, piMeasure, m=2000) 
  FATUpOutCall(Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
  FATUpOutPut(Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
}


# Option Price Comparison
#message("B=",B," K=",K)
#FATUpOutCall( Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
#GBMUpOutCall( Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, echo=TRUE)
#FATUpOutPut( Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
#GBMUpOutPut( Y=Y, S0=S0, K=K, B=B, sigma=sigma, r=r, echo=TRUE)
#FATEurCall( Y=Y, S0=S0, K=K, sigma=sigma, r=r, dTYprm=dTYprm, echo=TRUE)
#GBMEurCall( Y=Y, S0=S0, K=K, sigma=sigma, r=r, echo=TRUE)


