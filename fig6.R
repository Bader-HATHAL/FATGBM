# Figure 6 in the paper
# Load required libraries 
library("tictoc")
library(latex2exp)
library(invgamma)
source("RGamma_supOU_s.R")

ToPDF=FALSE
ToPDF=TRUE

# setup of parameters 
dt <- 1/252  # Time step
T <- 30      # Time horizon
a <- 1.5       # Shape parameter for inverse gamma 
b <- a-1     # Scale parameter for inverse gamma
epsilon <- 0.001
InvFeVector <- getInvCdf(a, b, epsilon)

# Set measures for supOU simulation
B <- 40 
alpha1 <- 0.5
piMeasure1 <- setPiMeasure(type='beta', prm=alpha1, B=B, T_WarmUp=2)
alpha2 <- 0.7
piMeasure2 <- setPiMeasure(type='beta', prm=alpha2, B=B, T_WarmUp=2)

# Plotting setup
va <- c(alpha1,alpha2)
txt <- c(paste("$alpha$=",alpha1,"  "),
         paste("$alpha$=",alpha2,"  "))
vat <- c(TeX(txt))
col <- c('red','blue')
vat4 <- c(vat,"true")
col4 <- c(col,"orange")

# Simulation of supOU processes 
seedn=134
set.seed(seedn)
tic()
supOU1 <- r_supOU_smart(InvFeVector,piMeasure1,dt=dt,T=T)
x1 <- supOU1$Y
toc()

set.seed(seedn)
print('simple algorithm for supOU')
tic()
supOU2 <- r_supOU_smart(InvFeVector,piMeasure2,dt=dt,T=T)
x2 <- supOU2$Y
toc()

# Plot of simulated supOU processes
if(ToPDF) pdf('fig5.pdf',width=12,height=3)
par(mar=c(2.1, 1.1, 1.1, 1.1))
plot(supOU1$t, pmin(x1,5)+0, type='l', col=col[1],
     xlab='', ylab='', ylim=c(0,10), xaxs="i", yaxt='n')
lines(supOU2$t, pmin(x2,5)+5, col=col[2])
abline(h=c(0,5,10), col='darkgrey')
grid()
if(ToPDF) dev.off()


# Autocorrelation Function (ACF) analysis

par(mar=c(2.1, 3.1, 1.1, 1.1))
LM <- 250
r1 <- acf(x1, lag.max=LM, plot = FALSE)
r2 <- acf(x2, lag.max=LM, plot = FALSE)

Est1=estimateAlpha(r1)
print(paste("true alpha=",alpha1," Est alpha=",round(Est1,3)))
Est2=estimateAlpha(r2)
print(paste("true alpha=",alpha2," Est alpha=",round(Est2,3)))

# Theoretical ACF functions
trueacf=function(t,a) {
  1/(1+1/1*t/a*B)^a # beta
  #1/(1+1*t/beta)^a # gamma2
}
trueacfL=function(t,lambda) {
  exp(-t*lambda)
}


# Plot: ACF Comparison (Empirical vs Theoretical)

if(ToPDF) pdf('RGamma-acf.pdf',width=6,height=2)

par(mar=c(2.1, 3.1, 1.1, 1.1))
plot(r1$lag,r1$acf,type='l',col=col[3],ylim=c(0,1),xaxs="i",xlab='',ylab='',xaxt='n', bty="n")

t <- seq(0,LM,by=1)
points(t,trueacf(t*dt,va[1]),col=col[1],pch=16,cex=0.6)
lines(r2$lag,r2$acf,col=col[2])
points(t,trueacf(t*dt,va[2]),col=col[2],pch=16,cex=0.6)
lines(r1$lag,r1$acf,col=col[1])
points(t,trueacf(t*dt,va[3]),col=col[3],pch=16,cex=0.6)

legend('topright',legend=vat,col=col,lwd=2, bty = "n")
grid()
at1 <- seq(0, LM, 1)
axis(side =1, at1, labels = TRUE)
mtext('correlation function',3,-1)
if(ToPDF) dev.off()


# Density Comparison (Empirical vs Theoretical Inverse Gamma)

if(2==2)
{
  if(ToPDF) pdf('RGamma-dens.pdf',width=6,height=2)
  par(mar=c(2.1, 2.1, 1.1, 1.1))
  d1=density(x1,adjust=0.2)
  d2=density(x2,adjust=0.2)
  plot(d1,col=col[1],ylim=c(0,1.03*max(d1$y)),xlim=c(0,0+1*quantile(x1,0.98)),main='',xlab='', bty='n')
  lines(d2,col=col[2])
  vv=seq(0,6,0.01)
  pv=dinvgamma(vv,a,b)
  lines(vv,pv,col='orange')
  legend('topright',legend=vat4,col=col4,lwd=2, bty = 'n')
  mtext('density',3,-1)
  grid()
  if(ToPDF) dev.off()
}


