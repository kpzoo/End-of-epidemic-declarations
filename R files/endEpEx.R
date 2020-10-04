######################################################################
## Compute end-of-epidemic declaration time from empirical data
# From: KV Parag et al., “An exact method for quantifying the reliability of end-of-epidemic 
# declarations in real time,” medRxiv, 2020.07.13.20152082, 2020.
######################################################################

# Method uses the serial interval distribution incidence curve to compute z,
# the probability that the epidemic has been eliminated at a given time
# Currently assumes perfect incidence - for practical apply corrections to
# overcomes biases in the incidence curve (e.g. backpropagation for delays)

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Packages for renewal model analysis
library(EpiEstim)
# Source all files in current folder
source('probElim.R')

# Compute over entire epidemic or from last 0
wholeEp = 1

# EpiEstim data for flu and SARS
data("Flu2009"); data("SARS2003")

# Choose data source 
#Iday = Flu2009$incidence[,2]; sidistr = Flu2009$si_distr
Iday = SARS2003$incidence; sidistr = SARS2003$si_distr

# Original data
ndayOrig = length(Iday); IdayOrig = Iday

# Artificial zeros to guarantee end links to SI support
suppSI = max(which(sidistr > 0))
Iday = c(Iday, rep(0, suppSI)); nday = length(Iday)

# Time of last case 
tlast = max(which(Iday > 0)) 
# Query times for 'is it over?' - set to 
if(wholeEp){
  tquery = 1:nday
}else{
  tquery = (tlast+1):nday
}
nquery = length(tquery)

# Set long window and prior on R
Rprior = c(100, 5); winLen = ndayOrig-2

# For each day beyond the last case obtain z
z = rep(0, nquery)
for(i in 1:nquery){
  # Elimination probability 
  zdata = probElim(Iday[1:tquery[i]], winLen, sidistr, Rprior)
  z[i] = zdata[[1]]
}
# Declaration time
tdec = min(which(z >= 0.95))
print(paste0(c('declaration time is ', tdec), collapse = ''))

# Plot incidence curve and elimination probability
if(wholeEp){
  quartz()
  par(mfrow=c(2,1))
  # Incidence curve
  plot(1:nday, Iday, type = 'h', bty = 'l', lwd = 2, col='red',
       xlab = 'time (days)', ylab = 'incidence')
  # Elimination probabilities over whole curve
  plot(1:nday, z, type = 'h', bty = 'l', lwd = 2, col = 'blue',
       xlab = 'time (days)', ylab = 'P(elimination)')
}