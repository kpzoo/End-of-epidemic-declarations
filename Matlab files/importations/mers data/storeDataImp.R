######################################################################
## Write data from EpiEstim and outbreaks packages
######################################################################

# Assumptions and modifications
# - save MERS data for end-of-epidemic import example
# - need example for under-reporting

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Main packages
library("EpiEstim")
library("caTools")
library("outbreaks")

# Boolean for plotting
wantPlot = 0

# Load data on MERS from 2014-2015
data("mers_2014_15")
Imers = mers_2014_15$incidence; genVals = mers_2014_15$si
# Ensure generation time long enough for appended zeros in Matlab
genmers = discr_si(0:1000, genVals$mean_si, genVals$std_si)
# Separate local vs imported
Iloc = Imers$local; Iimp = Imers$imported

# Append false zeros
nz = 100
Iloc = c(Iloc, rep(0, nz))
Iimp = c(Iimp, rep(0, nz))

# Total infectiousness from local + imported
Itot = Iloc + Iimp
Lmers = overall_infectivity(Itot, genmers)
Lmers[is.na(Lmers)] = 0

# Local and imported data storage
write.csv(Iloc, file = "Iloc.csv")
write.csv(Iimp, file = "Iimp.csv")
write.csv(Lmers, file = "Lmers.csv")
write.csv(genmers, file = "genmers.csv")