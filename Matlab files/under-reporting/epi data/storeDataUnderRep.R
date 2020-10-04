######################################################################
## Write data from EpiEstim and outbreaks packages
######################################################################

# Assumptions and modifications
# - save flu/SARS data for end-of-epidemic under-reporting example

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

# Boolean for epidemic choice
isFlu = 0
# Append false zeros
nz = 100

if(isFlu){
  # Load data on Flu 1918
  data("Flu1918"); Iday = Flu1918$incidence
  # Serial interval extended to more days
  sidist = Flu1918$si_distr
  
} else{
  # Load data on SARS2003
  data("SARS2003"); Iday = SARS2003$incidence
  # Serial interval extended to more days
  sidist = SARS2003$si_distr
}

# Append zeros to incidence and SI
Iday = c(Iday, rep(0, nz))
sidist = c(sidist, rep(0, 300))

# Total infectiousness from local cases
Lday = overall_infectivity(Iday, sidist)
Lday[is.na(Lday)] = 0

# Local and imported data storage
write.csv(Iday, file = "Iday.csv")
write.csv(Lday, file = "Lday.csv")
write.csv(sidist, file = "sidist.csv")