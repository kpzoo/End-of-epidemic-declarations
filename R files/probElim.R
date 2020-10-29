######################################################################
## Compute the probability of elimination z at any time
# From: KV Parag et al., “An exact method for quantifying the reliability of end-of-epidemic 
# declarations in real time,” medRxiv, 2020.07.13.20152082, 2020.
######################################################################

# Inputs - incidence curve (Iday), window length (winLen), serial interval distribution (sidistr),
# gamma prior R hyperparameters [a b] (Rprior)
# Output - elimination probability (z), sequence of probabilities (zseq), mean R (Rz)

probElim <- function(Iday, winLen, sidistr, Rprior){
  
  # Start computing beyond epi-curve
  idst = length(Iday) + 1
  # Append incidence with zeros
  Iz = c(Iday, rep(0, 50)); nz = length(Iz)
  # Times (ids) to interrogate
  idrange = idst:(nz-1); nrange = length(idrange)
  
  # Total infectiousness (same length as Iz)
  Lamz = overall_infectivity(Iz, sidistr)
  # Epi-curve sums over window lookbacks
  B = matrix(0, 1, nrange); A = B
  
  # At each future time use window to get epi-sums
  for (i in 1:nrange) {
    # ID of current time point
    idcurr = idrange[i] 
    # Truncate window if negative time
    idback = max(idcurr - winLen + 1, 1)
    # Range of times in window
    goback = seq(idcurr, idback, -1) 
    
    # Epi-sum computation
    A[i] = sum(Lamz[goback])
    B[i] = sum(Iz[goback])
  }
  
  # Parameters of posterior gamma on R
  alpha = Rprior[1] + B
  beta = 1./(1/Rprior[2] + A)
  # Posterior mean of R and its 99% 
  Rz = alpha*beta
  
  # Estimated probability of elimination sequence
  zseq = (1 +  Lamz[idrange+1]*Rz/alpha)^(-alpha)
  # Estimated probability of elimination
  z = prod(zseq)
  
  # Output elimination data list
  probElim = list(z, as.numeric(zseq), as.numeric(Rz))
}
