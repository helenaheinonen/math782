library(copula)
library(VineCopula)
library(VC2copula)


# empirical cdf of the sample evaluated at the point u 
Dn <- function(u, sample){
  return(F.n(u, sample))
}

# test statistic 
Sn <- function(sample){
  Dn_vals <- apply(sample, 1, Dn, sample = sample) 
  diffs   <- Dn_vals - sample[,1] * sample[,2]
  return(sum(diffs^2))
}

# computes the value at x = (u,v) of sqrt(n)(Dn(u,v) -  uv) where Dn is the empirical cdf of the inputted sample
process <- function(x, sample, n){
  return(sqrt(n)*(Dn(x, sample) - x[1]*x[2]))
}

# computes the Rosenblatt transform of a sample w.r.t. the given copula C
R <- function(sample, C){
  return(cbind(sample[,1], dduCopula(sample, C)))
}


# Scenario functions: 
# take in a sample and returns the transformed sample adhering to the specified scenario


# scenario 1: margins and transform known
scen1 <- function(sample, C){
  return(R(sample, C)) # margins and transform known 
}


# scenario 2: margins unknown, transform known 
scen2 <- function(sample, C){
  U.hat <- pobs(sample) # margins unknown so use ranks 
  return(R(U.hat, C)) # transform known 
}

# scenario 3: margins unknown, transform known, but use ranks  
scen3 <- function(sample, C){
  U.hat <- pobs(sample) # margins unknown
  R.U.hat <- R(U.hat, C) # transform known 
  # but then use ranks 
  return(cbind(R.U.hat[,1], rank(R.U.hat[,2])/((nrow(R.U.hat))+1))) 
}

# scenario 4: margins unknown, transform unknown
scen4 <- function(sample, C){
  U.hat <- pobs(sample) # margins unknown
  fit <- fitCopula(C, U.hat, method="itau")
  theta.hat <- attributes(fit)$estimate # transform unknown 
  family <- class(C)  # identifies the class of C so we can construct C.hat
  C.hat <- do.call(family, list(theta.hat))
  return(R(U.hat, C.hat)) 
}

#scenario 5: margins unknown, transform unknown, using ranks 
scen5 <- function(sample, C){
  U.hat <- pobs(sample) # margins unknown
  fit <- fitCopula(C, U.hat, method="itau")
  theta.hat <- attributes(fit)$estimate # transform unkown 
  family <- class(C)  # identifies the class of C so we can construct C.hat
  C.hat <- do.call(family, list(theta.hat))
  R.hat.U.hat <- R(U.hat, C.hat)
  return(cbind(R.hat.U.hat[,1], rank(R.hat.U.hat[,2])/((nrow(R.hat.U.hat))+1))) # use ranks 
  
}

# transforms the data to appropriate sample based on which scenario we are in 
transform <- function(sample, C, scenario){
  switch(scenario,
         "1" = scen1(sample, C),
         "2" = scen2(sample, C),
         "3" = scen3(sample, C),
         "4" = scen4(sample, C),
         "5" = scen5(sample, C))
}
