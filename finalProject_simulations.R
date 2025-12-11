library(copula)
library(VineCopula)
library(VC2copula)
source("simulationsSourceFile.R")
library(simsalapar)




x_points <- list(
  list(c(0.5, 0.5)),
  list(c(1, 0.5))
  )

varList <-
  varlist( 
    n.sim = list(type="N", expr = quote(m), value = 1000), # how many iterations we run 
    n = list(type="grid", expr = quote(n), value = c(250, 500)), # sample size 
    C = list(type="grid", expr = quote(C), value = c("Clayton","Gumbel")), # copula families 
    x = list(type="frozen", value = c(0.5,0.5)), # compute process at this point 
    tau = list(type="grid", expr = quote(tau), value = c(0.2, 0.5)), # different amounts of dependence 
    scenario = list(type="inner", expr = quote(scenario), value =c("1","2","3","4","5")) # our 5 scenarios (margins/transform known/unknown)
  )


varList2 <-
  varlist( 
    n.sim = list(type="N", expr = quote(m), value = 1000), # how many iterations we run 
    n = list(type="grid", expr = quote(n), value = c(250, 500)), # sample size 
    C = list(type="grid", expr = quote(C), value = c("Clayton","Gumbel")), # copula families 
    u = list(type="grid", expr = quote(u), value = c(0, 0.5, 1)), 
    v = list(type="grid", expr = quote(v), value = c(0, 0.5, 1)), # so we'll compute process point at all combos of (u,v)
    tau = list(type="frozen", value = 0.5), # different amounts of dependence 
    scenario = list(type="inner", expr = quote(scenario), value =c("1","2","3","4","5")) # our 5 scenarios (margins/transform known/unknown)
  )

varList3 <-
  varlist( 
    n.sim = list(type="N", expr = quote(m), value = 1000), # how many iterations we run 
    n = list(type="grid", expr = quote(n), value = c(250, 500)), # sample size 
    C = list(type="grid", expr = quote(C), value = c("Clayton","Gumbel")), # copula families 
    tau = list(type="frozen", value = 0.5), # different amounts of dependence 
    scenario = list(type="inner", expr = quote(scenario), value =c("1","2","3","4","5")) # our 5 scenarios (margins/transform known/unknown)
  )


# simsalapar will call doOne at each grid point and pass scenario in at each iteration(bc "inner")
# doOne spits out a numeric vector of length 5 with the process value for each of the 5 scenarios 

doOne <- function(n, C, u, v, tau, scenario){
  x <- c(u,v)
  cop <- switch(C,
                "Clayton" = claytonCopula(iTau(claytonCopula(),tau)),
                "Gumbel" = gumbelCopula(iTau(gumbelCopula(),tau)))
  U <- rCopula(n, cop) # simulates the original sample 
  
  results <- sapply(scenario, function(sc) {
    transformed <- transform(U, cop, sc)
    process(x, transformed, n)
  })
  results
}

doOne_Sn <- function(n, C, tau, scenario){
  cop <- switch(C,
                "Clayton" = claytonCopula(iTau(claytonCopula(),tau)),
                "Gumbel" = gumbelCopula(iTau(gumbelCopula(),tau)))
  U <- rCopula(n, cop) # simulates the original sample 
  
  results <- sapply(scenario, function(sc) {
    transformed <- transform(U, cop, sc)
    Sn(transformed)
  })
  results
}

# testing a single grid point
# hard coded scenario 1 (so we can compare)

set.seed(13)
manualTest <- numeric(1000)
cop <- claytonCopula(iTau(claytonCopula(),0.5))
for(i in 1:1000){
  set.seed(i)
  U <- rCopula(500, cop) # sample of size 500 from Clayton 
  R.U <- R(U, cop) # Rosenblatt transform 
  manualTest[i] <- process(x = c(0.5,0.5), R.U, n = 500) 
}
hist(manualTest)

test1 <- matrix(NA, nrow = 1000, ncol = 5)
for(i in 1:1000){
  set.seed(i)
  test1[i, ] <- doOne( # make each row be an iteration of doOne
    n = 500,
    C = "Clayton",
    u=0.5,
    v=0.5,
    tau = 0.5,
    scenario = c("1","2","3","4","5")
  )
}
colnames(test1) <- paste0("scenario", 1:5)
# for scenario 1 we need the first column 
hist(test1[,1]) 

manualSntest <- matrix(NA, nrow = 1000, ncol = 5)
for(i in 1:1000){
  set.seed(i)
  manualSntest[i, ] <- doOne_Sn(
    n=500,
    C = "Clayton",
    tau = 0.5,
    scenario = c("1", "2", "3", "4", "5")
  )
}
boxplot(manualSntest[,1])

# run simulation 
res <- doLapply(varList, doOne = doOne)
val <- getArray(res) # all the values 

# trying multiple x points
res_multipleX <- doLapply(varList2, doOne=doOne)
val_multipleX <- getArray(res_multipleX)

res_multipleX_itau <- doLapply(varList2, doOne=doOne)
val_multipleX_itau <- getArray(res_multipleX_itau)


# simulation for Sn values
res_Sn <- doLapply(varList3, doOne=doOne_Sn)
val_Sn <- getArray(res_Sn)

save(val, file = "processPoints_simulation1.Rdata") # from my first round of simulations 
save(val_multipleX, file="processPoints_simulation2.Rdata") # includes multiple x points
save(val_multipleX_itau, file="processPoints_simulation2_itau.Rdata")
save(val_Sn, file="Sn_simulation.Rdata")


# extracting a single grid point 
vals_scen1 <- val[
  1, # scenario 
  which(varList$n$value == 500), 
  which(varList$C$value=="Clayton"),
  which(varList$tau$value==0.5),
  1:1000
]
hist(vals_scen1) 
sd(vals_scen1)


vals_scen2 <- val[
  2, # scenario 
  which(varList$n$value == 500), 
  which(varList$C$value=="Clayton"),
  which(varList$tau$value==0.5),
  1:1000
]

vals_scen3 <- val[
  3, # scenario 
  which(varList$n$value == 500), 
  which(varList$C$value=="Clayton"),
  which(varList$tau$value==0.5),
  1:1000
]

vals_scen4 <- val[
  4, # scenario 
  which(varList$n$value == 500), 
  which(varList$C$value=="Clayton"),
  which(varList$tau$value==0.5),
  1:1000
]

vals_scen5 <- val[
  5, # scenario 
  which(varList$n$value == 500), 
  which(varList$C$value=="Clayton"),
  which(varList$tau$value==0.5),
  1:1000
]
hist(vals_scen5)
sd(vals_scen5)


# extracting Sn values 
Sn_scen1 <- val_Sn[
  1,
  which(varList3$n$value==500),
  which(varList3$C$value=="Clayton"),
  1:1000
]
boxplot(Sn_scen1)

Sn_scen2 <- val_Sn[
  2,
  which(varList3$n$value==500),
  which(varList3$C$value=="Clayton"),
  1:1000
]

Sn_scen3 <- val_Sn[
  3,
  which(varList3$n$value==500),
  which(varList3$C$value=="Clayton"),
  1:1000
]

Sn_scen4 <- val_Sn[
  4,
  which(varList3$n$value==500),
  which(varList3$C$value=="Clayton"),
  1:1000
]

Sn_scen5 <- val_Sn[
  5,
  which(varList3$n$value==500),
  which(varList3$C$value=="Clayton"),
  1:1000
]

boxplot(Sn_scen1)
boxplot(Sn_scen2)
boxplot(Sn_scen3)
boxplot(Sn_scen4)
boxplot(Sn_scen5)

clayton500_0.5_0.5 <- list(vals_scen1, vals_scen2, vals_scen3, vals_scen4, vals_scen5) 
# for n = 500, Clayton, process values at x = (0.5, 0.5)
stderrs_0.5_0.5 <- numeric(length=5)
for (i in 1:5){
  scen <- clayton500_0.5_0.5[[i]]
  stderrs_0.5_0.5[i] <- sd(scen)
} # why does std error go down with the estimation ? 


#investigate the std error from the histograms for different x points 

# extract it from the multiple x plot 
scen1_0.5_1 <- val_multipleX[
  1,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==1),
  1:1000
] 

scen1_1_0.5 <- val_multipleX[
  1,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==1),
  which(varList2$u$value==0.5),
  1:1000
] 


scen2_0.5_1 <- val_multipleX[
  2,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==1),
  1:1000
] 

scen2_1_0.5 <- val_multipleX[
  2,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==1),
  which(varList2$u$value==0.5),
  1:1000
] 

scen3_0.5_1 <- val_multipleX[
  3,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==1),
  1:1000
] 

scen3_1_0.5 <- val_multipleX[
  3,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==1),
  which(varList2$u$value==0.5),
  1:1000
] 


scen4_0.5_1 <- val_multipleX[
  4,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==1),
  1:1000
] 

scen4_1_0.5 <- val_multipleX[
  4,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==1),
  which(varList2$u$value==0.5),
  1:1000
] 


scen5_0.5_1 <- val_multipleX[
  5,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==1),
  1:1000
] 

scen5_1_0.5<- val_multipleX[
  2,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==1),
  which(varList2$u$value==0.5),
  1:1000
] 

scen1_itau_0.5_0.5 <- val_multipleX_itau[
  1,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==0.5),
  1:1000
]

scen2_itau_0.5_0.5 <- val_multipleX_itau[
  2,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==0.5),
  1:1000
]


scen3_itau_0.5_0.5 <- val_multipleX_itau[
  3,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==0.5),
  1:1000
]


scen4_itau_0.5_0.5 <- val_multipleX_itau[
  4,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==0.5),
  1:1000
]


scen5_itau_0.5_0.5 <- val_multipleX_itau[
  5,
  which(varList2$n$value==500),
  which(varList2$C$value=="Clayton"),
  which(varList2$u$value==0.5),
  which(varList2$u$value==0.5),
  1:1000
]

out <- list(sd(scen1_itau_0.5_0.5),sd(scen2_itau_0.5_0.5), sd(scen3_itau_0.5_0.5),sd(scen4_itau_0.5_0.5),sd(scen5_itau_0.5_0.5))

clayton500_1_0.5 <- list(scen1_1_0.5, scen2_1_0.5, scen3_1_0.5, scen4_1_0.5, scen5_1_0.5)
stderrs_1_0.5 <- numeric(length=5)
for (i in 1:5){
  scen <- clayton500_1_0.5[[i]]
  stderrs_1_0.5[i] <- sd(scen) / sqrt(length(scen)) 
} 

















# old cold
# scenario 1 

set.seed(123)
process_vals1 <- numeric(N)
Sn_vals1 <- numeric(N)
for(i in 1:N){
  U <- rCopula(n, C_known) # sample of size 500 from Clayton theta = 2  
  R.U <- cbind(U[,1], dduCopula(U, C_known)) # rosenblatt transform 
  process_vals1[i] <- process(x=x, R.U)
}
hist(process_vals1)
for(i in 1:N){
  U <- rCopula(n, C_known) # sample of size 500 from Clayton theta = 2  
  R.U <- cbind(U[,1], dduCopula(U, C_known)) # rosenblatt transform 
  Sn_vals1[i] <- Sn(R.U)
}


# scenario 2 and 3 

process_vals2 <- numeric(N)
Sn_vals2 <- numeric(N)

process_vals3 <- numeric(N)
Sn_vals3 <- numeric(N)

for(i in 1:N){
  U <- rCopula(n, C_known)
  U.hat <- pobs(U) # use scaled component wise ranks 
  R.U.hat <- cbind(U.hat[,1], dduCopula(U.hat, C_known)) # rosenblatt transform (known)
  process_vals2[i] <- process(x=x, R.U.hat)
  # can do scenario 3 in the same loop 
  R.U.hat.ranks <- cbind(R.U.hat[,1], rank(R.U.hat[,2])/(n+1))
  # the sample we are now dealing with is (\hat U_i, S_i) where S_i is the normalized rank
  # of the R(\hat V_i ; \hat U_i) from among R(\hat V_1 ; \hat U_1), … R(\hat V_n ; \hat U_n)
  process_vals3[i] <- process(x=x, R.U.hat.ranks)
}
for(i in 1:N){
  U <- rCopula(n, C_known)
  U.hat <- pobs(U) # use scaled component wise ranks 
  R.U.hat <- cbind(U.hat[,1], dduCopula(U.hat, C_known)) # rosenblatt transform 
  Sn_vals2[i] <- Sn(R.U.hat)
}
hist(process_vals2)

# for scenario 3 they are all the same value? 
# Either I interpreted your instructions wrong, or they don't model the situation appropriately?
hist(process_vals3)     


# scenario 4 and 5 

process_vals4 <- numeric(N)
process_vals5 <- numeric(N)

for(i in 1:N){
  U <- rCopula(n, C_known)
  U.hat <- pobs(U)
  # estimate theta
  fit <- fitCopula(claytonCopula(), U.hat, method="mpl")
  theta.hat <- attributes(fit)$estimate
  # now we need to do rosenblatt transform using this estimated theta 
  C.hat <- claytonCopula(theta.hat)
  R.hat.U.hat <- cbind(U.hat[,1], dduCopula(U.hat, C.hat)) # rosenblatt using the estimated theta
  process_vals4[i] <- process(x = x, R.hat.U.hat)
  R.hat.U.hat.ranks <- cbind(R.U.hat[,1], rank(R.hat.U.hat[,2])/(n+1)) # normalized ranks on the second column
  process_vals5[i] <- process(x = x, R.hat.U.hat.ranks)
}
hist(process_vals4)
hist(process_vals5) 





