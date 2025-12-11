library(copula)
library(VineCopula)
library(VC2copula)

data(uranium)
data <- uranium[,c(7, 3, 6)] # Titanium, Cobalt, Scandium 
head(data)
plot(data)
Uhat <- pobs(data)
pairs(Uhat)

fit <- fitCopula(normalCopula(dim = 3), Uhat, method = "mpl")
fit2 <- fitCopula(tCopula(dim=3, df = 8, df.fixed=TRUE), Uhat, method="mpl")
theta_hat <- coef(fit2)
Ctheta <- tCopula(param = theta_hat, df=8, dim = 3, df.fixed=TRUE)
Ctheta2 <- normalCopula(param=theta_hat, dim=3)
Ru <- pairwiseCcop(Uhat, Ctheta)
pairsRosenblatt(Ru, method="scatter")
pairsRosenblatt(Ru, method="QQchisq")


# generate from 3 dimensional clayton

C_true <- claytonCopula(param = 3, dim = 3)
n <- 500
U <- rCopula(n, C_true)

pairs(U, main = "True sample from 3D Clayton copula")

fit_wrong <- fitCopula(normalCopula(dim = 3), U, method = "mpl")
theta_wrong <- coef(fit_wrong)
C_wrong <- normalCopula(param = theta_wrong, dim = 3)

Ru <- pairwiseCcop(U, C_wrong)
pairsRosenblatt(Ru, method="scatter")
pairsRosenblatt(Ru, method="QQchisq")



# what about using a vine model?

# vine tree structure matrix 
Matrix <- c(3, 1, 2,
            0, 1, 2,
            0, 0, 2)
Matrix <- matrix(Matrix, 3, 3)

# copula family matrix
family <- c(0, 1, 3,   # 0 = independence, 1 = Gaussian, 3 = Clayton
            0, 0, 1,
            0, 0, 0)
family <- matrix(family, 3, 3)

#  parameter matrix
par <- c(0, 0.5, 1.2,
         0, 0, 0.8,
         0, 0, 0)
par <- matrix(par, 3, 3)

par2 <- matrix(0, 3, 3)

RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3"))

# simulate a sample of size 300 from the R-vine copula model
set.seed(123)
simdata <- RVineSim(300, RVM)
pairs(simdata)

fitNormal <- fitCopula(normalCopula(dim = 3), simdata, method = "mpl")
theta_hat <- coef(fitNormal)
Ctheta <- normalCopula(param=theta_hat, dim=3)
Ru <- pairwiseCcop(simdata, Ctheta)
pairsRosenblatt(Ru, method="scatter")
pairsRosenblatt(Ru, method="QQchisq")
