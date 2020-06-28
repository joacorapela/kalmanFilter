
require(ramcmc)
require(MARSS)
require(plotly)
require(mvtnorm)
require(gridExtra)
source("../src/squareRootKF.R")
source("../src/smoothLDS.R")
source("../src/estimateKFInitialCondPPCA.R")
# source("makeSymmetricMatrices.R")

processAll <- function() {
    nFactors <- 2
    simulationFilename <- "results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
    zsForPPCA <- t(as.matrix(zs))
    initialConds <- estimateKFInitialCondPPCA(z=zsForPPCA, nFactors=nFactors)

    A <- simRes$A
    C <- simRes$C
    Gamma <- simRes$Gamma
    Sigma <- simRes$Sigma
    xHat0 <- simRes$mu0
    V0 <- simRes$V0

    A0 <- initialConds$A
    # A0 <- A
    Gamma0 <- initialConds$Gamma
    C0 <- initialConds$C
    # C0 <- C
    Sigma0 <- initialConds$Sigma
    xHat00 <- xHat0
    V00 <- V0

    B1  <- matrix(list("b11", "b21", "b12", "b22"), nrow=2)
    U1  <- "unequal"
    Q1  <- "diagonal and equal"
    # Z1  <- factor(c(1,1,1,1,2,2,2,2))
    Z1  <- matrix(list("b11", "b21", "b31", "b41", "b51", "b61", "b71", "b81", "b12", "b22", "b32", "b42", "b52", "b62", "b72", "b82"), nrow=8)
    A1  <- "unequal"
    R1  <- "diagonal and equal"
    pi1 <- "unequal"
    V01 <- "zero"

    model.list <- list(B=B1, U=U1, Q=Q1, Z=Z1, A=A1, R=R1, x0=pi1, V0=V01)
    inits <- list(B=A0, Q=Gamma0, Z=C0, R=Sigma0, x0=xHat00, V0=V00)
    # kem = MARSS(zs, model=model.list)
    # kem = MARSS(zs, model=model.list, inits=inits)
    kem <- MARSS(zs, model=model.list, inits=list(x0=0), silent=2)

    browser()
}

processAll()
