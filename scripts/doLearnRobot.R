
library(R.matlab)
library(mvtnorm)
source("../src/filterLDS.R")
source("../src/smoothLDS.R")
source("../src/emEstimationKF.R")
source("makeSymmetricMatrices.R")

processAll <- function() {
    initialCondNoiseSD <- 0.00001
    nIter <- 10
    simulationFilename <- "~/dev/research/programs/github/python/pykalman/pykalman/datasets/data/robot.mat"

    simRes <- readMat(simulationFilename)
    x <- simRes$y
    A <- simRes$A
    Gamma <- 10*diag(5)
    C <- simRes$C
    Sigma <- 10*diag(2)
    mu0 <- simRes$x0
    V0 <- simRes$V.0

    A0 <- A
    Gamma0 <- Gamma
    C0 <- C
    Sigma0 <- Sigma
    mu00 <- mu0
    V00 <- V0

    res <- emEstimationKF(x=x, A0=A0, Gamma0=Gamma0, C0=C0, Sigma0=Sigma0, mu00=mu00, V00=V00, nIter=nIter)
    plot(unlist(res$logLikelihoods), type="b")
    browser()
}

processAll()
