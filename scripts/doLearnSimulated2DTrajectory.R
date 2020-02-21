
library(mvtnorm)
source("../src/filterLDS.R")
source("../src/smoothLDS.R")
source("../src/emEstimationKF.R")
source("makeSymmetricMatrices.R")

processAll <- function() {
    initialCondNoiseSD <- 0.00001
    nIter <- 100
    simulationFilename <- "results/simulation.RData"

    simRes <- get(load(simulationFilename))
    x <- simRes$x
    z <- simRes$z
    A <- simRes$A
    Gamma <- simRes$Gamma
    C <- simRes$C
    Sigma <- simRes$Sigma
    mu0 <- simRes$mu0
    V0 <- simRes$V0

    A0 <- A+matrix(rnorm(n=length(A), sd=initialCondNoiseSD), ncol=ncol(A))
    # A0 <- A
    # Gamma0 <- Gamma+makeSymmetricMatrix(m=matrix(rnorm(n=length(Gamma), sd=initialCondNoiseSD), ncol=ncol(Gamma)))
    Gamma0 <- Gamma
    C0 <- C+matrix(rnorm(n=length(C), sd=initialCondNoiseSD), ncol=ncol(C))
    # C0 <- C
    # Sigma0 <- Sigma+makeSymmetricMatrix(m=matrix(rnorm(n=length(Sigma), sd=sd(Sigma)/initialCondSNR), ncol=ncol(Sigma)))
    Sigma0 <- Sigma
    mu00 <- mu0+rnorm(n=length(mu0), sd=initialCondNoiseSD)
    # mu00 <- mu0
    V00 <- V0+makeSymmetricMatrix(m=matrix(rnorm(n=length(V0), sd=initialCondNoiseSD), ncol=ncol(V0)))
    # V00 <- V0

    res <- filterLDS(x=x, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    show(sprintf("True data log likelihood: %.4f", sum(log(res$c))))

    res <- emEstimationKF(x=x, A0=A0, Gamma0=Gamma0, C0=C0, Sigma0=Sigma0, mu00=mu00, V00=V00, nIter=nIter)
    plot(unlist(res$logLikelihoods), type="b")
    browser()
}

processAll()
