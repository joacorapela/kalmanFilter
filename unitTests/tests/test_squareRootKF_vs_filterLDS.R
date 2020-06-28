
require(mvtnorm)
require(ramcmc)
source("../../src/squareRootKF.R")
source("../../src/lsolve.R")
source("../../src/chol_downdate_higherOrder.R")
source("../../src/filterLDS.R")
source("~/dev/research/programs/src/R/math/l2Norm.R")

processAll <- function() {
    noiseSD <- 1e-4
    tol <- 1e-6
    initialStateType <- "init01"
    V0SD <- 1e4
    simulationFilename <- "../../scripts/results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    nLatents <- nrow(simRes$A)
    nNeurons <- nrow(simRes$C)
    nObs <- ncol(simRes$x)

    obs <- simRes$x
    us <- matrix(0, nrow=1, ncol=nObs)
    x0 <- simRes$mu0 + rnorm(n=nLatents, sd=noiseSD)
    V0 <- simRes$V0 + matrix(rnorm(n=nLatents^2, sd=noiseSD), nrow=nLatents)
    A0 <- simRes$A + matrix(rnorm(n=nLatents^2, sd=noiseSD), ncol=nLatents)
    Gamma0 <- simRes$Gamma + matrix(rnorm(n=nLatents^2, sd=noiseSD), ncol=nLatents)
    C0 <- simRes$C + matrix(rnorm(n=nNeurons*nLatents, sd=noiseSD), ncol=nLatents)
    Sigma0 <- simRes$Sigma + matrix(rnorm(n=nNeurons^2, sd=noiseSD), ncol=nNeurons)
    B0 <- matrix(0, nrow=nLatents, ncol=1)
    D0 <- matrix(0, nrow=nNeurons, ncol=1)
    x0 <- rep(0, times=nLatents)
    V0 <- V0SD^2*diag(rep(1, times=nLatents))

    fRes_squaareRootKF <- squareRootKF(A=A0, B=B0, C=C0, D=D0, x0=x0, initialStateType=initialStateType, SRSigmaX0=chol(x=V0), SRSigmaW=chol(x=Gamma0), SRSigmaV=chol(x=Sigma0), us=us, zs=obs)
    fRes_filterLDS <- filterLDS(x=obs, A=A0, Gamma=Gamma0, C=C0, Sigma=Sigma0, mu0=x0, V0=V0)

    for(n in 1:nObs) {
        diff <- fRes_squaareRootKF$x[,n]-fRes_filterLDS$z[,n]
        error <- l2Norm(diff)
        if(error>tol) {
            stop(sprintf("Error comparing states>%f", tol))
        }
    }
    show("Success")
    browser()
}

processAll()
