
require(mvtnorm)
require(ramcmc)
require(nlme)
source("../../src/squareRootKF.R")
source("../../src/lsolve.R")
source("../../src/chol_downdate_higherOrder.R")
source("../../src/filterLDS.R")
source("~/dev/research/programs/src/R/math/l2Norm.R")

processAll <- function() {
    noiseSD <- 1e-4
    tol <- 1e-2
    stateType0 <- "init10"
    V0SD <- 1e4
    simulationFilename <- "../../scripts/results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    nLatents <- nrow(simRes$A)
    nNeurons <- nrow(simRes$C)
    nObs <- ncol(simRes$x)

    obs <- simRes$x
    x0 <- simRes$mu0 + rnorm(n=nLatents, sd=noiseSD)
    V0 <- simRes$V0 + matrix(rnorm(n=nLatents^2, sd=noiseSD), nrow=nLatents)
    A0 <- simRes$A + matrix(rnorm(n=nLatents^2, sd=noiseSD), ncol=nLatents)
    Gamma0 <- simRes$Gamma + matrix(rnorm(n=nLatents^2, sd=noiseSD), ncol=nLatents)
    C0 <- simRes$C + matrix(rnorm(n=nNeurons*nLatents, sd=noiseSD), ncol=nLatents)
    Sigma0 <- simRes$Sigma + matrix(rnorm(n=nNeurons^2, sd=noiseSD), ncol=nNeurons)
    x0 <- rep(0, times=nLatents)
    V0 <- V0SD^2*diag(rep(1, times=nLatents))

    fRes_squareRootKF <- squareRootKF(B=A0, Z=C0, x0=x0, srV0=chol(x=V0), stateType0=stateType0, srQ=chol(x=Gamma0), srR=chol(x=Sigma0), ys=obs)
    fRes_filterLDS <- filterLDS(x=obs, A=A0, Gamma=Gamma0, C=C0, Sigma=Sigma0, mu0=x0, V0=V0)

    for(n in 1:nObs) {
        diff <- fRes_squareRootKF$xnn[, 1, n]-fRes_filterLDS$mu[,n]
        error <- l2Norm(diff)
        if(error>tol) {
            show(sprintf("Error comparing state[%d], error=%f>tol=%f", n, error, tol))
        } else {
            show(sprintf("Success comparing state[%d], error=%f<tol=%f", n, error, tol))
        }
    }
    show("Success")
    browser()
}

processAll()
