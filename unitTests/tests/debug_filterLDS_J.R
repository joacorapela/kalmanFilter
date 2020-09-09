
require(mvtnorm)
require(ramcmc)
require(nlme)
source("../../src/filterLDS_R.R")
source("../../src/lsolve.R")
source("../../shumwayAndStoffer11/Kfilter0.R")
source("~/dev/research/programs/src/R/math/l2Norm.R")

processAll <- function() {
    tol <- 1e-3
    stateType0 <- "init00"
    V0SD <- 1e4
    simulationFilename <- "../../scripts/results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    nLatents <- nrow(simRes$A)
    nNeurons <- nrow(simRes$C)
    nObs <- ncol(simRes$x)

    obs <- simRes$x
    x0 <- simRes$mu0
    V0 <- simRes$V0
    A0 <- simRes$A
    Gamma0 <- simRes$Gamma
    C0 <- simRes$C
    Sigma0 <- simRes$Sigma
    B0 <- matrix(0, nrow=nLatents, ncol=1)
    D0 <- matrix(0, nrow=nNeurons, ncol=1)
    x0 <- rep(0, times=nLatents)
    V0 <- V0SD^2*diag(rep(1, times=nLatents))

    fRes_filterLDS_R <- filterLDS_R(B=A0, Z=C0, x0=x0, srV0=chol(x=V0), stateType0=stateType0, Q=Q, R=Gamma0, ys=Sigma0)
    fRes_Kfilter0 <- Kfilter0(num=ncol(obs), y=t(obs), A=C0, mu0=x0, Sigma0=V0, Phi=A0, cQ=srQ, cR=srR)

    for(n in 1:nObs) {
        diff <- fRes_filterLDS_R$xnn[, 1, n]-fRes_Kfilter0$xf[, 1, n]
        error <- l2Norm(diff)
        if(error>tol) {
            stop(sprintf("Error comparing state[%d], error=%f>tol=%f", n, error, tol))
        }
    }
    show("Success")
    browser()
}

processAll()
