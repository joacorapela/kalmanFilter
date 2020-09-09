
require(mvtnorm)
require(ramcmc)
require(nlme)
source("../../src/filterLDS_R.R")
source("../../src/smoothLDS_R.R")
source("../../src/lsolve.R")
source("../../shumwayAndStoffer11/Kfilter0.R")
source("../../shumwayAndStoffer11/Ksmooth0.R")
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
    x0 <- rep(0, times=nLatents)
    V0 <- V0SD^2*diag(rep(1, times=nLatents))

    srQ <- chol(x=Gamma0)
    srR <- chol(x=Sigma0)
    fRes_filterLDS_R <- filterLDS_R(B=A0, Z=C0, x0=x0, V0=V0, stateType0=stateType0, Q=srQ%*%t(srQ), R=srR%*%t(srR), ys=obs)
    sRes_smoothLDS_R <- smoothLDS_R(B=A0, xnn=fRes_filterLDS_R$xnn, Vnn=fRes_filterLDS_R$Vnn, xnn1=fRes_filterLDS_R$xnn1, Vnn1=fRes_filterLDS_R$Vnn1, x00=x0, V00=V0)

    browser()
}

processAll()
