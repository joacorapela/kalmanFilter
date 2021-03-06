
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
    noiseSD <- 1e-8
    tol <- 1e-3
    stateType0 <- "init00"
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

    srQ <- chol(x=Gamma0)
    srR <- chol(x=Sigma0)
    fRes_filterLDS_R <- filterLDS_R(B=A0, Z=C0, x0=x0, V0=V0, stateType0=stateType0, Q=srQ%*%t(srQ), R=srR%*%t(srR), ys=obs)
    sRes_smoothLDS_R <- smoothLDS_R(B=A0, xnn=fRes_filterLDS_R$xnn, Vnn=fRes_filterLDS_R$Vnn, xnn1=fRes_filterLDS_R$xnn1, Vnn1=fRes_filterLDS_R$Vnn1, x00=x0, V00=V0)
    sRes_Ksmooth0 <- Ksmooth0(num=ncol(obs), y=t(obs), A=C0, mu0=x0, Sigma0=V0, Phi=A0, cQ=srQ, cR=srR)

    for(n in nObs:1) {
        diff <- sRes_smoothLDS_R$xnN[,1,n]-sRes_Ksmooth0$xs[,1,n]
        error <- l2Norm(diff)
        if(error>tol) {
            stop(sprintf("Error comparing state[%d], error=%f>tol=%f", n, error, tol))
        }
    }
    diff <- sRes_smoothLDS_R$x0N-sRes_Ksmooth0$x0n
    error <- l2Norm(diff)
    if(error>tol) {
        stop(sprintf("Error comparing state0, error=%f>tol=%f", n, error, tol))
    }
    show("Success")
    browser()
}

processAll()
