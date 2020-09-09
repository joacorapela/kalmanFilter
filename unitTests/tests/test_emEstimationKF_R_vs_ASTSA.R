
require(mvtnorm)
require(ramcmc)
require(nlme)
source("../../src/filterLDS_R.R")
source("../../src/smoothLDS_R.R")
source("../../src/emEstimationKF_R.R")
source("../../src/estimateKFInitialCondFA.R")
source("../../src/lsolve.R")
source("../../shumwayAndStoffer11/Kfilter0.R")
source("../../shumwayAndStoffer11/Ksmooth0.R")
source("../../shumwayAndStoffer11/EM0.R")
source("~/dev/research/programs/src/R/math/l2Norm.R")

processAll <- function() {
    noiseSD <- 1e-16
    tol <- 1e-3
    stateType0 <- "init00"
    V0SD <- 1e4
    maxIter <- 200
    tol <- 1e-3
    seed <- 1234567
    nFactors <- 2
    simulationFilename <- "../../scripts/results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    nLatents <- nrow(simRes$A)
    nNeurons <- nrow(simRes$C)
    nObs <- ncol(simRes$x)

    set.seed(seed)

    obs <- simRes$x
    obsForFA <- t(as.matrix(obs))
    initialConds <- estimateKFInitialCondFA(z=obsForFA, nFactors=nFactors)

    A0 <- initialConds$A
    Gamma0 <- 1e-3*diag(rep(1, ncol(A0)))
    C0 <- initialConds$C
    # C0 <- simRes$C
    Sigma0 <- diag(initialConds$sigmaDiag)
    x0 <- simRes$mu0
    V0 <- simRes$V0

    res_emEstimationKF_R <- emEstimationKF_R(ys=obs, B0=A0, Q0=Gamma0, Z0=simRes$C, R0=Sigma0, x00=x0, V00=V0, maxNIter=maxIter, tol=tol, stateType0="init00", varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=FALSE, observationCovariance=TRUE), checkSymmetry=FALSE, checkPD=FALSE)

    srQ0 <- chol(x=Gamma0)
    srR0 <- chol(x=Sigma0)
    res_EM0_ASTRA <- EM0(num=nObs, y=t(obs), A=simRes$C, mu0=x0, Sigma0=V0, Phi=A0, cQ=srQ0, cR=srR0, max.iter=maxIter, tol=tol)

    if(length(res_emEstimationKF_R$logLike)!=length(res_EM0_ASTRA$like)) {
           stop("Log likelihoods lengths differ")
    }
    for(n in 1:length(res_EM0_ASTRA$like)) {
        diff <- res_emEstimationKF_R$logLike[n]+res_EM0_ASTRA$like[n]
        error <- abs(diff)
        if(error>tol) {
            stop(sprintf("Error comparing log likelihoods at iteration=%d, error=%f>tol=%f", n, error, tol))
        }
    }
    show("Success")
    browser()
}

processAll()
