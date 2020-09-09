
require(MASS)
require(mvtnorm)
require(ramcmc)
require(nlme)
source("../../src/filterLDS_R.R")
source("../../src/smoothLDS_R.R")
source("../../src/emEstimationKF_R.R")
source("../../src/estimateKFInitialCondPPCA.R")
source("../../src/lsolve.R")
source("../../shumwayAndStoffer11/Kfilter0.R")
source("../../shumwayAndStoffer11/Ksmooth0.R")
source("../../shumwayAndStoffer11/EM0.R")
source("~/dev/research/programs/src/R/math/l2Norm.R")

processAll <- function() {
    dObs <- 2*1
    nFactors <- 2
    tol <- 1e-10
    stateType0 <- "init00"
    V0Var <- 1e1
    gamma0Var <- 1e1
    sigma0Var <- 1e1
    maxIter <- 4000
    seed <- 1234567
    simulationFilenamePattern <- "../../scripts/results/simulationCircleDObs%02d.RData"

    simulationFilename <- sprintf(simulationFilenamePattern, dObs)
    simRes <- get(load(simulationFilename))
    nLatents <- nrow(simRes$A)
    nNeurons <- nrow(simRes$C)
    nObs <- ncol(simRes$x)

    set.seed(seed)

    obs <- simRes$x
    obsForPPCA <- t(as.matrix(obs))
    initialConds <- estimateKFInitialCondPPCA(z=obsForPPCA, nFactors=nFactors)

    A0 <- initialConds$A
    Gamma0 <- gamma0Var*diag(rep(1, ncol(A0)))
    C0 <- initialConds$C
    # C0 <- simRes$C
    Sigma0 <- sigma0Var*diag(rep(1, times=dObs))
    x0 <- simRes$mu0
    V0 <- V0Var*diag(rep(1, times=nFactors))

    res_emEstimationKF_R <- emEstimationKF_R(ys=obs, B0=A0, Q0=Gamma0, Z0=C0, R0=Sigma0, x00=x0, V00=V0, maxNIter=maxIter, tol=tol, stateType0="init00", varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE), checkSymmetry=FALSE, checkPD=FALSE)
    # res_emEstimationKF_R <- emEstimationKF_R(ys=obs, B0=A0, Q0=Gamma0, Z0=simRes$C, R0=Sigma0, x00=x0, V00=V0, maxNIter=maxIter, tol=tol, stateType0="init00", varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=FALSE, observationCovariance=TRUE), checkSymmetry=FALSE, checkPD=FALSE)

    browser()
}

processAll()
