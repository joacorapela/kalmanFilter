
require(MARSS)
require(mvtnorm)
require(ramcmc)
source("../../src/filterLDS_R.R")
source("../../src/smoothLDS_R.R")
source("../../src/emEstimationKF_R.R")
source("../../src/estimateKFInitialCondPPCA.R")
source("../../src/lsolve.R")
source("../../src/chol_downdate_higherOrder.R")

processAll <- function() {
    dObs <- 2*1
    nFactors <- 2
    tol <- 1e-5
    noiseSD <- 1e-3
    stateType0 <- "init00"
    V0Var <- 1e1
    gamma0Var <- 1e1
    sigma0Var <- 1e1
    maxIter <- 100
    seed <- 1234567
    simulationFilenamePattern <- "../../scripts/results/simulationCircleDObs%02d.RData"

    simulationFilename <- sprintf(simulationFilenamePattern, dObs)
    simRes <- get(load(simulationFilename))
    nLatents <- nrow(simRes$A)
    nNeurons <- nrow(simRes$C)
    nObs <- ncol(simRes$x)

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

    # begin MARSS
    stateInputs <- NA
    obsInputs <- NA

    # begin create model
    B1List <- c()
    for(j in 1:nLatents) {
        for(i in 1:nLatents) {
            B1List <- c(B1List, list(sprintf("b%02d%02d", i, j)))
        }
    }
    B1  <- matrix(B1List, nrow=nLatents)
    U1  <- "zero"
    Q1  <- "unconstrained"
    Z1List <- c()
    for(j in 1:nLatents) {
        for(i in 1:nNeurons) {
            Z1List <- c(Z1List, list(sprintf("z%02d%02d", i, j)))
        }
    }
    Z1  <- matrix(Z1List, nrow=nNeurons)
    if(!is.na(stateInputs)) {
        C1  <- "unconstrained"
        c1 <- stateInputs
    } else {
        C1  <- "zero"
        c1 <- "zero"
    }
    A1  <- "zero"
    R1  <- "unconstrained"
    pi1 <- "unequal"
    V01 <- "unconstrained"
    if(!is.na(obsInputs)) {
        D1  <- "unconstrained"
        d1 <- obsInputs
    } else {
        D1  <- "zero"
        d1 <- "zero"
    }

    if(stateType0=="init00") {
        tinitx = 0
    } else {
        if(stateType0=="init10") {
            tinitx = 1
        } else {
            stop(sprintf("Invalid stateType0=%s", stateType0))
        }
    }

    model.list <- list(B=B1, U=U1, C=C1, c=c1, Q=Q1, Z=Z1, A=A1, D=D1, d=d1, R=R1, x0=pi1, V0=V01, tinitx=tinitx)

    B0_marss <- matrix(as.vector(A0), ncol=1)
    Q0_marss  <- matrix(Gamma0[lower.tri(Gamma0, diag=TRUE)], ncol=1)
    Z0_marss  <- matrix(as.vector(C0), ncol=1)
    R0_marss  <- matrix(Sigma0[lower.tri(Sigma0, diag=TRUE)], ncol=1)
    x0_marss  <- matrix(x0, ncol=1)
    V0_marss  <- matrix(V0[lower.tri(V0, diag=TRUE)], ncol=1)
    inits <- list(B=B0_marss, Q=Q0_marss, Z=Z0_marss , R=R0_marss, x0=x0_marss, V0=V0_marss)

    resMARSS <- MARSS(obs, model=model.list, inits=inits, fit=TRUE, silent=2)
    # end MARSS

    res_emEstimationKF_R <- emEstimationKF_R(ys=obs, B0=A0, Q0=Gamma0, Z0=C0, R0=Sigma0, x00=x0, V00=V0, maxNIter=maxIter, tol=tol, stateType0=stateType0, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE), checkSymmetry=FALSE, checkPD=FALSE)

    browser()
}

processAll()
