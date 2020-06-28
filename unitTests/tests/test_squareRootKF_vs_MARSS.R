
require(MARSS)
require(mvtnorm)
require(ramcmc)
source("../../src/squareRootKF.R")
source("../../src/lsolve.R")
source("../../src/chol_downdate_higherOrder.R")
source("../../src/filterLDS.R")
source("~/dev/research/programs/src/R/math/l2Norm.R")

processAll <- function() {
    tol <- 1e-5
    noiseSD <- 1e-4
    initialStateType <- "init01"
    tinitx <- 1
    V0SD <- 1e2
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
    Gamma0 <- simRes$Gamma + diag(rnorm(n=nLatents, sd=noiseSD))
    C0 <- simRes$C + matrix(rnorm(n=nNeurons*nLatents, sd=noiseSD), ncol=nLatents)
    Sigma0 <- simRes$Sigma + diag(rnorm(n=nNeurons, sd=noiseSD))
    B0 <- matrix(0, nrow=nLatents, ncol=1)
    D0 <- matrix(0, nrow=nNeurons, ncol=1)
    x0 <- rep(0, times=nLatents)
    V0 <- V0SD^2*diag(rep(1, times=nLatents))

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

    model.list <- list(B=B1, U=U1, C=C1, c=c1, Q=Q1, Z=Z1, A=A1, D=D1, d=d1, R=R1, x0=pi1, V0=V01, tinitx=tinitx)

    B0_marss <- matrix(as.vector(A0), ncol=1)
    Q0_marss  <- matrix(Gamma0[lower.tri(Gamma0, diag=TRUE)], ncol=1)
    Z0_marss  <- matrix(as.vector(C0), ncol=1)
    R0_marss  <- matrix(Sigma0[lower.tri(Sigma0, diag=TRUE)], ncol=1)
    x0_marss  <- matrix(x0, ncol=1)
    V0_marss  <- matrix(V0[lower.tri(V0, diag=TRUE)], ncol=1)
    inits <- list(B=B0_marss, Q=Q0_marss, Z=Z0_marss , R=R0_marss, x0=x0_marss, V0=V0_marss)

    MLEobj <- MARSS(obs, model=model.list, inits=inits, fit=FALSE, silent=2)
    MLEobj$par <- MLEobj$start
    kfList <- MARSSkfss(MLEobj)
    fRes_MARSS <- kfList$xtt
    # end MARSS

    fRes_squareRootKF <- squareRootKF(A=A0, B=B0, C=C0, D=D0, x0=x0, initialStateType=initialStateType, SRSigmaX0=chol(x=V0), SRSigmaW=chol(x=Gamma0), SRSigmaV=chol(x=Sigma0), us=us, zs=obs)

    for(n in 1:nObs) {
        diff <- fRes_squareRootKF$x[,n]-fRes_MARSS[,n]
        error <- l2Norm(diff)
        if(error>tol) {
            stop(sprintf("Error comparing state[%d], error=%f>tol=%f", n, error, tol))
        } else {
            show(sprintf("Success comparing state[%d], error=%f<tol=%f", n, error, tol))
        }
    }
    browser()
}

processAll()
