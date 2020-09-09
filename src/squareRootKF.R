squareRootKF <- function(B, Z, x0, srV0, stateType0, srQ, srR, ys, checkSymmetry=FALSE, symmetryTol=1e-6, checkPD=FALSE, pdTol=1e-6) {
    # from Kalman Filter course by Gregory Plet
    # thttp://mocha-java.uccs.edu/ECE5550/ECE5550-Notes05.pdf
    # x0 \in R^M
    # ys \in R^{N\times nObs}

    srGetUpdatedCov <- function(modelMat, srStateCov, srNoiseCov) {
        srUpdatedCov <- t(qr(t(cbind(modelMat%*%srStateCov, srNoiseCov)))$qr)
        srUpdatedCov <- srUpdatedCov[1:nrow(modelMat), 1:nrow(modelMat)]
        srUpdatedCov[upper.tri(x=srUpdatedCov)] <- 0
        return(srUpdatedCov)
    }

    M <- nrow(B)
    ys <- as.matrix(ys)
    nObs <- ncol(ys)
    N <- nrow(ys)
    xnn1 <- array(NA, dim=c(M, 1, nObs))
    Vnn1 <- array(NA, dim=c(M, M, nObs))
    xnn <- array(NA, dim=c(M, 1, nObs))
    Vnn <- array(NA, dim=c(M, M, nObs))
    innov <- array(NA, dim=c(N, 1, nObs))
    Kn <- array(NA, dim=c(M, N, nObs))
    Sn <- array(NA, dim=c(N, N, nObs))
    logLike <- 0
    x00 <- as.matrix(x0, nrow=M, ncol=1)
    srV00 = as.matrix(srV0, nrow=M, ncol=M)

    for(k in 1:nObs) {
        # predicted state mean and covariance
        if(k==1) {
            if(stateType0=="init00") {
                xnn1[,,k] <- B%*%x00
                srVnn1 <- srGetUpdatedCov(modelMat=B, srStateCov=srV0, srNoiseCov=srQ)
            } else {
                if(stateType0=="init10") {
                    xnn1[,,k] <- x00
                    srVnn1 <- srV0
                } else {
                    stop(sprintf("Invalid stateType0 value: %s", stateType0))
                }
            }
        } else {
            xnn1[,,k] <- B%*%xnn[,,k-1]
            srVnn1 <- srGetUpdatedCov(modelMat=B, srStateCov=srVnn, srNoiseCov=srQ)
        }
        Vnn1[,,k] <- srVnn1%*%t(srVnn1)
        # innovation
        innov[,,k] <- ys[,k]-Z%*%xnn1[,,k]
        # innovation covariance
        srS <- srGetUpdatedCov(modelMat=Z, srStateCov=srVnn1, srNoiseCov=srR)
        S <- srS%*%t(srS)
        Sn[,,k] <- S
        # log likelihood
        logLike <- logLike - logDet(pdSymm(S)) - t(innov[,,k])%*%solve(S, innov[,,k])
        # Kalman gain
        K <- lsolve(a=srS, b=lsolve(a=t(srS), b=Vnn1[,,k]%*%t(Z)))
        # filtered mean
        xnn[,,k] <- xnn1[,,k] + K%*%innov[,,k]
        # filtered covariance
        srVnn <- chol_downdate_higherOrder(L=srVnn1, U=K%*%srS)
        Vnn[,,k] <- srVnn%*%t(srVnn)
    }
    answer <- list(xnn1=xnn1, Vnn1=Vnn1, xnn=xnn, Vnn=Vnn, innov=innov, KN=K, Sn=Sn, logLike=logLike)
    return(answer)
}
