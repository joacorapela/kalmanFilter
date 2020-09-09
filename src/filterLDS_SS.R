getNearestPD <- function(S) {
    # Nicholas J. Higham, “Computing a Nearest Symmetric Positive Semidefinite
    # Matrix,” Linear Algebra and its Applications, vol. 103, pp. 103–118,
    # 1988.
    svdRes <- svd(x=S)
    H <- svdRes$v%*%diag(svdRes$d)%*%t(svdRes$v)
    answer <- (S+t(S)+H+t(H))/4
    return(answer)
}

filterLDS_SS <- function(B, Z, x0, V0, stateType0, Q, R, ys, checkSymmetry=FALSE, symmetryTol=1e-6, checkPD=FALSE, pdTol=1e-6) {
    M <- nrow(B)
    ys <- as.matrix(ys)
    nObs <- ncol(ys)
    N <- nrow(ys)
    xnn1 <- array(NA, dim=c(M, 1, nObs))
    Vnn1 <- array(NA, dim=c(M, M, nObs))
    xnn <- array(NA, dim=c(M, 1, nObs))
    Vnn <- array(NA, dim=c(M, M, nObs))
    innov <- array(NA, dim=c(N, 1, nObs))
    Sn <- array(NA, dim=c(N, N, nObs))
    logLike <- 0
    x00 <- as.matrix(x0, nrow=M, ncol=1)
    V00 = as.matrix(V0, nrow=M, ncol=M)

    for(k in 1:nObs) {
        # predicted state mean and covariance
        if(k==1) {
            if(stateType0=="init00") {
                xnn1[,,k] <- B%*%x00
                Vnn1[,,k] <- B%*%V00%*%t(B)+Q
            } else {
                if(stateType0=="init10") {
                    xnn1[,,k] <- x00
                    Vnn1[,,k] <- V0
                } else {
                    stop(sprintf("Invalid stateType0 value: %s", stateType0))
                }
            }
        } else {
            xnn1[,,k] <- B%*%xnn[,,k-1]
            Vnn1[,,k] <- B%*%Vnn[,,k-1]%*%t(B)+Q
        }
        # innovation
        innov[,,k] <- ys[,k]-Z%*%xnn1[,,k]
        # innovation covariance
        S <- Z%*%Vnn1[,,k]%*%t(Z)+R
        S <- (t(S)+S)/2
        Sn[,,k] <- S
        Sinv <- solve(S)
        # log likelihood
        detS <- det(S)
        if(detS<0) {
            warning(sprintf("det(S)=%f<0", detS))
            pdS <- getNearestPD(S)
            detPDS <- det(pdS)
            show(sprintf("det(pdS)=%f", detPDS))
            # browser()
            S <- pdS
            detS <- detPDS
        }
        logLike <- logLike - log(detS) - t(innov[,,k])%*%Sinv%*%innov[,,k] #SS16 (6.60)
        # Kalman gain
        K <- Vnn1[,,k]%*%t(Z)%*%Sinv
        # filtered mean
        xnn[,,k] <- xnn1[,,k]+K%*%innov[,,k]
        # filtered covariance
        Vnn[,,k] <- Vnn1[,,k]-K%*%Z%*%Vnn1[,,k]
    }
    logLike <- 0.5*logLike
    answer <- list(xnn1=xnn1, Vnn1=Vnn1, xnn=xnn, Vnn=Vnn, innov=innov, KN=K, Sn=Sn, logLike=logLike)
    return(answer)
}
