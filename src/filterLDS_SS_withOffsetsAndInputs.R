filterLDS_SS_withOffsetsAndInputs <- function(y, B, u, C, c, Q, x00, V00, stateType0, Z, a, D, d, R) {
    M <- nrow(B)
    nObs <- ncol(y)
    N <- nrow(y)
    xnn1 <- array(NA, dim=c(M, 1, nObs))
    Vnn1 <- array(NA, dim=c(M, M, nObs))
    xnn <- array(NA, dim=c(M, 1, nObs))
    Vnn <- array(NA, dim=c(M, M, nObs))
    innov <- array(NA, dim=c(N, 1, nObs))
    Sn <- array(NA, dim=c(N, N, nObs))
    x00 <- as.matrix(x00, nrow=M, ncol=1)
    V00 = as.matrix(V00, nrow=M, ncol=M)
    logLike <- 0

    for(k in 1:nObs) {
        # predicted state mean and covariance
        if(k==1) {
            if(stateType0=="init00") {
                xnn1[,,k] <- B%*%x00+u+C%*%c[1]
                Vnn1[,,k] <- B%*%V00%*%t(B)+Q
            } else {
                if(stateType0=="init10") {
                    xnn1[,,k] <- x00
                    Vnn1[,,k] <- V00
                } else {
                    stop(sprintf("Invalid stateType0 value: %s", stateType0))
                }
            }
        } else {
            xnn1[,,k] <- B%*%xnn[,,k-1]+u+C%*%c[k]
            Vnn1[,,k] <- B%*%Vnn[,,k-1]%*%t(B)+Q
        }
        # innovation
        ynn1 <- Z%*%xnn1[,,k]+a+D%*%d[k]
        innov[,,k] <- y[,k]-ynn1
        # innovation covariance
        S <- Z%*%Vnn1[,,k]%*%t(Z)+R
        S <- (t(S)+S)/2
        Sn[,,k] <- S
        Sinv <- solve(S)
        # log likelihood
        detS <- det(S)
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
