filterLDS_SS_withOffsetsAndInputs <- function(y, B, u, C, c, Q, m0, V0, Z, a, D, d, R, initStateAt=0) {
    # N: number of observations
    # M: dim state space
    # P: dim observations
    M <- nrow(B)
    N <- ncol(y)
    P <- nrow(y)
    xnn1 <- array(NA, dim=c(M, 1, N))
    Vnn1 <- array(NA, dim=c(M, M, N))
    xnn <- array(NA, dim=c(M, 1, N))
    Vnn <- array(NA, dim=c(M, M, N))
    innov <- array(NA, dim=c(P, 1, N))
    Sn <- array(NA, dim=c(P, P, N))
    m0 <- as.matrix(m0, nrow=M, ncol=1)
    V0 = as.matrix(V0, nrow=M, ncol=M)
    # logLike <- 0
    logLike <- -N*P*log(2*pi)

    for(k in 1:N) {
        # predicted state mean and covariance
        if(k==1) {
            if(initStateAt==0) {
                xnn1[,,k] <- B%*%m0+u+C%*%c[,,1]
                Vnn1[,,k] <- B%*%V0%*%t(B)+Q
            } else {
                if(initStateAt==1) {
                    xnn1[,,k] <- m0
                    Vnn1[,,k] <- V0
                } else {
                    stop(sprintf("Invalid initStateAt value: %s", initStateAt))
                }
            }
        } else {
            xnn1[,,k] <- B%*%xnn[,,k-1]+u+C%*%c[,,k]
            Vnn1[,,k] <- B%*%Vnn[,,k-1]%*%t(B)+Q
        }
        # innovation
        ynn1 <- Z%*%xnn1[,,k]+a+D%*%d[,,k]
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
