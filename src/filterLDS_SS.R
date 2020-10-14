
filterLDS_SS <- function(y, B, m0, V0, Q, Z, R) {
    # N: number of observations
    # M: dim state space
    # P: dim observations
    M <- nrow(B)
    # y <- as.matrix(y)
    N <- ncol(y)
    P <- nrow(y)
    xnn1 <- array(NA, dim=c(M, 1, N))
    Vnn1 <- array(NA, dim=c(M, M, N))
    xnn <- array(NA, dim=c(M, 1, N))
    Vnn <- array(NA, dim=c(M, M, N))
    innov <- array(NA, dim=c(P, 1, N))
    Sn <- array(NA, dim=c(P, P, N))
    m0 <- as.matrix(m0, nrow=M, ncol=1)
    V0 <- as.matrix(V0, nrow=M, ncol=M)

    # k==1
    xnn1[,,1] <- B%*%m0
    Vnn1[,,1] <- B%*%V0%*%t(B)+Q
    Stmp <- Z%*%Vnn1[,,1]%*%t(Z)+R
    Sn[,,1] <- (Stmp+t(Stmp))/2
    Sinv <- solve(Sn[,,1])
    K <- Vnn1[,,1]%*%t(Z)%*%Sinv
    innov[,,1] <- y[,1]-Z%*%xnn1[,,1]
    xnn[,,1] <- xnn1[,,1]+K%*%innov[,,1]
    Vnn[,,1] <- Vnn1[,,1]-K%*%Z%*%Vnn1[,,1]
    logLike <- -N*P*log(2*pi)-log(det(Sn[,,1]))-t(innov[,,1])%*%Sinv%*%innov[,,1]

    # k>1
    for(k in 2:N) {
        xnn1[,,k] <- B%*%xnn[,,k-1]
        Vnn1[,,k] <- B%*%Vnn[,,k-1]%*%t(B)+Q
        Stmp <- Z%*%Vnn1[,,k]%*%t(Z)+R
        Sn[,,k] <- (Stmp+t(Stmp))/2
        Sinv <- solve(Sn[,,k])
        K <- Vnn1[,,k]%*%t(Z)%*%Sinv
        innov[,,k] <- y[,k]-Z%*%xnn1[,,k]
        xnn[,,k] <- xnn1[,,k]+K%*%innov[,,k]
        Vnn[,,k] <- Vnn1[,,k]-K%*%Z%*%Vnn1[,,k]
        logLike <- logLike-log(det(Sn[,,k]))-t(innov[,,k])%*%Sinv%*%innov[,,k]
    }
    logLike <- 0.5*logLike
    answer <- list(xnn1=xnn1, Vnn1=Vnn1, xnn=xnn, Vnn=Vnn, innov=innov, KN=K, Sn=Sn, logLike=logLike)
    return(answer)
}
