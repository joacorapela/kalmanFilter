getNearestPD <- function(S) {
    # Nicholas J. Higham, “Computing a Nearest Symmetric Positive Semidefinite
    # Matrix,” Linear Algebra and its Applications, vol. 103, pp. 103–118,
    # 1988.
    svdRes <- svd(x=S)
    H <- svdRes$v%*%diag(svdRes$d)%*%t(svdRes$v)
    answer <- (S+t(S)+H+t(H))/4
    return(answer)
}

filterLDS_SS <- function(ys, B, mu0, V0, Q, Z, R, checkSymmetry=FALSE, symmetryTol=1e-6, checkPD=FALSE, pdTol=1e-6) {
    # N: number of observations
    # M: dim state space
    # P: dim observations
    M <- nrow(B)
    ys <- as.matrix(ys)
    N <- ncol(ys)
    P <- nrow(ys)
    xnn1 <- array(NA, dim=c(M, 1, N))
    Vnn1 <- array(NA, dim=c(M, M, N))
    xnn <- array(NA, dim=c(M, 1, N))
    Vnn <- array(NA, dim=c(M, M, N))
    innov <- array(NA, dim=c(P, 1, N))
    Sn <- array(NA, dim=c(P, P, N))
    logLike <- 0
    x00 <- as.matrix(x0, nrow=M, ncol=1)
    V00 = as.matrix(V0, nrow=M, ncol=M)
    xnn1[,,1] <- B%*%mu0
    Vnn1[,,1] <- B%*%V0%*%t(B)+Q
    Stmp <- Z%*%Vnn1[,,1]%*%t(Z)+R
    Sn[,,1] <- (Stmp+t(Stmp))/2
    Sinv <- solve(Sn[,,1])
    K <- Vnn1[,,1]%*%t(Z)%*%Sinv
    innov[,,1] <- y[1,]-Z%*%xnn1[,,1]
    xnn[,,1] <- xnn1[,,1]+K%*%innov[,,1]
    Vnn[,,1] <- Vnn1[,,1]-K%*%Z%*%Vnn1[,,1]
    logLike <- -N*P*log(2*pi)-log(det(Sn[,,1]))-t(innov[,,1])%*%Sinv%*%innov[,,1]
    for(k in 2:N) {
        xnn1[,,k] <- B%*%xnn[,,k-1]
        Vnn1[,,k] <- B%*%Vnn[,,k-1]%*%t(B)+Q
        Stmp <- Z%*%Vnn1[,,k]%*%t(Z)+R
        Sn[,,k] <- (Stmp+t(Stmp))/2
        Sinv <- solve(Sn[,,k])
        K <- Vnn1[,,k]%*%t(Z)%*TSinv
        innov[,,k] <- ys[k,]-Z%*%xnn1[,,k]
        xnn[,,k] <- xnn1[,,]+K%*%innov[,,k]
        Vnn[,,k] <- Vnn1[,,k]-K%*%Z%*%Vnn1[,,k]
        logLike <- logLike-log(det(Sn[,,k]))-t(innov[,,k])%*%Sinv%*%innov[,,k]
    }
    logLike <- 0.5*logLike
    answer <- list(xnn1=xnn1, Vnn1=Vnn1, xnn=xnn, Vnn=Vnn, innov=innov, KN=K, Sn=Sn, logLike=logLike)
    return(answer)
}
