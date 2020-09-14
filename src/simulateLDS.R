require(MASS)

simulateLDS <- function(N, B, Q, Z, R, m0, V0) {
    M <- nrow(B)
    P <- nrow(Z)
    # state noise
    w <- t(mvrnorm(n=N, mu=matrix(0, nrow=M, ncol=1), Sigma=Q))
    # measurement noise
    v <- t(mvrnorm(n=N, mu=matrix(0, nrow=P, ncol=1), Sigma=R))
    # initial state noise
    x <- matrix(NaN, nrow=M, ncol=N)
    y <- matrix(NaN, nrow=P, ncol=N)
    x0 <- matrix(mvrnorm(n=1, mu=m0, Sigma=V0), nrow=M)
    x[,1] <- B%*%x0+w[,1]
    for(n in 2:N) {
        x[,n] <- B%*%x[,n-1]+w[,n]
    }
    y <- Z%*%x+v
    answer <- list(x0=x0, x=x, y=y)
    return(answer)
}
