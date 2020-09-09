require(MASS)

simulateLDS <- function(nObs, B, Q, Z, R, x00, V00) {
    ssDim <- nrow(B)
    obsDim <- nrow(Z)
    # state noise
    w <- t(mvrnorm(n=nObs, mu=matrix(0, nrow=ssDim, ncol=1), Sigma=Q))
    # measurement noise
    v <- t(mvrnorm(n=nObs, mu=matrix(0, nrow=obsDim, ncol=1), Sigma=R))
    # initial state noise
    x <- matrix(NaN, nrow=ssDim, ncol=nObs)
    y <- matrix(NaN, nrow=obsDim, ncol=nObs)
    x[,1] <- t(mvrnorm(n=1, mu=x00, Sigma=V00))
    for(n in 2:nObs) {
        x[,n] <- B%*%x[,n-1]+w[,n]
    }
    y <- Z%*%x+v
    answer <- list(x=x, y=y)
    return(answer)
}
