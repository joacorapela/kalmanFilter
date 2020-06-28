require(MASS)

simulateLDS <- function(nObs, A, Gamma, C, Sigma, mu0) {
    ssDim <- nrow(A)
    obsDim <- nrow(C)
    # state noise
    w <- t(mvrnorm(n=nObs, mu=matrix(0, nrow=ssDim, ncol=1), Sigma=Gamma))
    # measurement noise
    v <- t(mvrnorm(n=nObs, mu=matrix(0, nrow=obsDim, ncol=1), Sigma=Sigma))
    # initial state noise
    z <- matrix(NaN, nrow=ssDim, ncol=nObs)
    x <- matrix(NaN, nrow=obsDim, ncol=nObs)
    z[,1] <- mu0
    for(n in 2:nObs) {
        z[,n] <- A%*%z[,n-1]+w[,n]
    }
    x <- C%*%z+v
    answer <- list(z=z, x=x)
    return(answer)
}
