require(MASS)

simulatePLDSwithOffsetsAndInputs <- function(B, u, C, c, Q, Z, a, D, d, x00, V00, eta=exp) {
    nObs <- ncol(c)
    ssDim <- nrow(B)
    obsDim <- nrow(Z)
    # state noise
    w <- t(mvrnorm(n=nObs, mu=matrix(0, nrow=ssDim, ncol=1), Sigma=Q))
    # initial state
    x0 <- mvrnorm(n=1, mu=x00, Sigma=V00)

    x <- matrix(NaN, nrow=ssDim, ncol=nObs)
    z <- matrix(NaN, nrow=obsDim, ncol=nObs)
    y <- matrix(NaN, nrow=obsDim, ncol=nObs)
    x[,1] <- B%*%x0+u+C%*%c[,1]+w[,1]
    z[,1] <- Z%*%x[,1]+a+D%*%d[,1]
    y[,1] <- rpois(1, eta(z[,1]))
    for(n in 2:nObs) {
        x[,n] <- B%*%x[,n-1]+u+C%*%c[,n]+w[,n]
        z[,n] <- Z%*%x[,n]+a+D%*%d[,n]
        y[,n] <- rpois(obsDim, eta(z[,n]))
    }
    answer <- list(x=x, z=z, y=y, stateType0="init00")
    return(answer)
}
