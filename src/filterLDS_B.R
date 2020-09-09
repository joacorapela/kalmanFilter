filterLDS_B <- function(x, A, Gamma, C, Sigma, mu0, V0) {
    # From Bishop 2006, Section 13.3.1
    nObs <- ncol(x)
    P <- list()
    K <- list()
    mu <- matrix(NaN, nrow=length(mu0), ncol=nObs)
    V <- list()
    c <- rep(NA, times=nObs)

    I <- diag(nrow=nrow(A))
    lhsKn <- C%*%V0%*%t(C)+Sigma
    rhsKn <- C%*%V0
    K[[1]] <- t(solve(a=lhsKn, b=rhsKn))
    mu[,1] <- mu0 + K[[1]]%*%(x[,1]-C%*%mu0)
    V[[1]] <- (I-K[[1]]%*%C)%*%V0
    mean_c1 <- C%*%mu0
    sigma_c1 <- C%*%V0%*%t(C)+Sigma
    # sigma_c1 <- makeSymmetricMatrix(m=sigmaC1)
    c[1] <- dmvnorm(x=x[,1], mean=mean_c1 , sigma=sigma_c1)
    for(n in 2:nObs) {
        P[[n-1]] <- A%*%V[[n-1]]%*%t(A)+Gamma

        lhsKn <- C%*%P[[n-1]]%*%t(C)+Sigma
        rhsKn <- C%*%P[[n-1]]
        K[[n]] <- t(solve(a=lhsKn, b=rhsKn))

        mu[,n] <- A%*%mu[,n-1] + K[[n]]%*%(x[,n]-C%*%A%*%mu[,n-1])
        V[[n]] <- P[[n-1]]-K[[n]]%*%C%*%P[[n-1]]

        mean_cn <- C%*%A%*%mu[,n-1]
        sigma_cn <- C%*%P[[n-1]]%*%t(C)+Sigma
        sigma_cn <- (sigma_cn+t(sigma_cn))/2
        c[n] <- dmvnorm(x=x[,n], mean=mean_cn, sigma=sigma_cn)
    }
    answer <- list(P=P, K=K, mu=mu, V=V, c=c)
    return(answer)
}
