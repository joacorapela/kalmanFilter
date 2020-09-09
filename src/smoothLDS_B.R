smoothLDS_B <- function(A, mu, V, P) {
    nObs <- ncol(mu)
    muHat <- matrix(NaN, nrow=nrow(mu), ncol=nObs)
    VHat <- list()
    # VHatLag1[[n]]=cov(z_{n+1}, z_{n})=J_{n}%*%VHat_{n+1}
    VHatLag1 <- list()

    muHat[,nObs] <- mu[,nObs]
    VHat[[nObs]] <- V[[nObs]]
    for(n in (nObs-1):1) {
        Jn <- t(solve(P[[n]], A%*%V[[n]]))
        muHat[,n] <- mu[,n] + Jn%*%(muHat[,n+1]-A%*%mu[,n])
        VHat[[n]] <- V[[n]]+Jn%*%(VHat[[n+1]]-P[[n]])%*%t(Jn)
        VHatLag1[[n]] <- Jn%*%VHat[[n+1]]
    }
    answer <- list(muHat=muHat, VHat=VHat, VHatLag1=VHatLag1)
    return(answer)
}
