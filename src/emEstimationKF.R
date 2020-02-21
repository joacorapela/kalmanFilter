
emEstimationKF <- function(x, A0, Gamma0, C0, Sigma0, mu00, V00, nIter, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE)) {

    A <- A0
    Gamma <- Gamma0
    C <- C0
    Sigma <- Sigma0
    mu0 <- mu00
    V0 <- V00
    N <- ncol(x)

    logLikelihoods <- c()
    for(iter in 1:nIter) {
        # start E-step
        res <- filterLDS(x=x, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
        P <- res$P
        K <- res$K
        mu <- res$mu
        V <- res$V
        c <- res$c

        logLikelihood <- sum(log(c)) # (13.63)
        print(sprintf("loglik(iter=%04d)=%.4f", iter, logLikelihood))
        logLikelihoods <- c(logLikelihoods, logLikelihood)

        res <- smoothLDS(A=A, mu=mu, V=V, P=P)
        muHat <- res$muHat
        VHat <- res$VHat
        VHatLag1 <- res$VHatLag1
        # end E-step

        # computer correlation matrices
        if(varsToEstimate$transitionMatrix || 
           varsToEstimate$transitionCovariance ||
           varsToEstimate$observationMatrix ||
           varsToEstimate$observationCovariance) {
            CHat <- list()
            for(n in 1:(N)) {
                CHat[[n]] <- VHat[[n]] + muHat[,n]%*%t(muHat[,n])
            }
        }
        CHatLag1 <- list()
        if(varsToEstimate$transitionMatrix || 
           varsToEstimate$transitionCovariance) {
            for(n in 1:(N-1)) {
                CHatLag1[[n]] <- VHatLag1[[n]] + muHat[,n+1]%*%t(muHat[,n])
            }
        }
        #

        # start M-step
        if(varsToEstimate$initialStateMean) {
            mu0 <- muHat[,1] # (13.110)
        } else {
            mu0 <- NA
        }
        if(varsToEstimate$initialStateCovariance) {
            V0 <- VHat[[1]] # (13.111)
        } else {
            V0 <- NA
        }

        # compute A
        if(varsToEstimate$transitionMatrix) {
            rhs <- Reduce("+", CHatLag1)
            lhs <- Reduce("+", CHat[1:(length(CHat)-1)])
            A <- t(solve(a=lhs, b=t(rhs))) # (13.113)
        } else {
            A <- NA
        }
        #

        # compute Gamma
        if(varsToEstimate$transitionCovariance) {
            sum <- matrix(data=0, nrow=nrow(Gamma), ncol=ncol(Gamma))
            for(n in 2:N) {
                sum <- sum+CHat[[n]]-A%*%t(CHatLag1[[n-1]])-CHatLag1[[n-1]]%*%t(A)+A%*%CHat[[n-1]]%*%t(A)
            }
            Gamma <- sum/(N-1)
        } else {
            Gamma <- NA
        }
        #

        # compute C
        if(varsToEstimate$observationMatrix) {
            rhs <- x[,1]%*%t(muHat[,1])
            for(n in 2:N) {
                rhs <- rhs + x[,n]%*%t(muHat[,n])
            }
            if(!varsToEstimate$transitionMatrix) {
                lhs <- Reduce("+", CHat[1:(length(CHat)-1)])
            }
            lhs <- lhs+CHat[[length(CHat)]]
            C <- t(solve(a=lhs, b=t(rhs)))
        } else {
            C <- NA
        }
        #

        # compute Sigma
        if(varsToEstimate$observationCovariance) {
            rhs <- x[,1]%*%t(muHat[,1])
            sum <- matrix(data=0, nrow=nrow(Sigma), ncol=ncol(Sigma))
            for(n in 1:N) {
                sum <- sum + x[,n]%*%t(x[,n])-C%*%muHat[,n]%*%t(x[,n])-x[,n]%*%t(muHat[,n])%*%t(C)+C%*%CHat[[n]]%*%t(C)
            }
            Sigma <- sum/N
        } else {
            Sigma <- NA
        }
        #

        # end M-step

# browser()
    }
    answer <- list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0, logLikelihoods=logLikelihoods)
    return(answer)
}
