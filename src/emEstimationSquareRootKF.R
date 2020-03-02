
emEstimationSquareRootKF <- function(zs, A0, SRSigmaW0, C0, SRSigmaV0, B0, D0, xHat00, SRSigmaX00, nIter=Inf, tol=-Inf, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE), logMsgPattern="loglik(iter=%%04d)=%%.4f") {

    if(!is.finite(nIter) && !is.finite(tol)) {
        stop("nIter and tol cannot be both infinite")
    }
    A <- A0
    SRSigmaW <- SRSigmaW0
    Gamma <- SRSigmaW%*%t(SRSigmaW)
    C <- C0
    SRSigmaV <- SRSigmaV0
    Sigma <- SRSigmaV%*%t(SRSigmaV)
    B <- B0
    D <- D0
    xHat0 <- xHat00
    SRSigmaX0 <- SRSigmaX00
    V0 <- SRSigmaX0%*%t(SRSigmaX0)
    us <- matrix(0, nrow=1, ncol=ncol(zs))
    N <- ncol(zs)

    logLikelihoods <- c()
    logLikelihood <- -Inf
    logLikelihoodIncrease <- Inf
    iter <- 0
    while(iter<=nIter && tol<abs(logLikelihoodIncrease)) {
        # start E-step
        res <- squareRootKF(A=A, B=B, C=C, D=D, xHat0=xHat0, SRSigmaX0=SRSigmaX0, SRSigmaW=SRSigmaW, SRSigmaV=SRSigmaV, us=us, zs=zs)
        P <- res$SigmaX[2:length(res$SigmaX)]
        mu <- res$xHat
        V <- res$SigmaXHat
        c <- res$c

        lastLogLikelihood <- logLikelihood
        logLikelihood <- sum(log(c)) # (13.63)
        if(is.finite(logLikelihood)) {
            logLikelihoodIncrease <- logLikelihood-lastLogLikelihood
        }
        print(sprintf(logMsgPattern, iter, logLikelihood))
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
            xHat0 <- muHat[,1] # (13.110)
        }
        if(varsToEstimate$initialStateCovariance) {
            V0 <- VHat[[1]] # (13.111)
            SRSigmaX0 <- chol(x=V0)
        }

        # compute A
        if(varsToEstimate$transitionMatrix) {
            rhs <- Reduce("+", CHatLag1)
            lhs <- Reduce("+", CHat[1:(length(CHat)-1)])
            A <- t(solve(a=lhs, b=t(rhs))) # (13.113)
        }
        #

        # compute Gamma
        if(varsToEstimate$transitionCovariance) {
            sum <- matrix(data=0, nrow=nrow(SRSigmaW0), ncol=ncol(SRSigmaW0))
            for(n in 2:N) {
                sum <- sum+CHat[[n]]-A%*%t(CHatLag1[[n-1]])-CHatLag1[[n-1]]%*%t(A)+A%*%CHat[[n-1]]%*%t(A)
            }
            Gamma <- sum/(N-1)
            SRSigmaW <- chol(x=Gamma)
        }
        #

        # compute C
        if(varsToEstimate$observationMatrix) {
            rhs <- zs[,1]%*%t(muHat[,1])
            for(n in 2:N) {
                rhs <- rhs + zs[,n]%*%t(muHat[,n])
            }
            if(!varsToEstimate$transitionMatrix) {
                lhs <- Reduce("+", CHat[1:(length(CHat)-1)])
            }
            lhs <- lhs+CHat[[length(CHat)]]
            C <- t(solve(a=lhs, b=t(rhs)))
        }
        #

        # compute Sigma
        if(varsToEstimate$observationCovariance) {
            rhs <- zs[,1]%*%t(muHat[,1])
            sum <- matrix(data=0, nrow=nrow(SRSigmaV0), ncol=ncol(SRSigmaV0))
            for(n in 1:N) {
                sum <- sum + zs[,n]%*%t(zs[,n])-C%*%muHat[,n]%*%t(zs[,n])-zs[,n]%*%t(muHat[,n])%*%t(C)+C%*%CHat[[n]]%*%t(C)
            }
            Sigma <- sum/N
            SRSigmaV <- chol(x=Sigma)
        }
        #

        # end M-step
        iter <- iter + 1
    }
    show(sprintf("Finished in %d iter with a log likelihood increase of %f", iter, logLikelihoodIncrease))
    answer <- list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, xHat0=xHat0, V0=V0, logLikelihoods=logLikelihoods)
    return(answer)
}
