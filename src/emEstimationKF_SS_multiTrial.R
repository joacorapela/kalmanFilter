
lag1CovSmootherLDS_SS <- function(Z, KN, B, Vnn, Jn, J0) {
    #SS16, Property 6.3
    M <- nrow(KN)
    N <- dim(Vnn)[3]
    Vnn1N <- array(NA, dim = c(M, M, N))
    eye <- diag(1, M)
    Vnn1N[, , N] <- (eye-KN%*%Z)%*%B%*%Vnn[,,N-1]
    for (k in N:3) {
        Vnn1N[,,k-1] <- Vnn[,,k-1]%*%t(Jn[,,k-2])+Jn[,,k-1]%*%(Vnn1N[,,k]-B%*%Vnn[,,k-1])%*%t(Jn[,,k-2])
    }
    Vnn1N[,,1] <- Vnn[,,1]%*%t(J0)+Jn[,,1]%*%(Vnn1N[,,2]-B%*%Vnn[,,1])%*%t(J0)
    return(Vnn1N)
}

emEstimationKF_SS_multiTrial <- function(y, B0, Q0, Z0, R0, m0, V0, maxIter=50, tol=1e-4, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE)) {
    B <- B0
    Q <- Q0
    Z <- Z0
    R <- R0
    V0 <- V0

    M <- nrow(B0)
    N <- nrow(Z0)
    nTrials <- dim(y)[1]
    cvg <- 1 + tol
    logLike <- matrix(0, maxIter, 1)
    for(iter in 1:maxIter) {
        logLike[iter] <- 0
        S11 <- matrix(0, nrow=M, ncol=M)
        S10 <- matrix(0, nrow=M, ncol=M)
        S00 <- matrix(0, nrow=N, ncol=M)
        Z <- matrix(0, nrow=M, ncol=M)
        for(r in 1:nTrials) {
            kf <- filterLDS_SS(y=y[r,,], B=B, m0=m0[r,], V0=V0[r,,], Q=Q, Z=Z, R=R)
            logLike[iter] <- logLike[iter]+kf$logLike
#         if (iter > 1) {
#             cvg <- (logLike[iter]-logLike[iter-1])/abs(logLike[iter])
#         }
#         if (cvg < 0) {
#             warning("Likelihood Not Increasing")
#         }
#         if (abs(cvg) < tol) {
#             break
#         }
            ks <- smoothLDS_SS(B=B, xnn=kf$xnn, Vnn=kf$Vnn, xnn1=kf$xnn1, Vnn1=kf$Vnn1, m0=m0[r,], V0=V0[r,,])
            # We want to first estimaate Z and then R, because R depends on Z
            Vnn1N <- lag1CovSmootherLDS_SS(Z=Z, KN=kf$KN, B=B, Vnn=kf$Vnn, Jn=ks$Jn, J0=ks$J0)
            S11 <- S11+ks$xnN[,,1]%*%t(ks$xnN[,,1])+ks$VnN[,,1]
            S10 <- S10+ks$xnN[,,1]%*%t(ks$x0N)+Vnn1N[,,1]
            S00 <- S00+ks$x0N%*%t(ks$x0N)+ks$V0N
            if(varsToEstimate$observationMatrix) {
                Z <- Z+y[r,,1]%*%t(ks$xnN[,,1])
            }
            for (i in 2:N) {
                S11 <- S11+ks$xnN[,,i]%*%t(ks$xnN[,,i])+ks$VnN[,,i]
                S10 <- S10+ks$xnN[,,i]%*%t(ks$xnN[,,i-1])+Vnn1N[,,i]
                S00 <- S00+ks$xnN[,,i-1]%*%t(ks$xnN[,,i-1])+ks$VnN[,,i-1]
                if(varsToEstimate$observationMatrix) {
                    Z <- Z+y[r,,i]%*%t(ks$xnN[,,i])
                }
            }
            if(varsToEstimate$initialStateMean) {
                m0[r,] <- ks$x0N
            }
            if(varsToEstimate$initialStateCovariance) {
                V0[r,,] <- ks$V0N
            }
        }
        show(sprintf("LogLike[%04d]=%f", iter, logLike[iter]))
        if(varsToEstimate$transitionMatrix) {
            B <- S10%*%solve(S00)
        }
        if(varsToEstimate$transitionCovariance) {
            Q <- (S11-B%*%t(S10))/(N*nTrials)
            Q <- (t(Q)+Q)/2
        }
        if(varsToEstimate$observationMatrix) {
            Z <- Z%*%solve(S11)
        }
        # Now that Z is estimated, lets estimate R, if requested
        if(varsToEstimate$observationCovariance) {
            R <- matrix(0, nrow=N, ncol=N)
            for(r in 1:nTrials) {
                for(n in 1:N) {
                    u <- y[r,,n]-Z%*%ks$xnN[,,n]
                    R <- R+u%*%t(u)+Z%*%ks$VnN[,,n]%*%t(Z)
                }
            }
            R <- R/(N*nTrials)
        }
    }
    answer <- list(B=B, Q=Q, Z=Z, R=R, m0=m0, V0=V0, logLike=logLike[1:iter], niter=iter, cvg=cvg)
    return(answer)
}
