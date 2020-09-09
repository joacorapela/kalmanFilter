
lag1CovSmootherLDS_SS <- function(Z, KN, B, Vnn, Jn, J0) {
    #SS16, Property 6.3
    M <- nrow(KN)
    N <- ncol(KN)
    nObs <- dim(Vnn)[3]
    Vnn1N <- array(NA, dim = c(M, M, nObs))
    eye <- diag(1, M)
    Vnn1N[, , nObs] <- (eye-KN%*%Z)%*%B%*%Vnn[,,nObs-1]
    for (k in nObs:3) {
        Vnn1N[,,k-1] <- Vnn[,,k-1]%*%t(Jn[,,k-2])+Jn[,,k-1]%*%(Vnn1N[,,k]-B%*%Vnn[,,k-1])%*%t(Jn[,,k-2])
    }
    Vnn1N[,,1] <- Vnn[,,1]%*%t(J0)+Jn[,,1]%*%(Vnn1N[,,2]-B%*%Vnn[,,1])%*%t(J0)
    return(Vnn1N)
}

emEstimationKF_SS <- function(ys, B0, Q0, Z0, R0, x00, V00, maxNIter=50, tol=1e-4, stateType0="init00", varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE), checkSymmetry=FALSE, symmetryTol=1e-6, checkPD=FALSE, pdTol=1e-6) {
    B <- as.matrix(B0)
    Q <- as.matrix(Q0)
    Z <- as.matrix(Z0)
    R <- as.matrix(R0)
    x0 <- x00
    V0 <- as.matrix(V00)

    M <- nrow(B0)
    N <- nrow(ys)
    nObs <- ncol(ys)
    ys <- as.matrix(ys)
    cvg <- 1 + tol
    logLike <- matrix(0, maxNIter, 1)
    for(iter in 1:maxNIter) {
        kf <- filterLDS_R(B=B, Z=Z, x0=x0, V0=V0, stateType0=stateType0, Q=Q, R=R, ys=ys, checkSymmetry=checkSymmetry, symmetryTol=symmetryTol, checkPD=checkPD, pdTol=pdTol)
        show(sprintf("LogLike[%04d]=%f", iter, kf$logLike))
        logLike[iter] <- kf$logLike
#         if (iter > 1) {
#             cvg <- (logLike[iter]-logLike[iter-1])/abs(logLike[iter])
#         }
#         if (cvg < 0) {
#             warning("Likelihood Not Increasing")
#         }
#         if (abs(cvg) < tol) {
#             break
#         }
        if(stateType0=="init00") {
            ks <- smoothLDS_R(B=B, xnn=kf$xnn, Vnn=kf$Vnn, xnn1=kf$xnn1, Vnn1=kf$Vnn1, stateType0=stateType0, x00=x0, V00=V0)
        } else {
            ks <- smoothLDS_R(B=B, xnn=kf$xnn, Vnn=kf$Vnn, xnn1=kf$xnn1, Vnn1=kf$Vnn1, stateType0=stateType0, stateType0=stateType0, x00=NA, V00=NA)
        }
        # We want to first estimaate Z and then R, because R depends on Z
        Vnn1N <- lag1CovSmootherLDS_R(Z=Z, KN=kf$KN, B=B, Vnn=kf$Vnn, Jn=ks$Jn, J0=ks$J0)
        S11 <- ks$xnN[,,1]%*%t(ks$xnN[,,1])+ks$VnN[,,1]
        S10 <- ks$xnN[,,1]%*%t(ks$x0N)+Vnn1N[,,1]
        S00 <- ks$x0N%*%t(ks$x0N)+ks$V0N
        if(varsToEstimate$observationMatrix) {
            Z <- ys[,1]%*%t(ks$xnN[,,1])
        }
        for (i in 2:nObs) {
            S11 <- S11+ks$xnN[,,i]%*%t(ks$xnN[,,i])+ks$VnN[,,i]
            S10 <- S10+ks$xnN[,,i]%*%t(ks$xnN[,,i-1])+Vnn1N[,,i]
            S00 <- S00+ks$xnN[,,i-1]%*%t(ks$xnN[,,i-1])+ks$VnN[,,i-1]
            if(varsToEstimate$observationMatrix) {
                Z <- Z+ys[,i]%*%t(ks$xnN[,,i])
            }
        }
        if(varsToEstimate$transitionMatrix) {
            B <- S10%*%solve(S00)
        }
        if(varsToEstimate$transitionCovariance) {
            Q <- (S11-B%*%t(S10))/nObs
            Q <- (t(Q)+Q)/2
        }
        if(varsToEstimate$observationMatrix) {
            Z <- Z%*%solve(S11)
        }
        # Now that Z is estimated, lets estimate R, if requested
        if(varsToEstimate$observationCovariance) {
            u <- ys[,1]-Z%*% ks$xnN[,,1]
            R <- u%*%t(u)+Z%*%ks$VnN[,,1]%*%t(Z)
            for (i in 2:nObs) {
                u <- ys[,i]-Z%*%ks$xnN[,,i]
                R <- R+u%*%t(u)+Z%*%ks$VnN[,,i]%*%t(Z)
            }
            R <- R/nObs
        }
        if(varsToEstimate$initialStateMean) {
            if(stateType0=="init00") {
                x0 <- ks$x0N
            } else {
                x0 <- ks$xnN[,,1]
            }
        }
        if(varsToEstimate$initialStateCovariance) {
            if(stateType0=="init00") {
                V0 <- ks$V0N
            } else {
                V0 <- ks$VnN[,,1]
            }
        }
    }
    answer <- list(B=B, Q=Q, Z=Z, R=R, x0=x0, V0=V0, logLike=logLike[1:iter], niter=iter, cvg=cvg)
    return(answer)
}
