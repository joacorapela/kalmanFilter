squareRootKF <- function(A, B, C, D, x0, initialStateType, SRSigmaX0, SRSigmaW, SRSigmaV, us, zs, checkSymmetry=FALSE, symmetryTol=1e-6, checkPD=FALSE, pdTol=1e-6) {
    # from Kalman Filter course by Gregory Plet
    # thttp://mocha-java.uccs.edu/ECE5550/ECE5550-Notes05.pdf
    # Translation between Plett's (this file) and Bishop's notation
    # (filterLDS.py)
    # x=zs
    # x0=mu0
    # SigmaX0%*%t(SigmaX0)=V0
    # SRSigmaW%*%t(SRSigmaW)=Gamma
    # SRSigmaV%*%t(SRSigmaV)=Sigma
    # us=0
    # B=0
    # d=0
    # Return values:
    # SigmaX[[n]] = P[[n-1]]
    # SigmaX[[2:n]] = P[[1:n-1]]
    # L = K
    # x = mu
    # SigmaXHat = V
    leftSolve <- function(a, b) {
        return(solve(t(t(a), t(b))))
    }
    Sigma <- SRSigmaV%*%t(SRSigmaV)
    nObs <- ncol(zs)
    SigmaX <- list()
    L <- list()
    x <- matrix(0, nrow=length(x0), ncol=nObs)
    SigmaXHat <- list()
    c <- rep(NA, times=nObs)
    mean_c1 <- C%*%x0
    sigma_c1 <- C%*%(SRSigmaX0%*%t(SRSigmaX0))%*%t(C)+SRSigmaV%*%t(SRSigmaV)
    # sigma_c1 <- makeSymmetricMatrix(m=sigmaC1)
    c[1] <- dmvnorm(x=zs[,1], mean=mean_c1 , sigma=sigma_c1)
    SRSigmaX <- SRSigmaX0
    for(k in 1:nObs) {
        z <- zs[, k]
        u <- us[, k]
        # SR-KF Step 1a: State estimate time update
        if(k==1) {
            if(initialStateType=="init00") {
                x[,k] <- A%*%x0 + B%*%u # use prior value of "u"
            } else {
                if(initialStateType=="init01") {
                    x[,k] <- x0
                } else {
                    error(sprintf("Invalid initial state type initialStateType=%s.\nUse initialStateType=init00 if the initial state corresponds to E[x_0|y_0] or initialStateType=init10 if it corresponds to E[x_1|y_0]", initialStateType))
                }
            }
        } else {
            x[,k] <- A%*%x[,k-1] + B%*%u # use prior value of "u"
        }
        # SR-KF Step 1b: Error covariance time update
        SRSigmaX <- t(qr(t(cbind(A%*%SRSigmaX, SRSigmaW)))$qr)
        SRSigmaX <- matrix(SRSigmaX[1:nrow(x),1:nrow(x)], nrow=nrow(x)) # using matrix in case dim(SRSigmaX)==(1,1)
        SRSigmaX[upper.tri(x=SRSigmaX)] <- 0
        SigmaX[[k]] <- SRSigmaX%*%t(SRSigmaX)
        if(checkSymmetry) {
            fDiff <- sum((SigmaX[[k]]-t(SigmaX[[k]]))^2)
            if(fDiff>symmetryTol) {
                stop(sprintf("SigmaX in iteration %d is not symmetric.\nThe Frobenius norm of the difference between SigmaX and t(SigmaX) is %f", k, fDiff))
            }
        }
        if(checkPD) {
            minEigenvalue <- min(eigen(x=SigmaX[[k]], only.values=TRUE)$values)
            if(minEigenvalue<pdTol) {
                stop(sprintf("SigmaX in iteration %d may not be positive definite. It minimum eigenvalue is %f", k, minEigenvalue))
            }
        }
        # SR-KF Step 1c: Estimate system output
        zHat <- C%*%x[,k]+D%*%u
        # SR-KF Step 2a: Compute Kalman gain matrix
        # Note: "help mrdivide" to see how "division" is implemented
        SRSigmaZ <- t(qr(t(cbind(C%*%SRSigmaX, SRSigmaV)))$qr)
        SRSigmaZ <- SRSigmaZ[1:nrow(zHat),1:nrow(zHat)]
        SRSigmaZ[upper.tri(x=SRSigmaZ)] <- 0
        if(checkPD) {
            Sn <- SRSigmaZ%*%t(SRSigmaZ)
            minEigenvalue <- min(eigen(x=Sn, only.values=TRUE)$values)
            show(sprintf("Minimum eigenvalue of Sn is %f", minEigenvalue))
            if(minEigenvalue<pdTol) {
                stop(sprintf("Sn in iteration %d may not be positive definite. It minimum eigenvalue is %f", k, minEigenvalue))
            }
        }
        L[[k]] <- lsolve(a=SRSigmaZ, b=lsolve(a=t(SRSigmaZ), b=(SRSigmaX%*%t(SRSigmaX))%*%t(C)))
        # SR-KF Step 2b: State estimate measurement update
        x[,k] <- x[,k]+L[[k]]%*%(z-zHat);
        # SR-KF Step 2c: Error covariance measurement update
        SRSigmaX <- chol_downdate_higherOrder(L=SRSigmaX, U=L[[k]]%*%SRSigmaZ)
        # [Store information for evaluation/plotting purposes]
        SigmaXHat[[k]] <- SRSigmaX%*%t(SRSigmaX)
        if(checkSymmetry) {
            fDiff <- sum((SigmaXHat[[k]]-t(SigmaXHat[[k]]))^2)
            if(fDiff>symmetryTol) {
                stop(sprintf("SigmaXHat in iteration %d is not symmetric.\nThe Frobenius norm of the difference between SigmaXHat and t(SigmaXHat) is %f", k, fDiff))
            }
        }
        if(checkPD) {
            minEigenvalue <- min(eigen(x=SigmaXHat[[k]], only.values=TRUE)$values)
            show(sprintf("Minimum eigenvalue of SigmaXHat is %f", minEigenvalue))
            if(minEigenvalue<pdTol) {
                stop(sprintf("SigmaXHat in iteration %d may not be positive definite. It minimum eigenvalue is %f", k, minEigenvalue))
            }
        }
        if(k>1) {
            mean_cn <- C%*%A%*%x[,k-1]
            sigma_cn <- C%*%SigmaX[[k]]%*%t(C)+Sigma
            # sigmaCn <- makeSymmetricMatrix(m=sigmaCn)
            c[k] <- dmvnorm(x=z, mean=mean_cn, sigma=sigma_cn)
        }
    }
    # answer <- list(P=P, K=K, mu=mu, V=V, c=c)
    answer <- list(SigmaX=SigmaX, L=L, x=x, SigmaXHat=SigmaXHat, c=c)
    return(answer)
}
