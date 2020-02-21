squareRootKF <- function(A, B, C, D, xHat0, SRSigmaX0, SRSigmaW, SRSigmaV, us, zs) {
#  filterLDS <- function(x, A, Gamma, C, Sigma, mu0, V0) {
    # from Kalman Filter course by Gregory Plet
    # thttp://mocha-java.uccs.edu/ECE5550/ECE5550-Notes05.pdf
    # Translation between Plett's (this file) and Bishop's notation
    # (filterLDS.py)
    # x=zs
    # xHat0=mu0
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
    # xHat = mu
    # SigmaXHat = V
    leftSolve <- function(a, b) {
        return(solve(t(t(a), t(b))))
    }
    nObs <- ncol(zs)
    SigmaX <- list()
    L <- list()
    xHat <- matrix(0, nrow=length(xHat0), ncol=nObs)
    SigmaXHat <- list()
    c <- rep(NA, times=nObs)
    mean_c1 <- C%*%xHat0
    sigma_c1 <- C%*%(SRSigmaX0%*%t(SRSigmaX0))%*%t(C)+SRSigmaV%*%t(SRSigmaV)
    # sigma_c1 <- makeSymmetricMatrix(m=sigmaC1)
    c[1] <- dmvnorm(x=zs[,1], mean=mean_c1 , sigma=sigma_c1)
    SRSigmaX <- SRSigmaX0
    for(k in 1:nObs) {
        z <- zs[, k]
        u <- us[, k]
        # SR-KF Step 1a: State estimate time update
        if(k==1) {
            xHat[,k] <- A%*%xHat0 + B%*%u; # use prior value of "u"
        } else {
            xHat[,k] <- A%*%xHat[,k-1] + B%*%u; # use prior value of "u"
        }
        # SR-KF Step 1b: Error covariance time update
        SRSigmaX <- t(qr(t(cbind(A%*%SRSigmaX, SRSigmaW)))$qr)
        SRSigmaX <- SRSigmaX[1:nrow(xHat),1:nrow(xHat)]
        SRSigmaX[upper.tri(x=SRSigmaX)] <- 0
        SigmaX[[k]] <- SRSigmaX%*%t(SRSigmaX)
        # SR-KF Step 1c: Estimate system output
        zhat <- C%*%xHat[,k]+D%*%u
        # SR-KF Step 2a: Compute Kalman gain matrix
        # Note: "help mrdivide" to see how "division" is implemented
        SRSigmaZ <- t(qr(t(cbind(C%*%SRSigmaX, SRSigmaV)))$qr)
        SRSigmaZ <- SRSigmaZ[1:nrow(zhat),1:nrow(zhat)]
        SRSigmaZ[upper.tri(x=SRSigmaZ)] <- 0
        L[[k]] <- lsolve(a=SRSigmaZ, b=lsolve(a=t(SRSigmaZ), b=(SRSigmaX%*%t(SRSigmaX))%*%t(C)))
        # SR-KF Step 2b: State estimate measurement update
        xHat[,k] <- xHat[,k]+L[[k]]%*%(z-zhat);
        # SR-KF Step 2c: Error covariance measurement update
        SRSigmaX <- chol_downdate_higherOrder(L=SRSigmaX, U=L[[k]]%*%SRSigmaZ)
        # [Store information for evaluation/plotting purposes]
        SigmaXHat[[k]] <- SRSigmaX%*%t(SRSigmaX)
        if(k>1) {
            mean_cn <- C%*%A%*%xHat[,k-1]
            sigma_cn <- C%*%SigmaX[[k]]%*%t(C)+SRSigmaV%*%t(SRSigmaV)
            # sigmaCn <- makeSymmetricMatrix(m=sigmaCn)
            c[k] <- dmvnorm(x=z, mean=mean_cn, sigma=sigma_cn)
        }
    }
    # answer <- list(P=P, K=K, mu=mu, V=V, c=c)
    answer <- list(SigmaX=SigmaX, L=L, xHat=xHat, SigmaXHat=SigmaXHat, c=c)
    return(answer)
}
