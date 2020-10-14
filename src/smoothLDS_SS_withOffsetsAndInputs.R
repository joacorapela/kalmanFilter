smoothLDS_SS_withOffsetsAndInputs <- function(B, u, C, c, Q, xnn, Vnn, xnn1, Vnn1, initStateAt=0, m0=NA, V0=NA) {
    if(initStateAt==1 && (!is.na(m0) || !is.na(V0)))
        warning("m0 and V0 are not used when initStateAt==1")
    if(initStateAt==0 && (is.na(m0) || is.na(V0)))
        stop("m0 and V0 are needed when initStateAt==0")
    nObs <- dim(xnn)[3]
    M <- nrow(B)
    xnN <- array(NA, dim=c(M, 1, nObs))
    VnN <- array(NA, dim=c(M, M, nObs))
    Jn <- array(NA, dim=c(M, M, nObs))

    xnN[,,nObs] <- xnn[,,nObs]
    VnN[,,nObs] <- Vnn[,,nObs]
    for(n in nObs:2) {
        Jn[,,n-1] <- t(solve(Vnn1[,,n], B%*%Vnn[,,n-1]))
        xnN[,,n-1] <- xnn[,,n-1]+Jn[,,n-1]%*%(xnN[,,n]-xnn1[,,n])-(Vnn[,,n-1]-Jn[,,n-1]%*%Vnn1[,,n]%*%t(Jn[,,n-1]))%*%t(B)%*%solve(Q,u+C%*%c[,,n])
        VnN[,,n-1] <- Vnn[,,n-1]+Jn[,,n-1]%*%(VnN[,,n]-Vnn1[,,n])%*%t(Jn[,,n-1])
    }
    if(initStateAt==1) {
        # initial state x01 and V01
        # no need to return the smooth estimates of the state at time 0: x0N and V0N
        answer <- list(xnN=xnN, VnN=VnN, Jn=Jn)
        return(answer)
    } else {
        if(initStateAt==0) {
            # initial state m0 and V0
            # return the smooth estimates of the state at time 0: x0N and V0N
            J0 <- t(solve(Vnn1[,,1], B%*%V0))
            x0N <- m0+J0%*%(xnN[,,1]-xnn1[,,1])-(V0-J0%*%Vnn1[,,1]%*%t(J0))%*%t(B)%*%solve(Q, u+C%*%c[,,1])
            V0N <- V0+J0%*%(VnN[,,1]-Vnn1[,,1])%*%t(J0)
            answer <- list(xnN=xnN, VnN=VnN, Jn=Jn, x0N=x0N, V0N=V0N, J0=J0)
            # browser()
            return(answer)
        } else {
            stop(sprintf("Invalid initialStateAt=%d", initStateAt))
        }
    }
}
