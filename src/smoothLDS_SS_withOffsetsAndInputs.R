smoothLDS_SS_withOffsetsAndInputs <- function(B, u, C, c, Q, xnn, Vnn, xnn1, Vnn1, stateType0, x00=NA, V00=NA, c0=NA) {
    if(stateType0=="init01" && (!is.na(x00) || !is.na(V00) || !is.na(c0)))
        warning("x00, V00 and c0 are not used when stateType0==init01")
    if(stateType0=="init00" && (is.na(x00) || is.na(V00)))
        stop("x00, V00 and c0 are needed when stateType0==init00")
    nObs <- dim(xnn)[3]
    M <- nrow(B)
    xnN <- array(NA, dim=c(M, 1, nObs))
    VnN <- array(NA, dim=c(M, M, nObs))
    Jn <- array(NA, dim=c(M, M, nObs))

    xnN[,,nObs] <- xnn[,,nObs]
    VnN[,,nObs] <- Vnn[,,nObs]
    for(n in nObs:2) {
        Jn[,,n-1] <- t(solve(Vnn1[,,n], B%*%Vnn[,,n-1]))
        xnN[,,n-1] <- xnn[,,n-1]+Jn[,,n-1]%*%(xnN[,,n]-xnn1[,,n])-(Vnn[,,n-1]-Jn[,,n-1]%*%Vnn1[,,n]%*%t(Jn[,,n-1]))%*%t(B)%*%solve(Q,u+C%*%c[n-1])
        VnN[,,n-1] <- Vnn[,,n-1]+Jn[,,n-1]%*%(VnN[,,n]-Vnn1[,,n])%*%t(Jn[,,n-1])
    }
    if(stateType0=="init01") {
        # initial state x01 and V01
        # no need to return the smooth estimates of the state at time 0: x0N and V0N
        answer <- list(xnN=xnN, VnN=VnN, Jn=Jn)
        return(answer)
    } else {
        if(stateType0=="init00") {
            # initial state x00 and V00
            # return the smooth estimates of the state at time 0: x0N and V0N
            J0 <- t(solve(Vnn1[,,1], B%*%V00))
            V0N <- V00+J0%*%(VnN[,,1]-Vnn1[,,1])%*%t(J0)
            b = solve(Q, u+C%*%c0)
            e = t(B)%*%b
            f = solve(V0N, e)
            x0N <- x00+J0%*%(xnN[,,1]-xnn1[,,n])-f
            answer <- list(xnN=xnN, VnN=VnN, Jn=Jn, x0N=x0N, V0N=V0N, J0=J0)
            # browser()
            return(answer)
        } else {
            stop(sprintf("Invalid initiaal state type %s", stateType0))
        }
    }
}
