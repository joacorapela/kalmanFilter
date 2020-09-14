smoothLDS_SS <- function(B, xnn, Vnn, xnn1, Vnn1, m0, V0) {
    N <- dim(xnn)[3]
    M <- nrow(B)
    xnN <- array(NA, dim=c(M, 1, N))
    VnN <- array(NA, dim=c(M, M, N))
    Jn <- array(NA, dim=c(M, M, N))

    xnN[,,N] <- xnn[,,N]
    VnN[,,N] <- Vnn[,,N]
    for(n in N:2) {
        Jn[,,n-1] <- Vnn[,,n-1]%*%t(B)%*%solve(Vnn1[,,n])
        xnN[,,n-1] <- xnn[,,n-1] + Jn[,,n-1]%*%(xnN[,,n]-xnn1[,,n])
        VnN[,,n-1] <- Vnn[,,n-1]+Jn[,,n-1]%*%(VnN[,,n]-Vnn1[,,n])%*%t(Jn[,,n-1])
    }
    # initial state x00 and V00
    # return the smooth estimates of the state at time 0: x0N and V0N
    J0 <- V0%*%t(B)%*%solve(Vnn1[,,1])
    x0N <- m0+J0%*%(xnN[,,1]-xnn1[,,1])
    V0N <- V0+J0%*%(VnN[,,1]-Vnn1[,,1])%*%t(J0)
    answer <- list(xnN=xnN, VnN=VnN, Jn=Jn, x0N=x0N, V0N=V0N, J0=J0)
    # browser()
    return(answer)
}
