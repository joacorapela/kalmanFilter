require(ramcmc)
source("../src/chol_downdate_higherOrder.R")

processAll <- function() {
    N <- 10
    R <- 3
    aux <- matrix(rnorm(n=N*N), nrow=N)
    A = aux%*%t(aux) # A positive semidefinite
    L <- chol(x=A)

    for(i in 1:R) {
        u <- 1e-3*rnorm(n=N)
        U <- u%*%t(u)
    }
    trueLRank1Downdate <- chol(x=A-u%*%t(u))
    estLRankd1Downdate <- chol_downdate(L=L, u=u)
    trueLRankRDowndate <- chol(x=A-U)
    estLRankdRDowndate <- chol_downdate_higherOrder(L=L, U=U)

    browser()
}

processAll()
