require(ramcmc)
source("../../src/chol_downdate_higherOrder.R")

processAll <- function() {
    tol <- 1e-4
    uScale <- 1e-2
    N <- 10
    R <- 3
    aux <- matrix(rnorm(n=N*N), nrow=N)
    A = aux%*%t(aux) # A positive semidefinite
    L <- chol(x=A)

    u <- matrix(uScale*rnorm(n=N), nrow=N)
    U <- matrix(uScale*rnorm(n=N*N), nrow=N)
    trueLRank1Downdate <- chol(x=A-u%*%t(u))
    estLRankd1Downdate <- chol_downdate(L=L, u=u)
    trueLRankRDowndate <- chol(x=A-U%*%t(U))
    estLRankdRDowndate <- chol_downdate_higherOrder(L=L, U=U)

    errorRank1 <- mean((trueLRank1Downdate-estLRankd1Downdate)^2)
    if(errorRank1>tol) {
        stop(sprintf("Failed chol_downdate test with error=%f>tol=%f", errorRank1, tol))
    } else {
        show(sprintf("Passed chold_downdate test with error=%f<tol=%f", errorRank1, tol))
    }
    errorRankR <- mean((trueLRankRDowndate-estLRankdRDowndate)^2)
    if(errorRankR>tol) {
        stop(sprintf("Failed chol_downdate_highOrder test with error=%f>tol=%f", errorRankR, tol))
    } else {
        show(sprintf("Passed chold_downdate_highOrder test with error=%f<tol=%f", errorRankR, tol))
    }
    browser()
}

processAll()
