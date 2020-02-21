chol_downdate_higherOrder <- function(L, U) {
    for(j in 1:ncol(U)) {
        L <- chol_downdate(L, U[,j])
    }
    return(L)
}
