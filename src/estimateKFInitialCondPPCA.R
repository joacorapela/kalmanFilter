estimateKFInitialCondPPCA <- function(z, nFactors) {
    # res <- factanal(~z, factors=nFactors, rotation="varimax", na.action=na.exclude)
    res <- pca(z, method="ppca", nPcs=nFactors)
    C <- loadings(res)
    latents <- scores(res)
    Sigma <- diag(nrow(C))*1e-2
    XnextT <- latents[2:nrow(latents),]
    XcurT <- latents[1:(nrow(latents)-1),]
    # lmRes <- lm(XnextT~0+XcurT)
    # A <- solve(a=t(XcurT)%*%XcurT, b=t(XcurT)%*%XnextT)
    A <- t(ginv(XcurT)%*%XnextT)
    Gamma <- diag(nrow(A))*1e-2
    return(list(A=A, Gamma=Gamma, C=C, Sigma=Sigma))
}
