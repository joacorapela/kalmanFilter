require(pcaMethods)
require(MASS)

estimateKFInitialCondPPCA <- function(z, nFactors) {
    # res <- factanal(~z, factors=nFactors, rotation="varimax", na.action=na.exclude)
    res <- pca(z, method="ppca", nPcs=nFactors)
    Z <- loadings(res)
    latents <- scores(res)
    XnextT <- latents[2:nrow(latents),]
    XcurT <- latents[1:(nrow(latents)-1),]
    # lmRes <- lm(XnextT~0+XcurT)
    # A <- solve(a=t(XcurT)%*%XcurT, b=t(XcurT)%*%XnextT)
    B <- t(ginv(XcurT)%*%XnextT)
    return(list(B=B, Z=Z))
}
