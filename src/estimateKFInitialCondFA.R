estimateKFInitialCondFA <- function(z, nFactors, rotation="varimax", scores="regression", control=list(trace=TRUE, nstart=5)) {
    res <- factanal(z, factors=nFactors, rotation=rotation, scores=scores, control=control)
    Z <- as(res$loadings, "matrix")
    latents <- res$scores
    RDiag <- res$uniqueness
    XnextT <- latents[2:nrow(latents),]
    XcurT <- latents[1:(nrow(latents)-1),]
    # lmRes <- lm(XnextT~0+XcurT)
    # B <- solve(a=t(XcurT)%*%XcurT, b=t(XcurT)%*%XnextT)
    B <- t(ginv(XcurT)%*%XnextT)
    return(list(B=B, Z=Z, RDiag=RDiag))
}
