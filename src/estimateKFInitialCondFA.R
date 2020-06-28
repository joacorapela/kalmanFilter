estimateKFInitialCondFA <- function(z, nFactors, rotation="varimax", scores="regression") {
    res <- factanal(z, factors=nFactors, rotation=rotation, scores=scores, trace=TRUE)
    C <- as(res$loadings, "matrix")
    latents <- res$scores
    sigmaDiag <- res$uniqueness
    XnextT <- latents[2:nrow(latents),]
    XcurT <- latents[1:(nrow(latents)-1),]
    # lmRes <- lm(XnextT~0+XcurT)
    # A <- solve(a=t(XcurT)%*%XcurT, b=t(XcurT)%*%XnextT)
    A <- t(ginv(XcurT)%*%XnextT)
    return(list(A=A, C=C, sigmaDiag=sigmaDiag))
}
