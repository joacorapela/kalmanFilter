computeAIC <- function(dsSSM) {
    computeNParams <- function(dsSSM) {
        nParams <- 0
        nParams <- nParams+length(dsSSM$B)
        nParams <- nParams+length(dsSSM$u)
        nParams <- nParams+length(dsSSM$C)
        if(dsSSM$covsConstraints$Q=="unconstrained") {
            nParams <- nParams+length(dsSSM$Q)
        } else {
            if(dsSSM$covsConstraints$Q=="diagonal and unequal") {
                nParams <- nParams+length(diag(dsSSM$Q))
            } else {
                stop(sprintf("Invalid covConstraint for Q: %s", dsSSM$covsConstraints$Q))
            }
        }
        nParams <- nParams+length(dsSSM$Z)
        nParams <- nParams+length(dsSSM$a)
        nParams <- nParams+length(dsSSM$D)
        if(dsSSM$covsConstraints$R=="unconstrained") {
            nParams <- nParams+length(dsSSM$R)
        } else {
            if(dsSSM$covsConstraints$R=="diagonal and unequal") {
                nParams <- nParams+length(diag(dsSSM$R))
            } else {
                stop(sprintf("Invalid covConstraint for R: %s", dsSSM$covsConstraints$R))
            }
        }
        nParams <- nParams+length(dsSSM$m0)
        if(dsSSM$covsConstraints$V0=="unconstrained") {
            nParams <- nParams+length(dsSSM$V0)
        } else {
            if(dsSSM$covsConstraints$V0=="diagonal and unequal") {
                nParams <- nParams+length(diag(dsSSM$V0))
            } else {
                stop(sprintf("Invalid covConstraint for V0: %s", dsSSM$covsConstraints$V0))
            }
        }
        return(nParams)
    }

    logLik <- dsSSM$logLik[length(dsSSM$logLik)]
    k <- computeNParams(dsSSM=dsSSM)
    AIC <- 2*(k-logLik)
    return(AIC)
}
