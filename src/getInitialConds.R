getInitialConds <- function(initialValues, stateDim, obsDim, stateInputMemorySamples, obsInputMemorySamples, nStateInputs, nObsInputs) {
    if(tolower(initialValues$m0)=="random uniform") {
        m0Min <- as.double(initialValues$m0Min)
        m0Max <- as.double(initialValues$m0Max)
        m0 <- runif(n=stateDim, min=m0Min, max=m0Max)
    } else {
        if(tolower(initialValues$m0)=="unconstrained") {
            m0 <- eval(parse(text=initialValues$m0))
        } else {
            stop(sprintf("Invalid initial value for m0: %s", initialValues$m0))
        }
    }
    if(tolower(initialValues$V0)=="diagonal and equal") {
        value <- as.double(initialValues$V0Value)
        V0 <- diag(rep(value, times=stateDim))
    } else {
        if(tolower(initialValues$V0)=="unconstrained") {
            V0 <- eval(parse(text=initialValues$V0))
        } else {
            stop(sprintf("Invalid initial value for V0: %s", initialValues$V0))
        }
    }
    if(tolower(initialValues$u0)=="equal") {
        value <- as.double(initialValues$u0Value)
        u0 <- rep(value, times=stateDim)
    } else {
        if(tolower(initialValues$u0)=="unconstrained") {
            u0 <- eval(parse(text=initialValues$u0))
        } else {
            stop(sprintf("Invalid initial value for u0: %s", initialValues$u0))
        }
    }
    if(!is.nan(stateInputMemorySamples)) {
        if(tolower(initialValues$C0)=="equal") {
            value <- as.double(initialValues$C0Value)
            C0 <- matrix(value, nrow=stateDim, ncol=(1+stateInputMemorySamples)*nStateInputs)
        } else {
            if(tolower(initialValues$C0)=="unconstrained") {
                C0 <- eval(parse(text=initialValues$C0))
            } else {
                stop(sprintf("Invalid initial value for C0: %s", initialValues$C0))
            }
        }
    }
    if(tolower(initialValues$Q0)=="diagonal and equal") {
        value <- as.double(initialValues$Q0Value)
        Q0 <- diag(rep(value, times=stateDim))
    } else {
        if(tolower(initialValues$Q0)=="unconstrained") {
            Q0 <- eval(parse(text=initialValues$Q0))
        } else {
            stop(sprintf("Invalid initial value for Q0: %s", initialValues$Q0))
        }
    }
    if(tolower(initialValues$a0)=="equal") {
        value <- as.double(initialValues$a0Value)
        a0 <- rep(value, times=obsDim)
    } else {
        if(tolower(initialValues$a0)=="unconstrained") {
            a0 <- eval(parse(text=initialValues$a0))
        } else {
            stop(sprintf("Invalid initial value for a0: %s", initialValues$a0))
        }
    }
    if(!is.nan(obsInputMemorySamples)) {
        if(tolower(initialValues$D0)=="equal") {
            value <- as.double(initialValues$D0Value)
            D0 <- matrix(value, nrow=obsDim, ncol=(1+obsInputMemorySamples)*nObsInputs)
        } else {
            if(tolower(initialValues$D0)=="unconstrained") {
                D0 <- eval(parse(text=initialValues$D0))
            } else {
                stop(sprintf("Invalid initial value for D0: %s", initialValues$D0))
            }
        }
    }
    if(tolower(initialValues$R0)=="diagonal and equal") {
        value <- as.double(initialValues$R0Value)
        R0 <- diag(rep(value, times=obsDim))
    } else {
        if(tolower(initialValues$R0)=="unconstrained") {
            R0 <- eval(parse(text=initialValues$R0))
        } else {
            stop(sprintf("Invalid initial value for R0: %s", initialValues$R0))
        }
    }
    answer <- list(m0=m0, V0=V0, u=u0, Q=Q0, a=a0, R=R0)
    if(!is.nan(stateInputMemorySamples)) {
        answer <- c(answer, list(C=C0))
    }
    if(!is.nan(obsInputMemorySamples)) {
        answer <- c(answer, list(D=D0))
    }
    return(answer)
}
