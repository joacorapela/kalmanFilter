
source("../src/filterLDS_R.R")

getParamsVecFromParamsList <- function(paramsList) {
    paramsVec <- c()
    for(i in 1:length(paramsList)) {
        paramsVec <- c(paramsVec, as.vector(paramsList[[i]]))
    }
    return(paramsVec)
}
getParamsListFromParamsVec <- function(paramsVec, paramsList0) {
    paramsList <- list()
    index <- 0
    for(i in 1:length(paramsList0)) {
        item0 <- paramsList0[[i]]
        lenItem0 <- length(as.vector(item0))
        listItem <- list()
        listItem[[names(paramsList0)[i]]] <- matrix(paramsVec[index+(1:lenItem0)], ncol=ncol(item0))
        index <- index+lenItem0
        paramsList <- c(paramsList, listItem)
    }
    return(paramsList)
}

processAll <- function() {
    dObs <- 10
    stateType0 <- "init00"
    simulationFilenamePattern <- "results/simulationSepiDObs%02d.RData"
    estResultsFilenamePattern <- "results/sepiDObs%02dMyEM_estimation.RData"

    simulationFilename <- sprintf(simulationFilenamePattern, dObs)
    estResultsFilename <- sprintf(estResultsFilenamePattern, dObs)
    simRes <- get(load(simulationFilename))
    ys <- simRes$x
    loadRes <- get(load(estResultsFilename))
    srQ <- chol(loadRes$Q)
    srR <- chol(loadRes$R)
    srV0 <- chol(loadRes$V0)
    paramsList0 <- list(B=loadRes$B, srQ=srQ, Z=loadRes$Z, srR=srR, x0=loadRes$x0, srV0=srV0)
    paramsVec0 <- getParamsVecFromParamsList(paramsList=paramsList0)
    f <- function(paramsVec) {
        paramsList <- getParamsListFromParamsVec(paramsVec=paramsVec, paramsList0=paramsList0)
        Q <- t(paramsList$srQ)%*%paramsList$srQ
        R <- t(paramsList$srR)%*%paramsList$srR
        V0 <- t(paramsList$srV0)%*%paramsList$srV0
        kfRes <- filterLDS_R(B=paramsList$B, Z=paramsList$Z, x0=paramsList$x0, V0=V0, stateType0=stateType0, Q=Q, R=R, ys=ys)
        answer <- -kfRes$logLike
        return(answer)
    }
    est <- optim(par=paramsVec0, fn=f, gr=NULL, method='BFGS', hessian=TRUE, control=list(trace=1, REPORT=1))

    browser()
}

processAll()

