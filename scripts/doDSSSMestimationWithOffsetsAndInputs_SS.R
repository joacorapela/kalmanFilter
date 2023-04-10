
require(ini)
# source("../src/buildStimHist.R")
source("../src/getInitialConds.R")
source("../src/computeAIC.R")
source("../src/emEstimationKF_SS_withOffsetsAndInputs.R")
source("../src/filterLDS_SS_withOffsetsAndInputs.R")
source("../src/smoothLDS_SS.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/estimateKFInitialCondPPCA.R")

processAll <- function() {
    stateDim <- 2
    # simResNumber <- 32201051
    simResNumber <- 82902398
    stateInputMemorySamples <- 0
    obsInputMemorySamples <- 0
    initialCondMethod <- "FA"
    nStartFA <- 5
    propTrain <- 0.75
    estConfigFilename <- "data/00000001_estimation_metaData.ini"
    modelsLogFilename <- "log/circle.log"

    simResMetaDataFilenamePattern <- "results/%08d_simulation.ini"
    simResFilenamePattern <- "results/%s_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"

    simResMetaDataFilename <- sprintf(simResMetaDataFilenamePattern, simResNumber)
    simResMetaData <- read.ini(simResMetaDataFilename)
    simConfigNumber <- as.numeric(simResMetaData$simulation_info$simConfigNumber)
    exit <- FALSE
    while(!exit) {
        estNumber <- sample(1e8, 1)
        estResFilename <- sprintf(estResFilenamePattern, estNumber)
        if(!file.exists(estResFilename)) {
            exit <- TRUE
        }
    }
    estResMetaDataFilename <- sprintf(estResMetaDataFilenamePattern, estNumber)
    show(sprintf("Saving estimation results in: %s", estResFilename))
    show(sprintf("Saving estimation meta data in: %s", estResMetaDataFilename))

    estConfig <- read.ini(estConfigFilename)

    stateCovType <- estConfig$covariance_type$states
    initialStateCovType <- estConfig$covariance_type$initialStates
    obsCovType <- estConfig$covariance_type$observations
    covsConstraints <- list(V0=initialStateCovType, Q=stateCovType, R=obsCovType)                                                            
    # EM convergence tolerance
    tol <- as.double(estConfig$control_variables$tol)
    maxIter <- as.numeric(estConfig$control_variables$maxIter)
    maxIter <- maxIter+1 # to make a fair comparison with MARSS
    minIter <- as.numeric(estConfig$control_variables$minIter)

    simResFilename <- sprintf(simResFilenamePattern, simResNumber)
    simRes <- get(load(simResFilename))
    y <- simRes$y
    c <- simRes$c
    d <- simRes$d
    obsDim <- nrow(y)
    nStateInputs <- nrow(c)
    nObsInputs <- nrow(d)
    
    initialConds <- getInitialConds(initialValues=estConfig$initial_values, stateDim=stateDim, obsDim=obsDim, stateInputMemorySamples=stateInputMemorySamples, obsInputMemorySamples=obsInputMemorySamples, nStateInputs=nStateInputs, nObsInputs=nObsInputs)

    stateMem <- ncol(initialConds$C)
    # cHist <- buildStimHist(stim=c, memory=stateMem)
    cHist <- c
    dim(cHist) <- c(nrow(cHist), 1, ncol(cHist))
    obsMem <- ncol(initialConds$D)
    # dHist <- buildStimHist(stim=d, memory=obsMem)
    dHist <- d 
    dim(dHist) <- c(nrow(dHist), 1, ncol(dHist))

    N <- ncol(y)
    nTrain <- round(N*propTrain)
    yTrain <- y[,1:nTrain]
    yValidate <- y[,(nTrain+1):N]
    cTrain <- cHist[,,1:nTrain, drop=FALSE]
    cValidate <- cHist[,,(nTrain+1):N, drop=FALSE]
    dTrain <- dHist[,,1:nTrain, drop=FALSE]
    dValidate <- dHist[,,(nTrain+1):N, drop=FALSE]

    dataForEstInitialCond <- t(as.matrix(yTrain))
    if(initialCondMethod=="FA") {
        controlFA <- list(trace=TRUE, nstart=nStartFA)
        estRes <- estimateKFInitialCondFA(z=dataForEstInitialCond, nFactors=stateDim, control=controlFA)
        initialConds <- c(list(B=estRes$B, Z=estRes$Z, R=diag(estRes$RDiag)), initialConds)
    } else {
        if(initialCondMethod=="PPCA") {
            estRes <- estimateKFInitialCondPPCA(z=dataForEstInitialCond, nFactors=stateDim)
            initialConds <- c(list(B=estRes$B, Z=estRes$Z), initialConds)
        } else {
            stop(sprintf("Invalid initialCondMethod=%s", initialCondMethod))
        }
    }

    startTime <- proc.time()[3]
    dsSSM <- emEstimationKF_SS_withOffsetsAndInputs(y=yTrain, c=cTrain, d=dTrain, B0=initialConds$B, u0=initialConds$u, C0=initialConds$C, Q0=initialConds$Q, Z0=initialConds$Z, a0=initialConds$a, D0=initialConds$D, R0=initialConds$R, m0=initialConds$m0, V0=initialConds$V0, minIter=minIter, maxIter=maxIter, tol=tol, varsToEstimate=list(m0=TRUE, V0=TRUE, B=TRUE, u=TRUE, C=TRUE, Q=TRUE, Z=TRUE, a=TRUE, D=FALSE, R=TRUE), covsConstraints=covsConstraints)
    elapsedTime <- proc.time()[3]-startTime
    elapsedTime <- unname(elapsedTime)

    show(sprintf("Elapsed time: %f", elapsedTime))

    AIC <- computeAIC(dsSSM=dsSSM)
    kf <- filterLDS_SS_withOffsetsAndInputs(y=yValidate, c=cValidate, d=dValidate, B=dsSSM$B, u=dsSSM$u, C=dsSSM$C, Q=dsSSM$Q, m0=dsSSM$xNN, V0=dsSSM$VNN, Z=dsSSM$Z, a=dsSSM$a, D=dsSSM$D, R=dsSSM$R)
    cvLogLike <- kf$logLike

    ks <- smoothLDS_SS(B=dsSSM$B, xnn=kf$xnn, Vnn=kf$Vnn, xnn1=kf$xnn1, Vnn1=kf$Vnn1, m0=dsSSM$xNN, V0=dsSSM$VNN)
    estMeanYValidate <- dsSSM$Z%*%ks$xnN[,1,]+as.vector(dsSSM$a)+dsSSM$D%*%dValidate[,1,]
    mses <- rowMeans((yValidate-estMeanYValidate)^2)
    yValidateVars <- apply(yValidate, 1, var)
    pevs <- 1-(mse/yValidateVars)
    meanPEVs <- mean(pevs)
    show(sprintf('Mean percentage of explained variance %f', meanPEVs));

    logMessage <- sprintf("%d, %d, %d, %d, %d, %s, %f, %f, %f, %f, %f\n", estNumber, stateDim, obsDim, stateInputMemorySamples, obsInputMemorySamples, initialCondMethod, dsSSM$logLik[length(dsSSM$logLik)], AIC, cvLogLike, meanPEVs, elapsedTime)
    show(logMessage)
    cat(logMessage, file=modelsLogFilename, append=TRUE)

    estRes <- list(dsSSM=dsSSM, y=yTrain, c=cTrain, d=dTrain, initialConds=initialConds, AIC=AIC, cvLogLike=cvLogLike, pevs=pevs)
    save(file=estResFilename, estRes)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simResNumber=simResNumber)
    metaData[["estimation_config_info"]] <- list(estConfigFilename=estConfigFilename)
    metaData[["estimation_summary"]] <- list(logLik=dsSSM$logLik[length(dsSSM$logLik)], cvLogLike=cvLogLike, elapsedTime=elapsedTime)
    write.ini(x=metaData, filepath=estResMetaDataFilename)

    browser()
}

processAll()
