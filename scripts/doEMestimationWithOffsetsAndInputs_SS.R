
require(MASS)
require(ini)
# require(reshape2)
source("../src/emEstimationKF_SS_withOffsetsAndInputs.R")
source("../src/filterLDS_SS_withOffsetsAndInputs.R")
source("../src/smoothLDS_SS_withOffsetsAndInputs.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/estimateKFInitialCondPPCA.R")

processAll <- function() {
    estConfigNumber <- 3
    simResNumber <- 58388369
    simConfigFilenamePattern <- "data/%08d_simulation_metaData.ini"
    simResMetaDataFilenamePattern <- "results/%08d_simulation.ini"
    simResFilenamePattern <- "results/%s_simulation.RData"

    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"
    estFilenamePattern <- "results/%08d_estimation.RData"

    simResMetaDataFilename <- sprintf(simResMetaDataFilenamePattern, simResNumber)
    simResMetaData <- read.ini(simResMetaDataFilename)
    simConfigNumber <- as.numeric(simResMetaData$simulation_info$simConfigNumber)
    simConfigFilename <- sprintf(simConfigFilenamePattern, simConfigNumber)
    simConfig <- read.ini(simConfigFilename)

    sRate <- as.double(simConfig$control_variables$sRate)
    dt <- 1/sRate
    c <- as.matrix(read.table(simConfig$state_variables$cFilename))
    c <- c*dt
    dim(c) <- c(1, 1, length(c))
    d <- as.matrix(read.table(simConfig$measurements_variables$dFilename))
    dim(d) <- c(1, 1, length(d))

    simResFilename <- sprintf(simResFilenamePattern, simResNumber)

    exit <- FALSE
    while(!exit) {
        estResNumber <- sample(1e8, 1)
        estResFilename <- sprintf(estResFilenamePattern, estResNumber)
        if(!file.exists(estResFilename)) {
            exit <- TRUE
        }
    }
    estResMetaDataFilename <- sprintf(estResMetaDataFilenamePattern, estResNumber)
    show(sprintf("Saving estimation results in: %s", estResFilename))
    show(sprintf("Saving estimation meta data in: %s", estResMetaDataFilename))

    estConfigFilename <- sprintf(estConfigFilenamePattern, estConfigNumber)
    estConfig <- read.ini(estConfigFilename)

    # EM convergence tolerance
    tol <- as.double(estConfig$control_variables$tol)
    maxIter <- as.numeric(estConfig$control_variables$maxIter)
    maxIter <- maxIter+1 # to make a fair comparison with MARSS

    V0 <- eval(parse(text=estConfig$initial_values$V0))
    u0 <- eval(parse(text=estConfig$initial_values$u0))
    C0 <- eval(parse(text=estConfig$initial_values$C0))
    Q0 <- eval(parse(text=estConfig$initial_values$Q0))
    a0 <- eval(parse(text=estConfig$initial_values$a0))
    D0 <- eval(parse(text=estConfig$initial_values$D0))
    R0 <- eval(parse(text=estConfig$initial_values$R0))

    M <- nrow(V0)

    simRes <- get(load(simResFilename))
    y <- simRes$y
    P <- nrow(y)
    if(tolower(estConfig$initial_values$m0)=="simulated") {
        m0 <- simRes$m0
    } else {
        if(tolower(estConfig$initial_values$m0)=="randomuniform") {
            m0Min <- as.double(estConfig$initial_value$m0Min)
            m0Max <- as.double(estConfig$initial_value$m0Max)
            m0 <- runif(n=M, min=m0Min, max=m0Max)
        } else {
            m0 <- eval(parse(text=estConfig$initial_values$m0))
        }
    }

    # initialConds <- estimateKFInitialCondFA(z=t(as.matrix(y)), nFactors=M)
    initialConds <- estimateKFInitialCondPPCA(z=t(as.matrix(y)), nFactors=M)
    initialConds <- c(initialConds, list(u0=u0, C0=C0, Q0=Q0, a0=a0, D0=D0, R0=R0, m0=m0, V0=V0))

    startTime <- proc.time()[3]
    dsSSM <- emEstimationKF_SS_withOffsetsAndInputs(y=y, c=c, d=d, B0=initialConds$B, u0=initialConds$u0, C0=initialConds$C0, Q0=initialConds$Q0, Z0=initialConds$Z, a0=initialConds$a0, D0=initialConds$D0, R0=initialConds$R0, m0=initialConds$m0, V0=initialConds$V0, maxIter=maxIter, tol=tol, varsToEstimate=list(m0=TRUE, V0=TRUE, B=TRUE, u=TRUE, C=TRUE, Q=TRUE, Z=TRUE, a=TRUE, D=FALSE, R=TRUE))
    elapsedTime <- proc.time()[3]-startTime
    elapsedTime <- unname(elapsedTime)

    show(sprintf("Elapsed time: %f", elapsedTime))

    estRes <- list(dsSSM=dsSSM, y=y, c=c, d=d, initialConds=initialConds)
    save(file=estResFilename, estRes)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simResNumber=simResNumber)
    metaData[["estimation_config_info"]] <- list(estConfigNumber=estConfigNumber)
    metaData[["estimation_summary"]] <- list(logLik=dsSSM$logLik[length(dsSSM$logLik)], elapsedTime=elapsedTime)
    write.ini(x=metaData, filepath=estResMetaDataFilename)

    browser()
}

processAll()
