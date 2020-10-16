
require(MASS)
require(reshape2)
require(ini)
source("../src/emEstimationKF_SS.R")
source("../src/filterLDS_SS.R")
source("../src/smoothLDS_SS.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/estimateKFInitialCondPPCA.R")

processAll <- function() {
    estConfigNumber <- 2
    # simResNumber <- 95498373
    simResNumber <- 21919562
    simulationFilenamePattern <- "results/%s_simulation.RData"

    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"
    estFilenamePattern <- "results/%08d_estimation.RData"

    simulationFilename <- sprintf(simulationFilenamePattern, simResNumber)

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

    V0 <- eval(parse(text=estConfig$initial_values$V0))
    M <- nrow(V0)
    Q0 <- eval(parse(text=estConfig$initial_values$Q0))
    R0 <- eval(parse(text=estConfig$initial_values$R0))

    simRes <- get(load(simulationFilename))
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
    initialConds <- c(initialConds, list(Q0=Q0, R0=R0, m0=m0, V0=V0))

    estRes <- emEstimationKF_SS(y=y, B0=initialConds$B, Q0=initialConds$Q0, Z0=initialConds$Z, R0=initialConds$R0, m0=initialConds$m0, V0=initialConds$V0, maxIter=maxIter, tol=tol, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE))

    estRes <- c(estRes, list(initialConds=initialConds))
    save(file=estResFilename, estRes)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simResNumber=simResNumber)
    metaData[["estimation_config_info"]] <- list(estConfigNumber=estConfigNumber)
    write.ini(x=metaData, filepath=estResMetaDataFilename)

    browser()
}

processAll()
