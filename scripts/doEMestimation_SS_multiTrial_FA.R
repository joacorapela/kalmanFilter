
require(MASS)
require(ini)
source("../src/emEstimationKF_SS_multiTrial.R")
source("../src/filterLDS_SS.R")
source("../src/smoothLDS_SS.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/estimateKFInitialCondPPCA.R")

processAll <- function() {
    estConfigNumber <- 3
    simResNumber <- 79772839
    simulationFilenamePattern <- "results/%s_simulation.RData"

    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResFilenamePattern <- "results/%08d_estimation.RData"
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
    nTrials <- as.numeric(estConfig$control_variables$nTrials)

    # initial state covariance
    V0Trial1 <- eval(parse(text=estConfig$initial_values$V0Trial1))
    M <- nrow(V0Trial1)
    V0 <- array(dim=c(nTrials, M, M))
    V0[1,,] <- V0Trial1
    for(r in 2:nTrials) {
        V0[r,,] <- eval(parse(text=estConfig$initial_values[[sprintf("V0Trial%d", r)]]))
    }

    # initial state mean
    if(tolower(estConfig$initial_values$m0Type)=="simulated") {
        m0 <- simRes$m0
    } else {
        if(tolower(estConfig$initial_values$m0Type)=="randomuniform") {
            m0Min <- as.double(estConfig$initial_value$m0Min)
            m0Max <- as.double(estConfig$initial_value$m0Max)
            m0 <- array(dim=c(nTrials, M))
            for(r in 1:nTrials) {
                m0[r,] <- runif(n=M, min=m0Min, max=m0Max)
            }
        } else {
            if(tolower(estConfig$initial_values$m0Type)=="given") {
                m0 <- array(dim=c(nTrials, M))
                for(r in 1:nTrials) {
                    m0[r,] <- eval(parse(text=estConfig$initial_values[[sprintf("m0Trial%d", r)]]))
                }
            } else {
                stop(sprintf("Invalid m0Type=%s", estConfig$initial_values$m0Type))
            }
        }
    }

    Q0 <- eval(parse(text=estConfig$initial_values$Q0))
    R0 <- eval(parse(text=estConfig$initial_values$R0))

    simRes <- get(load(simulationFilename))
    y <- simRes$y

    P <- dim(y)[2]
    N <- dim(y)[3]
    yConcatenated <- matrix(NA, nrow=P, ncol=nTrials*N)
    for(r in 1:nTrials) {
        yConcatenated[,(r-1)*N+(1:N)] <- y[r,,]
    }
    initialConds <- estimateKFInitialCondPPCA(z=t(as.matrix(yConcatenated)), nFactors=M)
    initialConds <- c(initialConds, list(Q0=Q0, R0=R0, m0=m0, V0=V0))

    estRes <- emEstimationKF_SS_multiTrial(y=y, B0=initialConds$B, Q0=initialConds$Q0, Z0=initialConds$Z, R0=initialConds$R0, m0=initialConds$m0, V0=initialConds$V0, maxIter=maxIter, tol=tol, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE))
    # estRes <- emEstimationKF_SS_multiTrial(y=y, B0=simRes$B, Q0=simRes$Q, Z0=simRes$Z, R0=simRes$R, m0=simRes$m0, V0=simRes$V0, maxIter=maxIter, tol=tol, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE))

    estRes <- c(estRes, list(initialConds=initialConds))
    save(file=estResFilename, estRes)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simResNumber=simResNumber)
    metaData[["estimation_config_info"]] <- list(estConfigNumber=estConfigNumber)
    write.ini(x=metaData, filepath=estResMetaDataFilename)

    browser()
}

processAll()
