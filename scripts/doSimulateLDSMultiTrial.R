require(ini)
source("../src/simulateLDS.R")

processAll <- function() {
    simConfigNumber <- 1
    simConfigFilenamePattern <- "data/%08d_simulation_metaData.ini"
    simResMetaDataFilenamePattern <- "results/%08d_simulation.ini"
    simResFilenamePattern <- "results/%08d_simulation.RData"

    simConfigFilename <- sprintf(simConfigFilenamePattern, simConfigNumber)
    simConfig <- read.ini(simConfigFilename)

    exit <- FALSE
    while(!exit) {
        simResNumber <- sample(1e8, 1)
        simResFilename <- sprintf(simResFilenamePattern, simResNumber)
        if(!file.exists(simResFilename)) {
            exit <- TRUE
        }
    }
    simResMetaDataFilename <- sprintf(simResMetaDataFilenamePattern, simResNumber)
    show(sprintf("Simulation results in: %s", simResFilename))

    # sampling rate
    sRate <- as.double(simConfig$control_variables$sRate)
    dt <- 1/sRate

    # number of observations
    N <- as.numeric(simConfig$control_variables$N)

    # number of trials
    nTrials <- as.numeric(simConfig$control_variables$nTrials)

    # state transition
    Btmp <- eval(parse(text=simConfig$state_variables$B))
    B <- dt*Btmp + diag(nrow(Btmp))

    # state noise covariance
    Q <- eval(parse(text=simConfig$state_variables$Q))

    # initial state mean
    m0Trial1 <- eval(parse(text=simConfig$initial_state_variables$m0Trial1))
    stateDim <- length(m0Trial1)
    m0 <- array(dim=c(nTrials, stateDim))
    m0[1,] <- m0Trial1
    for(i in 2:nTrials) {
        m0[i,] <- eval(parse(text=simConfig$initial_state_variables[[sprintf("m0Trial%d", i)]]))
    }

    # initial state covariance
    V0Trial1 <- eval(parse(text=simConfig$initial_state_variables$V0Trial1))
    V0 <- array(dim=c(nTrials, stateDim, stateDim))
    V0[1,,] <- V0Trial1
    for(i in 2:nTrials) {
        V0[i,,] <- eval(parse(text=simConfig$initial_state_variables[[sprintf("V0Trial%d", i)]]))
    }

    # state-measurement transfer
    Z <- eval(parse(text=simConfig$measurements_variables$Z))
    
    # measurements noise covariance
    R <- eval(parse(text=simConfig$measurements_variables$R))

    P <- nrow(Z)
    M <- ncol(Z)
    x <- array(dim=c(nTrials, M, N))
    y <- array(dim=c(nTrials, P, N))
    for(i in 1:nTrials) {
        res <- simulateLDS(N=N, B=B, Q=Q, m0=m0[i,], V0=V0[i,,], Z=Z, R=R)
        x[i,,] <- res$x
        y[i,,] <- res$y
    }
    simRes <- list(x=x, y=y, B=B, Q=Q, m0=m0, V0=V0, Z=Z, R=R)
    save(simRes, file=simResFilename)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simConfigNumber=simConfigNumber)
    write.ini(x=metaData, filepath=simResMetaDataFilename)

    browser()
}

processAll()
