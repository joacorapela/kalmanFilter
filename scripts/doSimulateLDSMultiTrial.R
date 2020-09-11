require(ini)
source("../src/simulateLDS.R")

processAll <- function() {
    simConfigFilename <- "data/00000001_simulation_metaData.ini"
    simFilenamePattern <- "results/%08d_simulation.RData"

    simConfig <- read.ini(simConfigFilename)

    exit <- FALSE
    while(!exit) {
        simResNumber <- sample(1e8, 1)
        simFilename <- sprintf(simFilenamePattern, simResNumber)
        if(!file.exists(simFilename)) {
            exit <- TRUE
        }
    }
    show(sprintf("Simulation results in: %s", simFilename))

    # sampling rate
    sRate <- as.double(simConfig$control_variables$sRate)
    dt <- 1/sRate

    # number of observations
    nObs <- as.numeric(simConfig$control_variables$nObs)

    # number of trials
    nTrials <- as.numeric(simConfig$control_variables$nTrials)

    # state transition
    Btmp <- eval(parse(text=simConfig$state_variables$B))
    B <- dt*Btmp + diag(nrow(Btmp))

    # state noise covariance
    Q <- eval(parse(text=simConfig$state_variables$Q))

    # initial state mean
    x00Trial1 <- eval(parse(text=simConfig$initial_state_variables$x00Trial1))
    stateDim <- length(x00Trial1)
    x00 <- array(dim=c(nTrials, stateDim))
    x00[1,] <- x00Trial1
    for(i in 2:nTrials) {
        x00[i,] <- eval(parse(text=simConfig$initial_state_variables[[sprintf("x00Trial%d", i)]]))
    }

    # initial state covariance
    V00Trial1 <- eval(parse(text=simConfig$initial_state_variables$V00Trial1))
    V00 <- array(dim=c(nTrials, stateDim, stateDim))
    V00[1,,] <- V00Trial1
    for(i in 2:nTrials) {
        V00[i,,] <- eval(parse(text=simConfig$initial_state_variables[[sprintf("V00Trial%d", i)]]))
    }

    # state-measurement transfer
    Z <- eval(parse(text=simConfig$measurements_variables$Z))
    
    # measurements noise covariance
    R <- eval(parse(text=simConfig$measurements_variables$R))

    dObs <- nrow(Z)
    dLat <- ncol(Z)
    x <- array(dim=c(nTrials, dLat, nObs))
    y <- array(dim=c(nTrials, dObs, nObs))
    for(i in 1:nTrials) {
        res <- simulateLDS(nObs=nObs, B=B, Q=Q, x00=x00[i,], V00=V00[i,,], Z=Z, R=R)
        x[i,,] <- res$x
        y[i,,] <- res$y
    }
    simRes <- list(x=x, y=y, B=B, Q=Q, x00=x00, V00=V00, Z=Z, R=R)
    save(simRes, file=simFilename)

    browser()
}

processAll()
