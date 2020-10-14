require(plotly)
require(htmlwidgets)
require(R.matlab)
require(ini)

source("../src/simulateLDSwithOffsetsAndInputs.R")

processAll <- function() {
    simConfigNumber <- 4
    xlab <- "x"
    ylab <- "y"
    simConfigFilenamePattern <- "data/%08d_simulation_metaData.ini"
    simResMetaDataFilenamePattern <- "results/%08d_simulation.ini"
    simResFilenamePattern <- "results/%08d_simulation.RData"
    simFigFilenamePattern <- "figures/%08d_simulation.%s"

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

    # state transition
    Btmp <- eval(parse(text=simConfig$state_variables$B))
    B <- dt*Btmp + diag(nrow(Btmp))

    # state offset
    u <- eval(parse(text=simConfig$state_variables$u))
    u <- u*dt

    # state input matrix
    C <- eval(parse(text=simConfig$state_variables$C))

    # state inputs
    c <- as.matrix(read.table(simConfig$state_variables$cFilename))
    c <- c*dt

    # state noise covariance
    Q <- eval(parse(text=simConfig$state_variables$Q))

    # initial state mean
    x00 <- eval(parse(text=simConfig$initial_state_variables$x00))

    # initial state covariance
    V00 <- eval(parse(text=simConfig$initial_state_variables$V00))

    # state-measurement transfer
    Z <- eval(parse(text=simConfig$measurements_variables$Z))

    # measurements offset
    a <- eval(parse(text=simConfig$measurements_variables$a))

    # measurments input matrix
    D <- eval(parse(text=simConfig$measurements_variables$D))

    # measurements inputs
    d <- as.matrix(read.table(simConfig$measurements_variables$dFilename))

    # measurements noise covariance
    R <- eval(parse(text=simConfig$measurements_variables$R))

    res <- simulateLDSwithOffsetsAndInputs(B=B, u=u, C=C, c=c, Q=Q, x00=x00, V00=V00, Z=Z, a=a, D=D, d=d, R=R)
    simRes <- c(res, list(B=B, u=u, C=C, c=c, Q=Q, x00=x00, V00=V00, Z=Z, a=a, D=D, d=d, R=R))
    save(simRes, file=simResFilename)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simConfigNumber=simConfigNumber)
    write.ini(x=metaData, filepath=simResMetaDataFilename)

    nObs <- ncol(res$x)
    hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, res$x[1,], res$x[2,])
    hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, res$y[1,], res$y[2,])
    df <- data.frame(t(cbind(res$x, res$y)))
    df <- cbind(df, c(rep("latent", nObs), rep("measurement", nObs)))
    df <- cbind(df, c(hoverTextLatents, hoverTextObservations))
    colnames(df) <- c("x", "y", "type", "hoverText")
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
    fig <- fig %>% add_annotations(x=c(res$x[1,1], res$x[1,nObs]), y=c(res$x[2,1], res$x[2,nObs]), text=c("start", "end"))
    fig <- fig %>% add_annotations(x=c(res$x[1,750], res$x[1,1250]), y=c(res$x[2,750], res$x[2,1250]), text=c("stim start", "stim end"))
    simPNGFilename <- sprintf(simFigFilenamePattern, simResNumber, "png")
    simHTMLFilename <- sprintf(simFigFilenamePattern, simResNumber, "html")
    orca(p=fig, file=simPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(simHTMLFilename)),basename(simHTMLFilename)))
    print(fig)

    browser()
}

processAll()
