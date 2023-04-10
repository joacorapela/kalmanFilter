require(plotly)
require(htmlwidgets)
require(ini)
require(R.matlab)

source("../src/buildTwoStimHist.R")
source("../src/simulatePLDSwithOffsetsAndInputs.R")

processAll <- function() {
    xlab <- "x"
    ylab <- "y"
    simConfigFilename <- "data/00000010_simulation_metaData.ini"
    simResMetaDataFilenamePattern <- "results/%08d_simulation.ini"
    simResFilenamePattern <- "results/%08d_simulation.%s"
    # simFigFilenamePattern <- "figures/%08d_simulation.%s"

    simConfig <- read.ini(simConfigFilename)

    exit <- FALSE
    while(!exit) {
        simResNumber <- sample(1e8, 1)
        rSimResFilename <- sprintf(simResFilenamePattern, simResNumber, "RData")
        matSimResFilename <- sprintf(simResFilenamePattern, simResNumber, "mat")
        if(!file.exists(rSimResFilename) & !file.exists(matSimResFilename)) {
            exit <- TRUE
        }
    }
    simResMetaDataFilename <- sprintf(simResMetaDataFilenamePattern, simResNumber)
    show(sprintf("Simulation results in: %s and %s", rSimResFilename, matSimResFilename))

    # sampling rate
    sRate <- as.double(simConfig$control_variables$sRate)
    stateMem <- as.double(simConfig$control_variables$stateMem)
    obsMem <- as.double(simConfig$control_variables$obsMem)
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
    cHist <- buildTwoStimHist(stim=c, memory=stateMem)

    # state noise covariance
    Q <- eval(parse(text=simConfig$state_variables$Q))

    # initial state mean
    x00 <- eval(parse(text=simConfig$initial_state_variables$m0))

    # initial state covariance
    V00 <- eval(parse(text=simConfig$initial_state_variables$V0))

    # state-measurement transfer
    Z <- eval(parse(text=simConfig$measurements_variables$Z))

    # measurements offset
    a <- eval(parse(text=simConfig$measurements_variables$a))

    # measurments input matrix
    D <- eval(parse(text=simConfig$measurements_variables$D))

    # measurements inputs
    d <- as.matrix(read.table(simConfig$measurements_variables$dFilename))
    dHist <- buildTwoStimHist(stim=d, memory=obsMem)

    res <- simulatePLDSwithOffsetsAndInputs(B=B, u=u, C=C, c=cHist, Q=Q, x00=x00, V00=V00, Z=Z, a=a, D=D, d=dHist)
    simRes <- c(res, list(B=B, u=u, C=C, c=c, Q=Q, x00=x00, V00=V00, Z=Z, a=a, D=D, d=d))
    save(simRes, file=rSimResFilename)
    writeMat(con=matSimResFilename, x=simRes$x, y=simRes$y, z=simRes$z, B=simRes$B, u=simRes$u, C=simRes$C, c=c, Q=simRes$Q, x00=simRes$x00, V00=simRes$V00, Z=simRes$Z, a=simRes$a, D=simRes$D, d=simRes$d)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simConfigFilename=simConfigFilename)
    write.ini(x=metaData, filepath=simResMetaDataFilename)

#     nObs <- ncol(res$x)
#     hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, res$x[1,], res$x[2,])
#     hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, res$y[1,], res$y[2,])
#     df <- data.frame(t(cbind(res$x, res$y)))
#     df <- cbind(df, c(rep("latent", nObs), rep("measurement", nObs)))
#     df <- cbind(df, c(hoverTextLatents, hoverTextObservations))
#     colnames(df) <- c("x", "y", "type", "hoverText")
#     fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
#     fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
#     fig <- fig %>% add_annotations(x=c(res$x[1,1], res$x[1,nObs]), y=c(res$x[2,1], res$x[2,nObs]), text=c("start", "end"))
#     fig <- fig %>% add_annotations(x=c(res$x[1,750], res$x[1,1250]), y=c(res$x[2,750], res$x[2,1250]), text=c("stim start", "stim end"))
#     simPNGFilename <- sprintf(simFigFilenamePattern, simResNumber, "png")
#     simHTMLFilename <- sprintf(simFigFilenamePattern, simResNumber, "html")
#     orca(p=fig, file=simPNGFilename)
#     saveWidget(widget=fig, file=file.path(normalizePath(dirname(simHTMLFilename)),basename(simHTMLFilename)))
#     print(fig)

    browser()
}

processAll()
