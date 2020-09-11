require(plotly)
require(ini)
require(htmlwidgets)

source("../src/simulateLDS.R")

processAll <- function() {
    simConfigNumber <- 2
    xlab <- "x"
    ylab <- "y"
    simConfigFilenamePattern <- "data/%08d_simulation_metaData.ini"
    simFilenamePattern <- "results/%08d_simulation.RData"
    simFigFilenamePattern <- "figures/%08d_simulation.%s"

    simConfigFilename <- sprintf(simConfigFilenamePattern, simConfigNumber)
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
    N <- as.numeric(simConfig$control_variables$N)

    # state transition
    Btmp <- eval(parse(text=simConfig$state_variables$B))
    B <- dt*Btmp + diag(nrow(Btmp))

    # state noise covariance
    Q <- eval(parse(text=simConfig$state_variables$Q))

    # initial state mean
    mu0 <- eval(parse(text=simConfig$initial_state_variables$mu0))

    # initial state covariance
    V0 <- eval(parse(text=simConfig$initial_state_variables$V0))

    # state-measurement transfer
    Z <- eval(parse(text=simConfig$measurements_variables$Z))

    # measurements noise covariance
    R <- eval(parse(text=simConfig$measurements_variables$R))

    res <- simulateLDS(N=N, B=B, Q=Q, mu0=mu0, V0=V0, Z=Z, R=R)
    simRes <- c(res, list(B=B, Q=Q, mu0=mu0, V0=V0, Z=Z, R=R))
    save(simRes, file=simFilename)

    hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:N, res$x[1,], res$x[2,])
    hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:N, res$y[1,], res$y[2,])
    df <- data.frame(t(cbind(res$x, res$y)))
    df <- cbind(df, c(rep("latent", N), rep("measurement", N)))
    df <- cbind(df, c(hoverTextLatents, hoverTextObservations))
    colnames(df) <- c("x", "y", "type", "hoverText")
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
    fig <- fig %>% add_annotations(x=c(res$x[1,1], res$x[1,N]), y=c(res$x[2,1], res$x[2,N]), text=c("start", "end"))
    simPNGFilename <- sprintf(simFigFilenamePattern, simResNumber, "png")
    simHTMLFilename <- sprintf(simFigFilenamePattern, simResNumber, "html")
    # orca(p=fig, file=simPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(simHTMLFilename)),basename(simHTMLFilename)))
    # print(fig)

    browser()
}

processAll()
