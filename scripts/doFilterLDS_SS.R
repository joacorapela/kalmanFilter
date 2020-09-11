require(ggplot2)
require(plotly)
require(htmlwidgets)
require(mvtnorm)
require(ini)

source("../src/filterLDS_SS.R")

processAll <- function() {
    simResNumber <- 73451420
    xlab <- "x"
    ylab <- "y"
    simFilenamePattern <- "results/%08d_simulation.RData"
    filterResMetaDataFilenamePattern <- "results/%08d_filteredSimulation.ini"
    filterResFilenamePattern <- "results/%08d_filteredSimulation.RData"
    filterResFigFilenamePattern <- "figures/%08d_filteredSimulation.%s"

    exit <- FALSE
    while(!exit) {
        filterResNumber <- sample(1e8, 1)
        filterResFilename <- sprintf(filterResFilenamePattern, filterResNumber)
        if(!file.exists(filterResFilename)) {
            exit <- TRUE
        }
    }
    filterResMetaDataFilename <- sprintf(filterResMetaDataFilenamePattern, filterResNumber)
    show(sprintf("Saving estimation results in: %s", filterResFilename))
    show(sprintf("Saving estimation meta data in: %s", filterResMetaDataFilename))

    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(simFilename))
    nObs <- ncol(simRes$x)

    filterRes <- filterLDS_SS(y=simRes$y, B=simRes$B, Q=simRes$Q, mu0=simRes$mu0, V0=simRes$V0, Z=simRes$Z, R=simRes$R)
    save(filterRes, file=filterResFilename)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simResNumber=simResNumber)
    write.ini(x=metaData, filepath=filterResMetaDataFilename)

    df <- data.frame(x=c(simRes$x[1,], simRes$y[1,], filterRes$xnn[1,1,]), y=c(simRes$x[2,], simRes$y[2,], filterRes$xnn[2,1,]), type=factor(c(rep("latent", nObs), rep("measurement", nObs), rep("filtered", nObs))))
    title <- sprintf("Log-Likelihood %f", filterRes$logLike)
    hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, simRes$x[1,], simRes$x[2,])
    hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, simRes$y[1,], simRes$y[2,])
    hoverTextFiltered <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, filterRes$xnn[1,1,], filterRes$xnn[2,1,])
    df <- cbind(df, c(hoverTextLatents, hoverTextObservations, hoverTextFiltered))
    colnames(df) <- c("x", "y", "type", "hoverText")
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
    fig <- fig %>% layout(title=title)
    filterResPNGFilename <- sprintf(filterResFigFilenamePattern, filterResNumber, "png")
    filterResHTMLFilename <- sprintf(filterResFigFilenamePattern, filterResNumber, "html")
    # orca(p=fig, file=filterResPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(filterResHTMLFilename)),basename(filterResHTMLFilename)))
    # print(fig)

    browser()
}

processAll()
