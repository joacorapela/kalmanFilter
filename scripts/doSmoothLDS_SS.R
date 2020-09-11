
library(plotly)

source("../src/smoothLDS_SS.R")

processAll <- function() {
    filterResNumber <- 05797720
    simFilenamePattern <- "results/%08d_simulation.RData"
    filterResMetaDataFilenamePattern <- "results/%08d_filteredSimulation.ini"
    filterResFilenamePattern <- "results/%08d_filteredSimulation.RData"
    smoothResFilenamePattern <- "results/%08d_smoothedSimulation.RData"
    smoothResMetaDataFilenamePattern <- "results/%08d_smoothedSimulation.ini"
    smoothResFigFilenamePattern <- "figures/%08d_smoothedSimulation.%s"

    exit <- FALSE
    while(!exit) {
        smoothResNumber <- sample(1e8, 1)
        smoothResFilename <- sprintf(smoothResFilenamePattern, smoothResNumber)
        if(!file.exists(smoothResFilename)) {
            exit <- TRUE
        }
    }
    smoothResMetaDataFilename <- sprintf(smoothResMetaDataFilenamePattern, smoothResNumber)
    show(sprintf("Saving estimation results in: %s", smoothResFilename))
    show(sprintf("Saving estimation meta data in: %s", smoothResMetaDataFilename))

    filterResMetaDataFilename <- sprintf(filterResMetaDataFilenamePattern, filterResNumber)
    ini <- read.ini(filepath=filterResMetaDataFilename)
    simResNumber <- as.integer(ini$simulation_info$simResNumber)
    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(simFilename))

    filterResFilename <- sprintf(filterResFilenamePattern, filterResNumber)
    filterRes <- get(load(filterResFilename))

    smoothRes <- smoothLDS_SS(B=simRes$B, xnn=filterRes$xnn, Vnn=filterRes$Vnn, xnn1=filterRes$xnn1, Vnn1=filterRes$Vnn1, mu0=simRes$mu0, V0=simRes$V0)
    save(smoothRes, file=smoothResFilename)

    metaData <- list()
    metaData[["filter_info"]] <- list(filterResNumber=filterResNumber)
    write.ini(x=metaData, filepath=smoothResMetaDataFilename)

    N <- dim(smoothRes$xnN)[3]
    df <- data.frame(x=c(simRes$x[1,], simRes$y[1,], filterRes$xnn[1,1,], smoothRes$xnN[1,1,]), y=c(simRes$x[2,], simRes$y[2,], filterRes$xnn[2,1,], smoothRes$xnN[2,1,]), type=factor(c(rep("latent", N), rep("measurement", N), rep("filtered", N), rep("smoothed", N))))
    hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:N, simRes$x[1,], simRes$x[2,])
    hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:N, simRes$y[1,], simRes$y[2,])
    hoverTextFiltered <- sprintf("sample %d, x %.02f, y %.02f", 1:N, filterRes$xnn[1,1,], filterRes$xnn[2,1,])
    hoverTextSmoothed <- sprintf("sample %d, x %.02f, y %.02f", 1:N, smoothRes$xnN[1,1,], smoothRes$xnN[2,1,])
    df <- cbind(df, c(hoverTextLatents, hoverTextObservations, hoverTextFiltered, hoverTextSmoothed))
    colnames(df) <- c("x", "y", "type", "hoverText")
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
    # smoothResPNGFilename <- sprintf(smoothResFigFilenamePattern, smoothResNumber, "png")
    smoothResHTMLFilename <- sprintf(smoothResFigFilenamePattern, smoothResNumber, "html")
    # orca(p=fig, file=smoothResPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(smoothResHTMLFilename)),basename(smoothResHTMLFilename)))
    # print(fig)

    browser()
}

processAll()
