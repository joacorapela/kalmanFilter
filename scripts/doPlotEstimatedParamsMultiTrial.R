
require(plotly)
require(ini)
require(RColorBrewer)
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    # estResNumber <- 99485472
    # estResNumber <- 67957061
    # estResNumber <- 93973929
    # estResNumber <- 66860571
    # estResNumber <- 54015032
    # estResNumber <- 60570612
    estResNumber <- 89084189
    simFilenamePattern <- "results/%08d_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"

    estResFilename <- sprintf(estResFilenamePattern, estResNumber)
    estRes <- get(load(estResFilename))
    initialConds <- estRes$initialConds

    estResMetaDataFilename <- sprintf(estResMetaDataFilenamePattern, estResNumber)
    ini <- read.ini(filepath=estResMetaDataFilename)
    simResNumber <- as.integer(ini$simulation_info$simResNumber)
    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(simFilename))

    figFilenamePattern <- "figures/%08d_logLike.%s"
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, "html")
    df <- data.frame(x=1:length(estRes$logLike), y=estRes$logLike)
    fig <- plot_ly(data=df, x=~x, y=~y, type='scatter', mode='lines+markers')
    fig <- fig %>% layout(xaxis = list(title="Iteration"), yaxis = list (title="Log Likelihood"))
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # orca(p=fig, file=pngFigFilename)
    # print(fig)

    figFilenamePattern <- "figures/%08d_B.%s"
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(true=simRes$B, initial=initialConds$B, estimated=estRes$B, title="B")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # orca(p=fig, file=pngFigFilename)
    # print(fig)

    figFilenamePattern <- "figures/%08d_Z.%s"
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(true=simRes$Z, initial=initialConds$Z, estimated=estRes$Z, title="Z")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # orca(p=fig, file=pngFigFilename)
    # print(fig)

    figFilenamePattern <- "figures/%08d_Q.%s"
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(true=simRes$Q, initial=initialConds$Q, estimated=estRes$Q, title="Q")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # orca(p=fig, file=pngFigFilename)
    # print(fig)

    figFilenamePattern <- "figures/%08d_R.%s"
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, "html")
    fig <- getPlotTrueInitialAndEstimatedMatrices(true=simRes$R, initial=initialConds$R, estimated=estRes$R, title="R")
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # orca(p=fig, file=pngFigFilename)
    # print(fig)

    figFilenamePattern <- "figures/%08d_trial%03d_m0.%s"
    nTrials <- nrow(estRes$m0)
    for(r in 1:nTrials) {
        pngFigFilename <- sprintf(figFilenamePattern, estResNumber, r, "png")
        htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, r, "html")
        fig <- getPlotTrueInitialAndEstimatedVectors(true=simRes$m0[r,], initial=initialConds$m0[r,], estimated=estRes$m0[r,], title=sprintf("m0 for trial %03d", r))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
        # orca(p=fig, file=pngFigFilename)
        # print(fig)
    }

    figFilenamePattern <- "figures/%08d_trial%03d_V0.%s"
    nTrials <- dim(estRes$V0)[1]
    for(r in 1:nTrials) {
        pngFigFilename <- sprintf(figFilenamePattern, estResNumber, r, "png")
        htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, r, "html")
        fig <- getPlotTrueInitialAndEstimatedMatrices(true=simRes$V0[r,,], initial=initialConds$V0[r,,], estimated=estRes$V0[r,,], title=sprintf("V0 for trial %03d", r))
        htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
        # orca(p=fig, file=pngFigFilename)
        # print(fig)
    }

    browser()
}

processAll()

