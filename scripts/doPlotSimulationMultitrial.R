
require(plotly)
require(htmlwigets)

processAll <- function() {
    simResNumber <- 78598469
    trialToPlot <- 8
    simFilenamePattern <- "results/%08d_simulation.RData"
    simFigFilenamePattern <- "figures/%08d_simulation_trial%02d.%s"

    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(file=simFilename))

    nObs <- dim(simRes$x)[3]
    title <- sprintf("Trial %d", trialToPlot)
    hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, simRes$x[trialToPlot,1,], simRes$x[trialToPlot,2,])
    hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, simRes$y[trialToPlot,1,], simRes$y[trialToPlot,2,])
    df <- data.frame(t(cbind(simRes$x[trialToPlot,,], simRes$y[trialToPlot,,])))
    df <- cbind(df, c(rep("latent", nObs), rep("measurement", nObs)))
    df <- cbind(df, c(hoverTextLatents, hoverTextObservations))
    colnames(df) <- c("x", "y", "type", "hoverText")
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
    fig <- fig %>% add_annotations(x=c(simRes$x[trialToPlot,1,1], simRes$x[trialToPlot,1,nObs]), y=c(simRes$x[trialToPlot,2,1], simRes$x[trialToPlot,2,nObs]), text=c("start", "end"))
    fig <- fig %>% layout(title=title)
    simPNGFilename <- sprintf(simFigFilenamePattern, simResNumber, trialToPlot, "png")
    simHTMLFilename <- sprintf(simFigFilenamePattern, simResNumber, trialToPlot, "html")
    # orca(p=fig, file=simPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(simHTMLFilename)),basename(simHTMLFilename)))
    print(fig)

}

processAll()
rm(processAll)

