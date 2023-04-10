
require(plotly)
require(htmlwidgets)

processAll <- function() {
    # simResNumber <- 68668073
    simResNumber <- 82902398
    # simResNumber <- 94887612
    # simResNumber <- 06798128
    # simResNumber <- 49312087
    # simResNumber <- 56876727
    # simResNumber <- 45051359
    # simResNumber <- 16514070
    # simResNumber <- 68430956
    # simResNumber <- 23386188
    # simResNumber <- 39298476
    # simResNumber <- 46195576
    # simResNumber <- 78973043
    # simResNumber <- 01879479
    simFilenamePattern <- "results/%08d_simulation.RData"
    latentsFigFilenamePattern <- "figures/%08d_simulation_latents.%s"
    observationsFigFilenamePattern <- "figures/%08d_simulation_observations.%s"

    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(file=simFilename))

    stim1Diffs <- diff(simRes$c[1,])
    stim1Onsets <- which(stim1Diffs>0)+1
    stim1Offsets <- which(stim1Diffs<0)

    stim2Diffs <- diff(simRes$c[2,])
    stim2Onsets <- which(stim2Diffs>0)+1
    stim2Offsets <- which(stim2Diffs<0)

    nObs <- dim(simRes$x)[2]
    obsIndices <- 1:nObs
    hoverText <- sprintf("sample %d, x1 %.02f, x2 %.02f", obsIndices, simRes$x[1,], simRes$x[2,])
    df <- data.frame(x1=simRes$x[1,], x2=simRes$x[2,])
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x1, y=~x2, text=~hoverText, hoverinfo="text")
    fig <- fig %>% add_annotations(x=c(simRes$x[1,stim1Onsets], simRes$x[1, stim1Offsets]), y=c(simRes$x[2, stim1Onsets], simRes$x[2, stim1Offsets]), text=c(rep("stim1 start", times=length(stim1Onsets)), rep("stim1 end", times=length(stim1Offsets))))                        
    fig <- fig %>% add_annotations(x=c(simRes$x[1,stim2Onsets], simRes$x[1, stim2Offsets]), y=c(simRes$x[2, stim2Onsets], simRes$x[2, stim2Offsets]), text=c(rep("stim2 start", times=length(stim2Onsets)), rep("stim2 end", times=length(stim2Offsets))))                        
    fig <- fig %>% layout(xaxis = list(title="Latent 1"), yaxis = list(title="Latent 2"))
    latentsPNGFilename <- sprintf(latentsFigFilenamePattern, simResNumber, "png")
    latentsHTMLFilename <- sprintf(latentsFigFilenamePattern, simResNumber, "html")
    # orca(p=fig, file=simPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(latentsHTMLFilename)), basename(latentsHTMLFilename)))
    # print(fig)

    # neuronsToPlot <- c(1,3,14,20)
    neuronsToPlot <- 1:7
    fig <- plot_ly(type="scatter", mode="lines+markers")
    for(neuronToPlot in neuronsToPlot) {
        fig <- fig %>% add_trace(x=obsIndices, y=simRes$z[neuronToPlot,], name=sprintf("z%d", neuronToPlot))
        fig <- fig %>% add_trace(x=obsIndices, y=simRes$y[neuronToPlot,], name=sprintf("y%d", neuronToPlot))
    }
    fig <- fig %>% add_trace(x=obsIndices, y=simRes$c[1,], name="stim1")
    fig <- fig %>% add_trace(x=obsIndices, y=simRes$c[2,], name="stim1")
    fig <- fig %>% layout(xaxis=list(title="Sample Number"), yaxis=list(title=""))
    observationsPNGFilename <- sprintf(observationsFigFilenamePattern, simResNumber, "png")
    observationsHTMLFilename <- sprintf(observationsFigFilenamePattern, simResNumber, "html")
    # orca(p=fig, file=simPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(observationsHTMLFilename)), basename(observationsHTMLFilename)))
    # print(fig)
    browser()
}

processAll()
rm(processAll)

