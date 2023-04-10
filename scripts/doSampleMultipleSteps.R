
require(plotly)
require(htmlwidgets)

processAll <- function(totalDurationSamples=5000, baselineValue=0, nSteps=50, stepValue=1, stepLength=5, cFilenamePattern="data/s_multiplSignedSteps_durationSamples%d_baselineValue%.02f_nSteps%d_stepValue%.02f_stepLength%d.csv", figFilenamePattern="figures/s_multiplSignedSteps_durationSamples%d_baselineValue%.02f_nSteps%d_stepValue%.02f_stepLength%d.%s") {
    cFilename <- sprintf(cFilenamePattern, totalDurationSamples, baselineValue, nSteps, stepValue, stepLength)
    pngFilename <- sprintf(figFilenamePattern, totalDurationSamples, baselineValue, nSteps, stepValue, stepLength, "png")
    htmlFilename <- sprintf(figFilenamePattern, totalDurationSamples, baselineValue, nSteps, stepValue, stepLength, "html")
    step_start_candidate = seq(from=1, to=totalDurationSamples-stepLength, by=stepLength)
    step_start <- sample(x=step_start_candidate, size=nSteps)
    cSample <- rep.int(0, times=totalDurationSamples)

    for(s in step_start) {
        signSample <- sample(c(-1, 1), 1)
        cSample[s+(0:(stepLength-1))] <- signSample*stepValue
    }
    write.table(cSample, cFilename, row.names=FALSE, col.names=FALSE)

    fig <- plot_ly(type='scatter', mode='markers+lines')
    fig <- fig%>%add_trace(x=1:length(cSample), y=cSample)
    fig <- fig %>% layout(xaxis=list(title="Sample Number"), yaxis=list(title="Stimulus Value"))
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # print(fig)

    browser()
}

processAll()

rm(processAll)
