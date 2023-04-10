
require(plotly)
require(htmlwidgets)

processAll <- function(totalDurationSamples=5000, baselineValue=0, nSteps=50, stepValue=1, stepLength=5, sFilenamePattern="data/s_multipleStim_steps_durationSamples%d_baselineValue%.02f_nSteps%d_stepValue%.02f_stepLength%d.csv", figFilenamePattern="figures/s_multipleStim_steps_durationSamples%d_baselineValue%.02f_nSteps%d_stepValue%.02f_stepLength%d.%s") {
    sFilename <- sprintf(sFilenamePattern, totalDurationSamples, baselineValue, nSteps, stepValue, stepLength)
    pngFilename <- sprintf(figFilenamePattern, totalDurationSamples, baselineValue, nSteps, stepValue, stepLength, "png")
    htmlFilename <- sprintf(figFilenamePattern, totalDurationSamples, baselineValue, nSteps, stepValue, stepLength, "html")
    step_start_candidate = seq(from=1, to=totalDurationSamples-stepLength, by=stepLength)
    step_start <- sample(x=step_start_candidate, size=nSteps)
    sSample <- matrix(rep.int(0, times=2*totalDurationSamples), nrow=2)

    for(s in step_start) {
        signSample <- sample(c(-1, 1), 1)
        if(signSample>0) {
            sSample[1, s+(0:(stepLength-1))] <- stepValue
        } else {
            sSample[2, s+(0:(stepLength-1))] <- stepValue
        }
    }
    write.table(sSample, sFilename, row.names=FALSE, col.names=FALSE)

    fig <- plot_ly(type='scatter', mode='markers+lines')
    fig <- fig%>%add_trace(x=1:totalDurationSamples, y=sSample[1,], name="stim 1")
    fig <- fig%>%add_trace(x=1:totalDurationSamples, y=sSample[2,], name="stim 2")
    fig <- fig %>% layout(xaxis=list(title="Sample Number"), yaxis=list(title="Stimulus Value"))
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    # print(fig)

    browser()
}

processAll()

rm(processAll)
