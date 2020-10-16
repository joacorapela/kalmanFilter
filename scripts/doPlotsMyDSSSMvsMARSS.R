
require(plotly)
require(RColorBrewer)
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")
source("../src/filterLDS_SS_withOffsetsAndInputs.R")
source("../src/smoothLDS_SS_withOffsetsAndInputs.R")

processAll <- function() {
    myEMestResNumber <- 48200862
    MARSSestResNumber <- 18421034
    simFilenamePattern <- "results/%08d_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"

    myEMestResFilename <- sprintf(estResFilenamePattern, myEMestResNumber)
    MARSSestResFilename <- sprintf(estResFilenamePattern, MARSSestResNumber)
    myEMestRes <- get(load(myEMestResFilename))
    MARSSestRes <- get(load(MARSSestResFilename))

    fig <- plot_ly(type='scatter', mode='markers')
    fig <- fig%>%add_trace(x=1:length(myEMestRes$logLike), y=myEMestRes$logLike, name="myEM")
    fig <- fig%>%add_trace(x=1:length(MARSSestRes$iter.record$logLik), y=MARSSestRes$iter.record$logLik, name="MARSS")
    fig <- fig %>% layout(xaxis=list(title="Iteration Number"), yaxis=list(title="Log Likelihood"))
    pngFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_logLikeVsIterNo.%s", myEMestResNumber, MARSSestResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_logLikeVsIterNo.%s", myEMestResNumber, MARSSestResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(fig)

    fig <- plot_ly(type='bar',
                   x=c("MyEM", "MARSS"),
                   y=c(myEMestRes$elapsedTime, 
                       MARSSestRes$elapsedTime))
    fig <- fig %>% layout(yaxis=list(title="Elapsed Time"))
    pngFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_elapsedTime.%s", myEMestResNumber, MARSSestResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_%08d_myEMvsMARSS_elapsedTime.%s", myEMestResNumber, MARSSestResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(fig)

    browser()
}

processAll()

