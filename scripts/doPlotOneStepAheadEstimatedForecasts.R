
require(plotly)
require(ini)
require(RColorBrewer)
source("../src/filterLDS_SS.R")
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    obsToPlot <- 1
    testSimResNumber <- -1 # use train simulation
    estResNumber <- 99485472
    simFilenamePattern <- "results/%08d_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"
    figFilenamePattern <- "figures//%.08d_%.08d_estimatedForecasts.%s"

    estResFilename <- sprintf(estResFilenamePattern, estResNumber)
    estRes <- get(load(estResFilename))
    initialConds <- estRes$initialConds

    estResMetaDataFilename <- sprintf(estResMetaDataFilenamePattern, estResNumber)
    ini <- read.ini(filepath=estResMetaDataFilename)
    simResNumber <- as.integer(ini$simulation_info$simResNumber)
    if(testSimResNumber<0) {
        testSimResNumber <- simResNumber
    }
    testSimFilename <- sprintf(simFilenamePattern, testSimResNumber)
    testSimRes <- get(load(testSimFilename))

    N <- ncol(testSimRes$y)
    P <- nrow(testSimRes$y)
    M <- nrow(testSimRes$x)

    eFRes <- filterLDS_SS(y=testSimRes$y, B=estRes$B, Z=estRes$Z, m0=estRes$m0, V0=estRes$V0, Q=estRes$Q, R=estRes$R)
    eMeanOneStepAheadForecast <- estRes$Z%*%eFRes$xnn1[,1,]
    eSTDOneStepAheadForecast <- matrix(NA, nrow=P, ncol=N)
    for(i in 1:N) {
        eSTDOneStepAheadForecast[,i] <- sqrt(diag(estRes$Z%*%eFRes$Vnn1[,,i]%*%t(estRes$Z)+estRes$R))
    }

    fig <- plot_ly(type='scatter', mode="markers")
    cols <- brewer.pal(max(3, P), "Set1")
    for(i in 1:P) {
        rgbValues <- col2rgb(cols[i])
        # observation
        fig <- fig%>%add_trace(x=1:ncol(testSimRes$y), y=testSimRes$y[i,], mode="markers", name=sprintf("observation[,%d]", i), marker=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.5)))
        # estimate
        fig <- fig%>%add_trace(x=1:ncol(eMeanOneStepAheadForecast), y=eMeanOneStepAheadForecast[i,], mode="lines", name=sprintf("estimate[,%d]", i), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash="solid"))
        fig <- fig%>%add_trace(x=1:ncol(eMeanOneStepAheadForecast), y=eMeanOneStepAheadForecast[i,]+1.96*eSTDOneStepAheadForecast[i,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE)
        fig <- fig%>%add_trace(x=1:ncol(eMeanOneStepAheadForecast), y=eMeanOneStepAheadForecast[i,]-1.96*eSTDOneStepAheadForecast[i,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE, fill="tonexty", fillcolor=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.2))
    }
    fig <- fig %>% layout(title=sprintf("Log Likelihood=%f", eFRes$logLike), xaxis=list(title="Sample"), yaxis=list (title="Observation"))
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, testSimResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, testSimResNumber, "html")
    # orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # print(fig)

    browser()
}

processAll()

