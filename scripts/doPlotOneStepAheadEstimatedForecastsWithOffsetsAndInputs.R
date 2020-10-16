
require(plotly)
require(ini)
require(RColorBrewer)
source("../src/filterLDS_SS_withOffsetsAndInputs.R")
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    obsToPlot <- 1
    # testSimResNumber <- 09685977
    testSimResNumber <- -1 # use train simulation
    # estResNumber <- 99485472
    # estResNumber <- 93973929
    # estResNumber <- 59934313
    estResNumber <- 27646949
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
    dim(testSimRes$c) <- c(1, 1, length(testSimRes$c))
    dim(testSimRes$d) <- c(1, 1, length(testSimRes$d))

    N <- ncol(testSimRes$y)
    P <- nrow(testSimRes$y)
    M <- nrow(testSimRes$x)

    eFRes <- filterLDS_SS_withOffsetsAndInputs(y=testSimRes$y, B=estRes$B, u=estRes$u, C=estRes$C, c=testSimRes$c, Q=estRes$Q, m0=estRes$m0, V0=estRes$V0, Z=estRes$Z, a=estRes$a, D=estRes$D, d=testSimRes$d, R=estRes$R)
    eMeanOneStepAheadForecast <- estRes$Z%*%eFRes$xnn1[,1,]+estRes$u[,1]%o%rep(1, times=N)+estRes$C%*%matrix(data=testSimRes$c, nrow=1)
    eSTDOneStepAheadForecast <- matrix(NA, nrow=P, ncol=N)
    for(i in 1:N) {
        eSTDOneStepAheadForecast[,i-1] <- sqrt(diag(estRes$Z%*%eFRes$Vnn1[,,i]%*%t(estRes$Z)+estRes$R))
    }

    fig <- plot_ly(type='scatter', mode="markers")
    cols <- brewer.pal(max(3, P), "Set1")
    for(i in 1:P) {
        rgbValues <- col2rgb(cols[i])
        # observation
        fig <- fig%>%add_trace(x=1:N, y=testSimRes$y[i,1:N], mode="markers", name=sprintf("observation[,%d]", i), marker=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.5)))
        # forecast
        fig <- fig%>%add_trace(x=1:N, y=eMeanOneStepAheadForecast[i,], mode="lines", name=sprintf("forecast[,%d]", i), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash="solid", width=0.5))
        fig <- fig%>%add_trace(x=1:N, y=eMeanOneStepAheadForecast[i,]+1.96*eSTDOneStepAheadForecast[i,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE)
        fig <- fig%>%add_trace(x=1:N, y=eMeanOneStepAheadForecast[i,]-1.96*eSTDOneStepAheadForecast[i,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE, fill="tonexty", fillcolor=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.2))
    }
    fig <- fig %>% layout(title=sprintf("Log Likelihood=%f", eFRes$logLike), xaxis=list(title="Sample"), yaxis=list (title="Observation/Forecast"))
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, testSimResNumber, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, testSimResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    print(fig)

    browser()
}

processAll()

