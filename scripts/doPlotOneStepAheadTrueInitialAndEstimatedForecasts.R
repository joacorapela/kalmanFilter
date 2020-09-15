
require(plotly)
require(ini)
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

    # true
    tFRes <- filterLDS_SS(y=testSimRes$y, B=testSimRes$B, Z=testSimRes$Z, m0=testSimRes$m0, V0=testSimRes$V0, Q=testSimRes$Q, R=testSimRes$R)
    tMeanOneStepAheadForecast <- testSimRes$Z%*%tFRes$xnn1[,1,]
    tSTDOneStepAheadForecast <- matrix(NA, nrow=P, ncol=N)
    for(i in 1:N) {
        tSTDOneStepAheadForecast[,i] <- sqrt(diag(testSimRes$Z%*%tFRes$Vnn1[,,i]%*%t(testSimRes$Z)+testSimRes$R))
    }

    # initial
    iFRes <- filterLDS_SS(y=testSimRes$y, B=initialConds$B, Z=initialConds$Z, m0=initialConds$m0, V0=initialConds$V0, Q=initialConds$Q, R=initialConds$R)
    iMeanOneStepAheadForecast <- initialConds$Z%*%iFRes$xnn1[,1,]
    iSTDOneStepAheadForecast <- matrix(NA, nrow=P, ncol=N)
    for(i in 1:N) {
        iSTDOneStepAheadForecast[,i] <- sqrt(diag(initialConds$Z%*%iFRes$Vnn1[,,i]%*%t(initialConds$Z)+initialConds$R))
    }

    # estimated
    eFRes <- filterLDS_SS(y=testSimRes$y, B=estRes$B, Z=estRes$Z, m0=estRes$m0, V0=estRes$V0, Q=estRes$Q, R=estRes$R)
    eMeanOneStepAheadForecast <- estRes$Z%*%eFRes$xnn1[,1,]
    eSTDOneStepAheadForecast <- matrix(NA, nrow=P, ncol=N)
    for(i in 1:N) {
        eSTDOneStepAheadForecast[,i] <- sqrt(diag(estRes$Z%*%eFRes$Vnn1[,,i]%*%t(estRes$Z)+estRes$R))
    }

    fig <- plot_ly(type='scatter', mode="markers")
    # observation
    fig <- fig%>%add_trace(x=1:ncol(testSimRes$y), y=testSimRes$y[obsToPlot,], mode="markers", name=sprintf("observation[,%d]", obsToPlot), marker=list(color="black"))

    # true
    fig <- fig%>%add_trace(x=1:ncol(tMeanOneStepAheadForecast), y=tMeanOneStepAheadForecast[obsToPlot,], mode="lines", name=sprintf("true[,%d]", obsToPlot), line=list(color="rgba(0,255,0,1)", dash="solid"))
    fig <- fig%>%add_trace(x=1:ncol(tMeanOneStepAheadForecast), y=tMeanOneStepAheadForecast[obsToPlot,]+1.96*tSTDOneStepAheadForecast[obsToPlot,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE)
    fig <- fig%>%add_trace(x=1:ncol(tMeanOneStepAheadForecast), y=tMeanOneStepAheadForecast[obsToPlot,]-1.96*tSTDOneStepAheadForecast[obsToPlot,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE, fill="tonexty", fillcolor="rgba(0,255,0,0.2")

    # initial
    fig <- fig%>%add_trace(x=1:ncol(iMeanOneStepAheadForecast), y=iMeanOneStepAheadForecast[obsToPlot,], mode="lines", name=sprintf("initial[,%d]", obsToPlot), line=list(color="rgba(0,0,255,1)", dash="solid"))
    fig <- fig%>%add_trace(x=1:ncol(iMeanOneStepAheadForecast), y=iMeanOneStepAheadForecast[obsToPlot,]+1.96*iSTDOneStepAheadForecast[obsToPlot,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE)
    fig <- fig%>%add_trace(x=1:ncol(iMeanOneStepAheadForecast), y=iMeanOneStepAheadForecast[obsToPlot,]-1.96*iSTDOneStepAheadForecast[obsToPlot,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE, fill="tonexty", fillcolor="rgba(0,0,255,0.2")

    # estimate
    fig <- fig%>%add_trace(x=1:ncol(eMeanOneStepAheadForecast), y=eMeanOneStepAheadForecast[obsToPlot,], mode="lines", name=sprintf("estimated[,%d]", obsToPlot), line=list(color="rgba(255,0,0,1)", dash="solid"))
    fig <- fig%>%add_trace(x=1:ncol(eMeanOneStepAheadForecast), y=eMeanOneStepAheadForecast[obsToPlot,]+1.96*eSTDOneStepAheadForecast[obsToPlot,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE)
    fig <- fig%>%add_trace(x=1:ncol(eMeanOneStepAheadForecast), y=eMeanOneStepAheadForecast[obsToPlot,]-1.96*eSTDOneStepAheadForecast[obsToPlot,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE, fill="tonexty", fillcolor="rgba(255,0,0,0.2")

    fig <- fig %>% layout(xaxis=list(title="Sample"), yaxis=list (title="Observation"))
    figFilenamePattern <- "figures//%.08d_%.08d_obs%dTrueInitialEstimatedForecasts.%s"
    pngFigFilename <- sprintf(figFilenamePattern, estResNumber, testSimResNumber, obsToPlot, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, testSimResNumber, obsToPlot, "html")
    # orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # print(fig)

    browser()
}

processAll()

