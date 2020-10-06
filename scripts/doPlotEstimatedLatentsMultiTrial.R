
require(plotly)
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")
source("../src/filterLDS_SS.R")
source("../src/smoothLDS_SS.R")

processAll <- function() {
    r <- 1
    # estResNumber <- 99485472
    # estResNumber <- 67957061
    # estResNumber <- 93973929
    estResNumber <- 66860571
    simFilenamePattern <- "results/%08d_simulation.RData"
    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"

    estResFilename <- sprintf(estResFilenamePattern, estResNumber)
    estRes <- get(load(estResFilename))
    # initialConds <- estRes$initialConds

    estResMetaDataFilename <- sprintf(estResMetaDataFilenamePattern, estResNumber)
    ini <- read.ini(filepath=estResMetaDataFilename)
    simResNumber <- as.integer(ini$simulation_info$simResNumber)
    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(simFilename))

    N <- dim(simRes$y)[3]
    P <- dim(simRes$y)[2]
    M <- dim(simRes$x)[2]

    # start run KF and KS using with initial parameters
    # fRes0 <- filterLDS_SS(y=simRes$y[r,,], B=initialConds$B, Z=initialConds$Z, m0=initialConds$m0[r,], V0=initialConds$V0[r,,], Q=initialConds$Q, R=initialConds$R)
    # sRes0 <- smoothLDS_SS(B=initialConds$B, xnn=fRes0$xnn, Vnn=fRes0$Vnn, xnn1=fRes0$xnn1, Vnn1=fRes0$Vnn1, m0=simRes$x0[r,], V0=simRes$V0[r,,])
    # end run KF and KS using with initial conditions

    # start run KF and KS using with estimated parameters
    fRes <- filterLDS_SS(y=simRes$y[r,,], B=estRes$B, Z=estRes$Z, m0=estRes$m0[r,], V0=estRes$V0[r,,], Q=estRes$Q, R=estRes$R)
    sRes <- smoothLDS_SS(B=estRes$B, xnn=fRes$xnn, Vnn=fRes$Vnn, xnn1=fRes$xnn1, Vnn1=fRes$Vnn1, m0=estRes$x0[r,], V0=estRes$V0[r,,])
    # end run KF and KS using with estimated parameters

    trueLinetype <- "solid"
    # initialLinetype <- "dot"
    estLinetype <- "dash"
    cols <- brewer.pal(max(3, M), "Set1")
    fig <- plot_ly(type='scatter', mode='lines')
    # simRes$x \in nTrials x nLatents x nSamples
    for(m in 1:M) {
        rgbValues <- col2rgb(cols[m])
        fig <- fig%>%add_trace(x=1:N, y=simRes$x[r,m,], name=sprintf("true[,%d]", m), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash=trueLinetype))
    }
    # xnN \in nLatents x 1 x nSamples
    # VnN \in nLatents x nLatents x nSamples
    eMeans <- matrix(data=NA, nrow=M,ncol=N)
    eSTDs <- matrix(data=NA, nrow=M,ncol=N)
    for(n in 1:N) {
        eMeans[,n] <- sRes$xnN[,1,n]
        eSTDs[,n] <- sqrt(diag(sRes$VnN[,,n]))
    }
    # browser()
    for(m in 1:M) {
        rgbValues <- col2rgb(cols[m])
        fig <- fig%>%add_trace(x=1:N, y=eMeans[m,], mode="lines", name=sprintf("estimated[,%d]", m), line=list(color=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 1), dash=estLinetype))
        fig <- fig%>%add_trace(x=1:N, y=eMeans[m,]-1.96*eSTDs[m,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE)
        fig <- fig%>%add_trace(x=1:N, y=eMeans[m,]+1.96*eSTDs[m,], mode="lines", line=list(color="rgba(0,0,0,0)"), showlegend=FALSE, fill="tonexty", fillcolor=sprintf("rgba(%d,%d,%d,%f)", rgbValues[1,1], rgbValues[2,1], rgbValues[3,1], 0.2))
    }
    fig <- fig %>% layout(title=sprintf("Trial %03d", r), xaxis = list(title="Sample"), yaxis = list (title="Value"))
    figFilenamePattern <- "figures//%.08d_trial%03d_latents.%s"
    # pngFigFilename <- sprintf(figFilenamePattern, estResNumber, r, "png")
    htmlFigFilename <- sprintf(figFilenamePattern, estResNumber, r, "html")
    # orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # print(fig)

    browser()
}

processAll()

