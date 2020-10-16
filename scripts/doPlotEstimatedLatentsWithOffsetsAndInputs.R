
require(plotly)
require(RColorBrewer)
require(ini)
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")
source("../src/filterLDS_SS_withOffsetsAndInputs.R")
source("../src/smoothLDS_SS_withOffsetsAndInputs.R")

processAll <- function() {
    # estResNumber <- 99485472
    # estResNumber <- 67957061
    # estResNumber <- 93973929
    # estResNumber <- 59934313
    estResNumber <- 43378880
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
    dim(simRes$c) <- c(1, 1, length(simRes$c))
    dim(simRes$d) <- c(1, 1, length(simRes$d))

    # start run KF and KS using with initial parameters
    # filterLDS_SS_withOffsetsAndInputs <- function(y, B, u, C, c, Q, m0, V0, Z, a, D, d, R, initStateAt=0) {
    fRes0 <- filterLDS_SS_withOffsetsAndInputs(y=simRes$y, B=initialConds$B, u=initialConds$u0, C=initialConds$C0, c=simRes$c, Q=initialConds$Q, m0=initialConds$m0, V0=initialConds$V0, Z=initialConds$Z, a=initialConds$a0, D=initialConds$D, d=simRes$d, R=initialConds$R)
    # smoothLDS_SS_withOffsetsAndInputs <- function(B, u, C, c, Q, xnn, Vnn, xnn1, Vnn1, initStateAt=0, m0=NA, V0=NA) {
    sRes0 <- smoothLDS_SS_withOffsetsAndInputs(B=initialConds$B, u=initialConds$u0, C=initialConds$C0, c=simRes$c, Q=initialConds$Q, xnn=fRes0$xnn, Vnn=fRes0$Vnn, xnn1=fRes0$xnn1, Vnn1=fRes0$Vnn1, m0=simRes$x0, V0=simRes$V0)
    # end run KF and KS using with initial conditions

    # start run KF and KS using with estimated parameters
    fRes <- filterLDS_SS_withOffsetsAndInputs(y=simRes$y, B=estRes$B, u=estRes$u, C=estRes$C, c=simRes$c, Q=estRes$Q, m0=estRes$m0, V0=estRes$V0, Z=estRes$Z, a=estRes$a, D=estRes$D, d=simRes$d, R=estRes$R)
    sRes <- smoothLDS_SS_withOffsetsAndInputs(B=estRes$B, u=estRes$u, C=estRes$C, c=simRes$c, Q=estRes$Q, xnn=fRes$xnn, Vnn=fRes$Vnn, xnn1=fRes$xnn1, Vnn1=fRes$Vnn1, m0=estRes$m0, V0=estRes$V0)
    # end run KF and KS using with estimated parameters

    trueLinetype <- "solid"
    initialLinetype <- "dot"
    estLinetype <- "dash"
    cols <- brewer.pal(max(3, nrow(simRes$x)), "Set1")
    fig <- plot_ly(type='scatter', mode='lines')
    for(i in 1:nrow(simRes$x)) {
        fig <- fig%>%add_trace(x=1:length(simRes$x[i,]), y=simRes$x[i,], name=sprintf("true[,%d]", i), line=list(color=cols[i], dash=trueLinetype))
    }
    for(i in 1:nrow(simRes$x)) {
        fig <- fig%>%add_trace(x=1:length(sRes0$xnN[i,1,]), y=sRes0$xnN[i,1,], name=sprintf("initial[,%d]", i), line=list(color=cols[i], dash=initialLinetype))
    }
    for(i in 1:nrow(simRes$x)) {
        fig <- fig%>%add_trace(x=1:length(sRes$xnN[i,1,]), y=sRes$xnN[i,1,], name=sprintf("estimated[,%d]", i), line=list(color=cols[i], dash=estLinetype))
    }
    fig <- fig %>% layout(xaxis = list(title="Sample"), yaxis = list (title="Latent Value"))
    pngFigFilename <- sprintf("figures//%.08d_latents.%s", estResNumber, "png")
    htmlFigFilename <- sprintf("figures//%.08d_latents.%s", estResNumber, "html")
    orca(p=fig, file=pngFigFilename)
    htmlwidgets::saveWidget(as_widget(fig), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))
    # print(fig)

    browser()
}

processAll()

