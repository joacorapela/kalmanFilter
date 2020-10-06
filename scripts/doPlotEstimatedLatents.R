
require(plotly)
source("../src/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/getPlotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    # estResNumber <- 99485472
    # estResNumber <- 67957061
    estResNumber <- 93973929
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

    # start run KF and KS using with initial parameters
    fRes0 <- filterLDS_SS(y=simRes$y, B=initialConds$B, Z=initialConds$Z, m0=initialConds$m0, V0=initialConds$V0, Q=initialConds$Q, R=initialConds$R)
    sRes0 <- smoothLDS_SS(B=initialConds$B, xnn=fRes0$xnn, Vnn=fRes0$Vnn, xnn1=fRes0$xnn1, Vnn1=fRes0$Vnn1, m0=simRes$x0, V0=simRes$V0)
    # end run KF and KS using with initial conditions

    # start run KF and KS using with estimated parameters
    fRes <- filterLDS_SS(y=simRes$y, B=estRes$B, Z=estRes$Z, m0=estRes$m0, V0=estRes$V0, Q=estRes$Q, R=estRes$R)
    sRes <- smoothLDS_SS(B=estRes$B, xnn=fRes$xnn, Vnn=fRes$Vnn, xnn1=fRes$xnn1, Vnn1=fRes$Vnn1, m0=estRes$x0, V0=estRes$V0)
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
    print(fig)

    browser()
}

processAll()

