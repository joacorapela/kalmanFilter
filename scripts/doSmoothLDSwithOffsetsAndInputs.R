library(ggplot2)
library(mvtnorm)
library(plotly)

source("../src/smoothLDS_SS_withOffsetsAndInputs.R")
source("getNormalEllipse.R")

processAll <- function() {
    filterResNumber <- 45380418
    simFilenamePattern <- "results/%08d_simulation.RData"
    filterResMetaDataFilenamePattern <- "results/%08d_filteredSimulation.ini"
    filterResFilenamePattern <- "results/%08d_filteredSimulation.RData"
    smoothResFilenamePattern <- "results/%08d_smoothedSimulation.RData"
    smoothResMetaDataFilenamePattern <- "results/%08d_smoothedSimulation.ini"
    smoothResFigFilenamePattern <- "figures/%08d_smoothedSimulation.%s"

    exit <- FALSE
    while(!exit) {
        smoothResNumber <- sample(1e8, 1)
        smoothResFilename <- sprintf(smoothResFilenamePattern, smoothResNumber)
        if(!file.exists(smoothResFilename)) {
            exit <- TRUE
        }
    }
    smoothResMetaDataFilename <- sprintf(smoothResMetaDataFilenamePattern, smoothResNumber)
    show(sprintf("Saving estimation results in: %s", smoothResFilename))
    show(sprintf("Saving estimation meta data in: %s", smoothResMetaDataFilename))

    filterResMetaDataFilename <- sprintf(filterResMetaDataFilenamePattern, filterResNumber)
    ini <- read.ini(filepath=filterResMetaDataFilename)
    simResNumber <- as.integer(ini$simulation_info$simResNumber)
    simFilename <- sprintf(simFilenamePattern, simResNumber)
    simRes <- get(load(simFilename))

    filterResFilename <- sprintf(filterResFilenamePattern, filterResNumber)
    filterRes <- get(load(filterResFilename))

    smoothRes <- smoothLDS_SS_withOffsetsAndInputs(B=simRes$B, u=simRes$u, C=simRes$C, c=simRes$c, Q=simRes$Q, xnn=filterRes$xnn, Vnn=filterRes$Vnn, xnn1=filterRes$xnn1, Vnn1=filterRes$Vnn1, stateType0="init00", x00=simRes$x00, V00=simRes$V00, c0=0.0)
    save(smoothRes, file=smoothResFilename)

    metaData <- list()
    metaData[["filter_info"]] <- list(filterResNumber=filterResNumber)
    write.ini(x=metaData, filepath=smoothResMetaDataFilename)

    nObs <- ncol(simRes$c)
    df <- data.frame(x=c(simRes$x[1,], simRes$y[1,], filterRes$xnn[1,1,], smoothRes$xnN[1,1,]), y=c(simRes$x[2,], simRes$y[2,], filterRes$xnn[2,1,], smoothRes$xnN[2,1,]), type=factor(c(rep("latent", nObs), rep("measurement", nObs), rep("filtered", nObs), rep("smoothed", nObs))))
    hoverTextLatents <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, simRes$x[1,], simRes$x[2,])
    hoverTextObservations <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, simRes$y[1,], simRes$y[2,])
    hoverTextFiltered <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, filterRes$xnn[1,1,], filterRes$xnn[2,1,])
    hoverTextSmoothed <- sprintf("sample %d, x %.02f, y %.02f", 1:nObs, smoothRes$xnN[1,1,], smoothRes$xnN[2,1,])
    df <- cbind(df, c(hoverTextLatents, hoverTextObservations, hoverTextFiltered, hoverTextSmoothed))
    colnames(df) <- c("x", "y", "type", "hoverText")
    fig <- plot_ly(data=df, type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=~x, y=~y, text=~hoverText, color=~type, hoverinfo="text")
    smoothResPNGFilename <- sprintf(smoothResFigFilenamePattern, smoothResNumber, "png")
    smoothResHTMLFilename <- sprintf(smoothResFigFilenamePattern, smoothResNumber, "html")
    orca(p=fig, file=smoothResPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(smoothResHTMLFilename)),basename(smoothResHTMLFilename)))
    print(fig)

    browser()
    dObs <- 2*4
    ellipseNPoints <- 100
    ellipseCriticalValue <- .95
    ellipseCol <- "red"
    propByEllipsePlot <- .10
    ellipsePointSize <- .1
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulationCircle.RData"
    filterResFilename <- "results/srFilterResSimulationCircle.RData"
    smoothResFilename <- "results/srSmootherResSimulationCircle.RData"
    smoothResFigFilename <- "figures/srSmootherResSimulationCircle.png"

    simRes <- get(load(simulationFilename))
    filterRes <- get(load(filterResFilename))
    x <- simRes$x
    z <- simRes$z
    A <- simRes$A
    P <- filterRes$SigmaX[2:length(filterRes$SigmaX)]
    mu <- filterRes$xHat
    V <- filterRes$SigmaXHat
    nObs <- ncol(x)
    byEllipsePlot <- propByEllipsePlot*nObs

    smoothRes <- smoothLDS(A=A, mu=mu, V=V, P=P)
    save(smoothRes, file=smoothResFilename)

    df <- data.frame(x=c(simRes$x[1,], simRes$z[1,], smoothRes$muHat[1,]), 
                     y=c(simRes$x[dObs,], simRes$z[2,], smoothRes$muHat[2,]),
                     type=factor(c(rep("measurement", nObs), rep("latent", nObs), rep("smoothed", nObs))))
    p <- ggplot()
    p <- p + geom_point(data=df, aes(x=x, y=y, colour=type))
    p <- p + scale_colour_manual(values=c("measurement"="blue", 
                                          "latent"="green", 
                                          "smoothed"="red"))
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)
    p <- p + theme(legend.title = element_blank()) 

    for(i in seq(from=1, to=nObs, by=byEllipsePlot)) {
        ellipse <- getNormalEllipse(mu=smoothRes$muHat[,i], 
                                    covar=smoothRes$VHat[[i]], 
                                    nPoints=ellipseNPoints, 
                                    criticalValue=ellipseCriticalValue)
        ellipseDF <- data.frame(x=ellipse[,1], y=ellipse[,2])
        p <- p + geom_point(data=ellipseDF, aes(x=x, y=y), colour=ellipseCol, 
                            size=ellipsePointSize)
    }

    ggsave(filename=smoothResFigFilename, plot=p)

    p <- ggplotly(p)
    print(p)

    browser()
}

processAll()
