library(ggplot2)
library(mvtnorm)
library(plotly)

source("../src/smoothLDS.R")
source("getNormalEllipse.R")

processAll <- function() {
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
