require(ggplot2)
require(mvtnorm)
require(ramcmc)                                                       

source("../src/smoothLDS.R")
source("getNormalEllipse.R")

processAll <- function() {
    ellipseNPoints <- 100
    ellipseCriticalValue <- .95
    ellipseCol <- "red"
    propByEllipsePlot <- .10
    ellipsePointSize <- .1
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulationPendulum.RData"
    filterResFilename <- "results/srFilterResPendulum.RData"
    smoothResFilename <- "results/srSmootherResPendulum.RData"
    smoothResFigFilename <- "figures/srSmootherResPendulum.png"

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

    df <- data.frame(x=c(simRes$t, simRes$t, simRes$t, simRes$t, simRes$t), y=c(simRes$x[1,], simRes$z[1,], simRes$z[2, ], smoothRes$muHat[1,], smoothRes$muHat[2,]), type=factor(c(rep("measurement", nObs), rep("latent acceleration", nObs), rep("latent position", nObs), rep("smoothed acceleration", nObs), rep("smoothed position", nObs))))
    p <- ggplot()
    p <- p + geom_point(data=df, aes(x=x, y=y, colour=type))
    p <- p + scale_colour_manual(values=c("measurement"="blue", "latent acceleration"="lightgreen", "latent position"="darkgreen", "smoothed acceleration"="orangered", "smoothed position"="darkred"))
    # p <- p + geom_point(data=df, aes(x=x, y=y))
    # p <- p + geom_hline(yintercept=0)
    # p <- p + geom_vline(xintercept=0)
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)
    p <- p + theme(legend.title = element_blank()) 

    # for(i in seq(from=1, to=nObs, by=byEllipsePlot)) {
    #     ellipse <- getNormalEllipse(mu=filterRes$xHat[,i], 
    #                                 covar=filterRes$SigmaXHat[[i]], 
    #                                 nPoints=ellipseNPoints, 
    #                                 criticalValue=ellipseCriticalValue)
    #     ellipseDF <- data.frame(x=ellipse[,1], y=ellipse[,2])
    #     p <- p + geom_point(data=ellipseDF, aes(x=x, y=y), colour=ellipseCol, 
    #                         size=ellipsePointSize)
    # }

    ggsave(filename=smoothResFigFilename, plot=p)

    print(p)

    browser()
}

processAll()
