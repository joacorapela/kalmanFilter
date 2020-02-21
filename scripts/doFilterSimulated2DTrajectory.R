library(ggplot2)
library(mvtnorm)

source("../src/filterLDS.R")
source("getNormalEllipse.R")

processAll <- function() {
    ellipseNPoints <- 100
    ellipseCriticalValue <- .95
    ellipseCol <- "red"
    propByEllipsePlot <- .10
    ellipsePointSize <- .1
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulation2DTrajectory.RData"
    filterResFilename <- "results/filterResSimulation2DTrajectory.RData"
    filterResFigFilename <- "figures/filterResSimulation2DTrajectory.png"

    simRes <- get(load(simulationFilename))
    x <- simRes$x
    z <- simRes$z
    A <- simRes$A
    Gamma <- simRes$Gamma
    C <- simRes$C
    Sigma <- simRes$Sigma
    mu0 <- simRes$mu0
    V0 <- simRes$V0
    nObs <- ncol(x)
    byEllipsePlot <- propByEllipsePlot*nObs

    filterRes <- filterLDS(x=x, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    save(filterRes, file=filterResFilename)

    df <- data.frame(x=c(simRes$x[1,], simRes$z[1,], filterRes$mu[1,]), 
                     y=c(simRes$x[2,], simRes$z[2,], filterRes$mu[2,]),
                     type=factor(c(rep("measurement", nObs), rep("latent", nObs), rep("filtered", nObs))))
    p <- ggplot()
    p <- p + geom_point(data=df, aes(x=x, y=y, colour=type))
    p <- p + scale_colour_manual(values=c("measurement"="blue", 
                                          "latent"="green", 
                                          "filtered"="red"))
    # p <- p + geom_point(data=df, aes(x=x, y=y))
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)
    p <- p + theme(legend.title = element_blank()) 

    for(i in seq(from=1, to=nObs, by=byEllipsePlot)) {
        ellipse <- getNormalEllipse(mu=filterRes$mu[,i], 
                                    covar=filterRes$V[[i]], 
                                    nPoints=ellipseNPoints, 
                                    criticalValue=ellipseCriticalValue)
        ellipseDF <- data.frame(x=ellipse[,1], y=ellipse[,2])
        p <- p + geom_point(data=ellipseDF, aes(x=x, y=y), colour=ellipseCol, 
                            size=ellipsePointSize)
    }

    ggsave(filename=filterResFigFilename, plot=p)

    print(p)

    browser()
}

processAll()
