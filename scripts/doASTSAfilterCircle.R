require(ggplot2)
require(plotly)
require(mvtnorm)
require(ramcmc)                                                       

source("../src/lsolve.R")
source("../src/squareRootKF.R")
source("../src/chol_downdate_higherOrder.R")
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
    filterResFigFilename <- "figures/srFilterResSimulationCircle.png"

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

    SRSigmaW=chol(x=Gamma)
    SRSigmaV=chol(x=Sigma)
    filterRes <- Kfilter0(num=ncol(x), y=t(x), A=C, mu0=mu0, Sigma0=V0, Phi=A, cQ=t(SRSigmaW), cR=t(SRSigmaV))

    save(filterRes, file=filterResFilename)

    df <- data.frame(x=c(simRes$x[1,], simRes$z[1,], filterRes$xf[1,1,]), 
                     y=c(simRes$x[dObs,], simRes$z[2,], filterRes$xf[2,1,]),
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
        ellipse <- getNormalEllipse(mu=filterRes$xf[,,i], 
                                    covar=matrix(filterRes$Pf[,,i], 
                                                 nrow=nrow(z)), 
                                    nPoints=ellipseNPoints, 
                                    criticalValue=ellipseCriticalValue)
        ellipseDF <- data.frame(x=ellipse[,1], y=ellipse[,2])
        p <- p + geom_point(data=ellipseDF, aes(x=x, y=y), colour=ellipseCol, 
                            size=ellipsePointSize)
    }

    ggsave(filename=filterResFigFilename, plot=p)

    p <- ggplotly(p)
    print(p)

    browser()
}

processAll()
