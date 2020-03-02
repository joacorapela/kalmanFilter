require(ggplot2)
require(plotly)
require(mvtnorm)
require(ramcmc)                                                       

source("../src/lsolve.R")
source("../src/squareRootKF.R")
source("../src/chol_downdate_higherOrder.R")
source("getNormalEllipse.R")

processAll <- function() {
    ellipseNPoints <- 100
    ellipseCriticalValue <- .95
    ellipseCol <- "red"
    propByEllipsePlot <- .01
    ellipsePointSize <- .1
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulation2DTrajectory.RData"
    filterResFilename <- "results/srFilterResSimulation2DTrajectory.RData"
    filterResFigFilename <- "figures/srFilterResSimulation2DTrajectory.png"

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

    # filterRes <- filterLDS(x=x, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    B <- matrix(0, nrow=nrow(A), ncol=1)
    D <- matrix(0, nrow=nrow(C), ncol=1)
    us <- matrix(0, nrow=1, ncol=ncol(x))
    SRSigmaX0=chol(x=V0)
    SRSigmaW=chol(x=Gamma)
    SRSigmaV=chol(x=Sigma)
    filterRes <- squareRootKF(A=A, B=B, C=C, D=D, xHat0=mu0, SRSigmaX0=SRSigmaX0, SRSigmaW=SRSigmaW, SRSigmaV=SRSigmaV, us=us, zs=x)
    save(filterRes, file=filterResFilename)

    df <- data.frame(x=c(simRes$x[1,], simRes$z[1,], filterRes$xHat[1,]), 
                     y=c(simRes$x[2,], simRes$z[2,], filterRes$xHat[2,]),
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
        ellipse <- getNormalEllipse(mu=filterRes$xHat[,i], 
                                    covar=filterRes$SigmaXHat[[i]], 
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
