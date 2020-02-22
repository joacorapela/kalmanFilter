require(ggplot2)
require(R.matlab)

source("../src/simulateLDS.R")

processAll <- function() {
    measurementNoiseSD <- 1e-3
    stateNoiseSD <- 1e-3
    dObs <- 4*4
    nObs <- 1000
    latentsVarPrior <- 10
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulation2DTrajectory.RData"
    simulationMatFilename <- "results/simulation2DTrajectory.mat"
    simulationFigFilename <- "figures/simulation2DTrajectory.png"

    # Modified simulation parameters from 
    # https://www.cs.utexas.edu/~teammco/misc/kalman_filter/

    # state transition
    A <-   t(matrix(c(1.0, 0.0, 0.01, 0.0, 
                      0.0, 1.0, 0.0,  0.01, 
                      0.0, 0.0, 1.0,  0.0, 
                      0.0, 0.0, 0.0,  1.0), nrow=4))
    # state noise
    Gamma <- diag(x=nrow(A))*stateNoiseSD
    # state-measurement transfer
    # 25% of the cells are tunned to position on x
    # 25% of the cells are tunned to position on y
    # 25% of the cells are tunned to velocity on x
    # 25% of the cells are tunned to velocity on y
    CPosXCol <- c(rep(1, times=1/4*dObs), rep(0, times=3/4*dObs))
    CPosYCol <- c(rep(0, times=1/4*dObs), rep(1, times=1/4*dObs), rep(0, times=2/4*dObs))
    CVelXCol <- c(rep(0, times=2/4*dObs), rep(1, times=1/4*dObs), rep(0, times=1/4*dObs))
    CVelYCol <- c(rep(0, times=3/4*dObs), rep(1, times=1/4*dObs))
    C <-     cbind(CPosXCol, CPosYCol, CVelXCol, CVelYCol)

    # measurement noise
    Sigma <- diag(x=dObs)*measurementNoiseSD

    nLatents <- nrow(A)
    mu0 <- rep(0, times=nrow(A))
    V0 <- diag(x=rep(latentsVarPrior, times=nrow(A)))

    res <- simulateLDS(nObs=nObs, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    simulationRes <- c(res, list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0))
    save(simulationRes, file=simulationFilename)
    writeMat(con=simulationMatFilename, x=res$x, z=res$z, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)

    df <- data.frame(t(cbind(res$x[c(1, dObs/4+1),], res$z[1:2,])))
    df <- cbind(df, c(rep("measurement", nObs), rep("latent", nObs)))
    colnames(df) <- c("x", "y", "type")
    p <- ggplot(data=df, aes(x=x, y=y, color=factor(type)))
    p <- p + geom_point()
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)
    p <- p + theme(legend.title = element_blank()) 
    ggsave(filename=simulationFigFilename, plot=p)

    print(p)

    browser()
}

processAll()
