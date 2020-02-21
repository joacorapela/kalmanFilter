require(ggplot2)
require(R.matlab)

source("../src/simulateLDS.R")

processAll <- function() {
    measurementNoiseSD <- 1e-1
    stateNoiseSD <- 1e-1
    nObs <- 1000
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulation2DTrajectory.RData"
    simulationMatFilename <- "results/simulation2DTrajectory.mat"
    simulationFigFilename <- "figures/simulation2DTrajectory.png"

    # Modified simulation parameters from 
    # https://www.cs.utexas.edu/~teammco/misc/kalman_filter/

    # state transition
      A <-   t(matrix(c(1.0, 0.0, 0.01, 0.0, 
                        0.0, 1.0, 0.0, 0.01, 
                        0.0, 0.0, 1.0, 0.0, 
                        0.0, 0.0, 0.0, 1.0), nrow=4))
    # state noise
    Gamma <- t(matrix(c(1.0, 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0, 0.0, 
                        0.0, 0.0, 1.0, 0.0, 
                        0.0, 0.0, 0.0, 1.0), nrow=4))*stateNoiseSD
    # state-measurement transfer
    C <-     t(matrix(c(1.0, 0.0, 0.0, 0.0, 
                        0.0, 1.0, 0.0, 0.0), ncol=2))

    # measurement noise
    Sigma <- t(matrix(c(1.0, 0.0,
                        0.0, 1.0), nrow=2))*measurementNoiseSD

    nLatents <- nrow(A)
    mu0 <- rep(0, times=nrow(A))
    V0 <- diag(x=rep(10, times=nrow(A)))

    res <- simulateLDS(nObs=nObs, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    simulationRes <- c(res, list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0))
    save(simulationRes, file=simulationFilename)
    writeMat(con=simulationMatFilename, x=res$x, z=res$z, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)

    df <- data.frame(t(cbind(res$x[1:2,], res$z[1:2,])))
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
