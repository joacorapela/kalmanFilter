
require(MASS)
require(ggplot2)
require(R.matlab)

source("../src/simulateLDS.R")

processAll <- function() {
    # solution \theta(t) = \theta_0 cos(\omega t), with \omega=\sqrt{g/l}
    initialFilterSD <- 1e-2
    stateNoiseSD <- 1e-1
    measurementNoiseSD <- 5e-1
    dt <- 1e-2
    omega2 <- 10
    nCells <- 20
    T <- 4*pi/sqrt(omega2)
    mu0 <- c(0, pi/6)
    V0 <- initialFilterSD*diag(c(1,1))
    A <- rbind(c(1, -omega2*dt), c(dt, 1))
    Gamma <- stateNoiseSD*diag(c(1,1))
    C <- matrix(c(rep(0, times=nCells), rep(1, times=nCells)), ncol=2)
    Sigma <- measurementNoiseSD*diag(rep(1, times=nCells))
    xlab <- "Time"
    ylab <- ""
    simulationFilename <- "results/simulationPendulum.RData"
    simulationMatFilename <- "results/simulationPendulum.mat"
    simulationFigFilename <- "figures/simulationPendulum.png"

    t <- seq(from=0, to=T, by=dt)
    nObs <- length(t)

    res <- simulateLDS(nObs=nObs, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    simulationRes <- c(res, list(t=t, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0))
    save(simulationRes, file=simulationFilename)
    writeMat(con=simulationMatFilename, t=t,  x=res$x, z=res$z, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)

    df <- data.frame(x=c(t, t, t), y=c(res$x[1,], res$z[1,], res$z[2,]), type=factor(c(rep("spike rate", nObs), rep("latent aceleration", nObs), rep("latent position", nObs))))
    p <- ggplot(data=df, aes(x=x, y=y, colour=type))
    p <- p + geom_point(data=df, aes(x=x, y=y))
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
