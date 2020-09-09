require(plotly)
require(ggplot2)
require(R.matlab)

source("../src/simulateLDS.R")

processAll <- function() {
    measurementNoiseSD <- 5*1e-2
    stateNoiseSD <- 1e-3
    dObs <- 2*1
    nObs <- 2000
    latentsVarPrior <- 10
    xlab <- "x"
    ylab <- "y"
    simulationFilenamePattern <- "results/simulationCircleDObs%02d.RData"
    simulationFigFilenamePattern <- "figures/simulationCircleDObs%02d.png"

    simulationFilename <- sprintf(simulationFilenamePattern, dObs)
    simulationFigFilename <- sprintf(simulationFigFilenamePattern, dObs)

    # From Fig 4.5b of http://www.staff.city.ac.uk/o.castro-alvaredo/dynamical/dynamicalsystems.pdf

    # state transition
    Atmp <- matrix(c(-128, 80, -272, 128), ncol=2)
    dt <- 1e-3
    x0 <- c(1,1)
    A <- dt*Atmp + diag(2)

    # state noise
    Gamma <- diag(x=nrow(A))*stateNoiseSD

    # state-measurement transfer
    CLatent1Col <- c(rep(1, times=1/2*dObs), rep(0, times=1/2*dObs))
    # CLatent1Col <- runif(dObs)
    CLatent2Col <- 1-CLatent1Col
    C <- cbind(CLatent1Col, CLatent2Col)

    # measurement noise
    Sigma <- diag(x=dObs)*measurementNoiseSD

    nLatents <- nrow(A)
    mu0 <- c(5, 5)

    res <- simulateLDS(nObs=nObs, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0)
    simulationRes <- c(res, list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0))
    save(simulationRes, file=simulationFilename)

    df <- data.frame(t(cbind(res$x[c(1, dObs/2+1),], res$z)))
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

    p <- ggplotly(p)
    print(p)

    browser()
}

processAll()
