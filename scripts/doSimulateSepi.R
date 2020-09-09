require(plotly)
require(ggplot2)
require(R.matlab)

source("../src/simulateLDS.R")

processAll <- function() {
    measurementNoiseSD <- 1e-3
    stateNoiseSD <- 1e-3
    dObs <- 10
    f <- 1.0
    dt <- 1e-3
    T <- 1.0
    xlabLatents <- "Time (sec)"
    ylabLatents <- "Latents"
    xlabNeurons <- "Time (sec)"
    ylabNeurons <- "Firing Rate"
    simulationFilenamePattern <- "results/simulationSepiDObs%02d.RData"
    simulationFigFilenamePattern <- "figures/simulationSepiDObs%02d.png"

    simulationFilename <- sprintf(simulationFilenamePattern, dObs)
    simulationFigFilename <- sprintf(simulationFigFilenamePattern, dObs)

    nObs <- round(T/dt)

    # state transition
    A <- matrix(c(1.0, -(2*pi*f)^2*dt, 1*dt, 1.0), ncol=2)

    # state noise
    Gamma <- diag(x=nrow(A))*stateNoiseSD

    # state-measurement transfer
    C <- matrix(runif(n=2*dObs, min=-1, max=1), ncol=2)

    # measurement noise
    Sigma <- diag(x=dObs)*measurementNoiseSD

    nLatents <- nrow(A)
    x0 <- c(0, 2*pi*f)

    res <- simulateLDS(nObs=nObs, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=x0)
    simulationRes <- c(res, list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=x0))
    save(simulationRes, file=simulationFilename)

    t <- ((1:nObs)-1)*dt

    df <- data.frame(t=t, vel=res$z[1,], acel=res$z[2,])
    dfM <- melt(df, id.vars="t")
    p <- ggplot(data=dfM, aes(x=t, y=value, color=variable))
    p <- p + geom_point()
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + xlab(xlabLatents)
    p <- p + ylab(ylabLatents)
    p <- p + theme(legend.title = element_blank()) 
    ggsave(filename=simulationFigFilename, plot=p)
    p <- ggplotly(p)
    print(p)

    df <- data.frame(t(rbind(t, res$x)))
    colnames(df) <- c("t", sprintf("neuron %02d", 1:(ncol(df)-1)))
    dfM <- melt(df, id.vars="t")
    p <- ggplot(data=dfM, aes(x=t, y=value, color=variable))
    p <- p + geom_point()
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + xlab(xlabNeurons)
    p <- p + ylab(ylabNeurons)
    p <- p + theme(legend.title = element_blank()) 
    ggsave(filename=simulationFigFilename, plot=p)
    p <- ggplotly(p)
    print(p)

    browser()
}

processAll()
