library(ggplot2)
library(gridExtra)
library(R.matlab)

source("../src/simulateLDS.R")

processAll <- function() {
    V0std <- 0.1
    nObs <- 2000
    xlab <- "x"
    ylab <- "y"
    simulationFilename <- "results/simulationRobot.RData"
    simulationFigFilename <- "figures/simulationRobot.png"
    dataFilename <- "~/dev/research/programs/github/python/pykalman/pykalman/datasets/data/robot.mat"

    simRes <- readMat(dataFilename)
    A <- simRes$A
    Gamma <- 10*diag(5)
    C <- simRes$C
    Sigma <- 10*diag(2)
    mu0 <- simRes$x0
    V0 <- simRes$V.0

    nLatents <- nrow(A)

    res <- simulateLDS(nObs=nObs, A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0)
    simulationRes <- c(res, list(A=A, Gamma=Gamma, C=C, Sigma=Sigma, mu0=mu0, V0=V0))
    save(simulationRes, file=simulationFilename)

    df1 <- data.frame(x=res$x[1,], y=res$x[2,])
    p1 <- ggplot(data=df1, aes(x=x, y=y))
    p1 <- p1 + geom_point(data=df1, aes(y=y))
    p1 <- p1 + geom_line(data=df1, aes(y=y))
    p1 <- p1 + geom_hline(yintercept=0)
    p1 <- p1 + geom_vline(xintercept=0)
    p1 <- p1 + xlab(xlab)
    p1 <- p1 + ylab(ylab)
    p1 <- p1 + ggtitle("Observations")
    p1 <- p1 + theme(legend.title = element_blank()) 

    latentID <- c()
    allZ <- c()
    for(i in 1:nrow(res$z)) {
        latentID <- c(latentID, rep(i, length(res$z[i,])))
        allZ <- c(allZ, res$z[i,])
    }
    df2 <- data.frame(x=1:ncol(res$z), y=allZ, latentID=latentID)
    p2 <- ggplot(data=df2, aes(x=x, y=y, color=factor(latentID)))
    p2 <- p2 + geom_point(data=df2, aes(x=x, y=y))
    p2 <- p2 + geom_hline(yintercept=0)
    p2 <- p2 + geom_vline(xintercept=0)
    p2 <- p2 + xlab("Latent Value")
    p2 <- p2 + ylab("Time")
    p2 <- p2 + ggtitle("Latents")
    p2 <- p2 + theme(legend.title = element_blank()) 

    grid.arrange(p1, p2, ncol = 1)
    p <- arrangeGrob(p1, p2, ncol=1)
    ggsave(filename=simulationFigFilename, plot=p)

    browser()
}

processAll()
