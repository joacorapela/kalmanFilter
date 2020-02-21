
require(mvtnorm)
require(gridExtra)
source("../src/filterLDS.R")
source("../src/smoothLDS.R")
source("../src/emEstimationKF.R")
source("makeSymmetricMatrices.R")

processAll <- function() {
    initialCondNoiseSD <- 100.0
    nIter <- 5
    simulationFilename <- "results/simulationRobot.RData"
    figFilename <- "results/learningRobot.png"

    simRes <- get(load(simulationFilename))
    x <- simRes$x
    z <- simRes$z
    A <- simRes$A
    Gamma <- simRes$Gamma
    C <- simRes$C
    Sigma <- simRes$Sigma
    mu0 <- simRes$mu0
    V0 <- simRes$V0

    # A0 <- diag(ncol(A))
    # A0 <- A+matrix(rnorm(n=length(A), sd=initialCondNoiseSD), ncol=ncol(A))
    A0 <- A
    # Gamma0 <- Gamma+makeSymmetricMatrix(m=matrix(rnorm(n=length(Gamma), sd=initialCondNoiseSD), ncol=ncol(Gamma)))
    Gamma0 <- Gamma
    # C0 <- C+matrix(rnorm(n=length(C), sd=initialCondNoiseSD), ncol=ncol(C))
    C0 <- C
    # Sigma0 <- Sigma+makeSymmetricMatrix(m=matrix(rnorm(n=length(Sigma), sd=sd(Sigma)/initialCondSNR), ncol=ncol(Sigma)))
    Sigma0 <- Sigma
    mu00 <- mu0+rnorm(n=length(mu0), sd=initialCondNoiseSD)
    # mu00 <- mu0
    # V00 <- V0+makeSymmetricMatrix(m=matrix(rnorm(n=length(V0), sd=initialCondNoiseSD), ncol=ncol(V0)))
    V00 <- V0

    eRes <- emEstimationKF(x=x, A0=A0, Gamma0=Gamma0, C0=C0, Sigma0=Sigma0, mu00=mu00, V00=V00, nIter=nIter)

    df <- data.frame(x=1:length(eRes$logLikelihoods), 
                     y=eRes$logLikelihoods)
    p1 <- ggplot(df, aes(x=x, y=y))
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_point()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("Log Likelihood")

    fRes <- filterLDS(x=x, A=eRes$A, Gamma=eRes$Gamma, C=eRes$C, Sigma=eRes$Sigma, mu0=eRes$mu0, V0=eRes$V0)
    sRes <- smoothLDS(A=eRes$A, mu=fRes$mu, V=fRes$V, P=fRes$P)

    data <- data.frame()
    for(i in 1:nrow(z)) {
        dataBlock <- data.frame(sample=1:length(z[i,]), 
                                latent=z[i,],
                                latentID=rep(i, length(z[i,])),
                                latentType=rep("true", 
                                               length(z[i,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:nrow(sRes$muHat)) {
        dataBlock <- data.frame(sample=1:length(sRes$muHat[i,]), 
                                latent=sRes$muHat[i,],
                                latentID=rep(i, length(sRes$muHat[i,])),
                                latentType=rep("estimated", 
                                               length(sRes$muHat[i,])))
        data <- rbind(data, dataBlock)
    }
    p2 <- ggplot(data, aes(x=sample, y=latent, 
                           color=factor(latentID),
                           linetype=factor(latentType)))
    p2 <- p2 + geom_line()
    p2 <- p2 + geom_hline(yintercept=0)
    p2 <- p2 + geom_vline(xintercept=0)
    p2 <- p2 + ylab("Latent Value")
    p2 <- p2 + xlab("Time")
    p2 <- p2 + theme(legend.title = element_blank()) 

    grid.arrange(p1, p2, ncol = 1)
    p <- arrangeGrob(p1, p2, ncol=1)
    ggsave(filename=figFilename, plot=p)

    browser()
}

processAll()
