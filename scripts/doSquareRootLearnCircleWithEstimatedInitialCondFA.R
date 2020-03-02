
require(pcaMethods)
require(ramcmc)
require(plotly)
require(mvtnorm)
require(gridExtra)
source("../src/squareRootKF.R")
source("../src/smoothLDS.R")
source("../src/emEstimationSquareRootKF.R")
source("../src/estimateKFInitialCondFA.R")
# source("makeSymmetricMatrices.R")

processAll <- function() {
    initialCondNoiseSD <- 1e-2
    nFactors <- 2
    # nIter <- 10
    nIter <- 1000
    tol <- 1e-3
    simulationFilename <- "results/simulationCircle.RData"
    figFilename <- "figures/srLearningCircle.png"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
    zsForFA <- t(as.matrix(zs))
    initialConds <- estimateKFInitialCondFA(z=zsForFA, nFactors=nFactors, rotation="none")


    A <- simRes$A
    C <- simRes$C
    B <- matrix(0, nrow=nrow(A), ncol=1)
    D <- matrix(0, nrow=nrow(C), ncol=1)
    us <- matrix(0, nrow=1, ncol=ncol(zs))
    Gamma <- simRes$Gamma
    SRSigmaW <- chol(x=Gamma)
    Sigma <- simRes$Sigma
    SRSigmaV <- chol(x=Sigma)
    xHat0 <- simRes$mu0
    V0 <- simRes$V0
    SRSigmaX0 <- chol(x=V0)

    A0 <- initialConds$A
    # A0 <- A
    Gamma0 <- initialConds$Gamma
    SRSigmaW0 <- chol(x=Gamma0)
    # C0 <- C+matrix(rnorm(n=length(C), sd=initialCondNoiseSD), ncol=ncol(C))
    # C0 <- initialConds$C
    C0 <- C
    B0 <- B
    D0 <- D
    Sigma0 <- initialConds$Sigma
    SRSigmaV0 <- chol(x=Sigma0)
    xHat00 <- xHat0
    V00 <- V0
    SRSigmaX00 <- chol(x=V00)

    res <- squareRootKF(A=A, B=B, C=C, D=D, xHat0=xHat0, SRSigmaX0=SRSigmaX0, SRSigmaW=SRSigmaW, SRSigmaV=SRSigmaV, us=us, zs=zs)
    logMsgPattern <- sprintf("loglik(iter=%%04d)=%%.4f (%.4f)", sum(log(res$c)))

    eRes <- emEstimationSquareRootKF(zs=zs, A0=A0, SRSigmaW0=SRSigmaW0, C0=C0, SRSigmaV0=SRSigmaV0, B0=B0, D0=D0, xHat00=xHat00, SRSigmaX0=SRSigmaX00, nIter=nIter, tol=tol, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=FALSE, observationCovariance=TRUE), logMsgPattern=logMsgPattern)

    df <- data.frame(x=1:length(eRes$logLikelihoods), 
                     y=eRes$logLikelihoods)
    p1 <- ggplot(df, aes(x=x, y=y))
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_point()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("Log Likelihood")

    fRes <- squareRootKF(A=eRes$A, B=B, C=eRes$C, D=D, xHat0=eRes$xHat0, SRSigmaX0=chol(x=eRes$V0), SRSigmaW=chol(x=eRes$Gamma), SRSigmaV=chol(x=eRes$Sigma), us=us, zs=zs)
    sRes <- smoothLDS(A=eRes$A, mu=fRes$xHat, V=fRes$SigmaXHat, P=fRes$SigmaX[2:length(fRes$SigmaX)])
    fRes0 <- squareRootKF(A=A0, B=B0, C=C0, D=D0, xHat0=xHat00, SRSigmaX0=chol(x=V00), SRSigmaW=chol(x=Gamma0), SRSigmaV=chol(x=Sigma0), us=us, zs=zs)
    sRes0 <- smoothLDS(A=A0, mu=fRes0$xHat, V=fRes0$SigmaXHat, P=fRes0$SigmaX[2:length(fRes0$SigmaX)])

    data <- data.frame()
    for(i in 1:nrow(simRes$z)) {
        dataBlock <- data.frame(sample=1:length(simRes$z[i,]), 
                                latent=simRes$z[i,],
                                latentID=rep(i, length(simRes$z[i,])),
                                latentType=rep("true", 
                                               length(simRes$z[i,])))
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
    for(i in 1:nrow(sRes0$muHat)) {
        dataBlock <- data.frame(sample=1:length(sRes0$muHat[i,]), 
                                latent=sRes0$muHat[i,],
                                latentID=rep(i, length(sRes0$muHat[i,])),
                                latentType=rep("initial", 
                                               length(sRes0$muHat[i,])))
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

    # browser()
    # print(p1)
    # browser()
    # print(p2)
    # browser()
    # grid.arrange(p1, p2, ncol = 1)
    # p <- arrangeGrob(p1, p2, ncol=1)
    # ggsave(filename=figFilename, plot=p)
    p2 <- ggplotly(p2)
    print(p2)
    # plot(unlist(res$logLikelihoods), type="b")
    browser()
}

processAll()
