
require(ramcmc)
require(astsa)
require(plotly)
require(mvtnorm)
require(gridExtra)
source("../src/squareRootKF.R")
source("../src/smoothLDS.R")
source("../src/estimateKFInitialCondPPCA.R")
# source("makeSymmetricMatrices.R")

processAll <- function() {
    nFactors <- 2
    maxIter <- 30000
    tol <- 1e-8
    simulationFilename <- "results/simulationCircle.RData"
    figFilename <- "figures/astsaLearningCircle.png"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
    zsForPPCA <- t(as.matrix(zs))
    initialConds <- estimateKFInitialCondPPCA(z=zsForPPCA, nFactors=nFactors)

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
    C0 <- initialConds$C
    # C0 <- C
    B0 <- B
    D0 <- D
    Sigma0 <- initialConds$Sigma
    SRSigmaV0 <- chol(x=Sigma0)
    xHat00 <- xHat0
    V00 <- V0
    SRSigmaX00 <- chol(x=V00)

    emRes = EM0(num=ncol(zs), y=t(zs), A=C, mu0=xHat00, Sigma0=V00, Phi=A0, cQ=SRSigmaW0, cR=SRSigmaV0, max.iter=maxIter, tol=tol)

    df <- data.frame(x=1:length(emRes$like), 
                     y=emRes$like)
    p1 <- ggplot(df, aes(x=x, y=y))
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_point()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("Log Likelihood")

    fRes <- squareRootKF(A=emRes$Phi, B=B, C=C0, D=D, xHat0=emRes$mu0, SRSigmaX0=chol(x=emRes$Sigma0), SRSigmaW=chol(emRes$Q), SRSigmaV=chol(emRes$R), us=us, zs=zs)
    sRes <- smoothLDS(A=emRes$Phi, mu=fRes$xHat, V=fRes$SigmaXHat, P=fRes$SigmaX[2:length(fRes$SigmaX)])
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
