
require(astsa)
require(MASS)
require(ramcmc)
require(plotly)
require(mvtnorm)
require(gridExtra)
require(reshape2)
source("../src/squareRootKF.R")
source("../src/smoothLDS.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/plotTrueInitialAndEstimatedMatrices.R")
source("../src/plotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    nFactors <- 2
    maxIter <- 100
    tol <- 1e-8
    simulationFilename <- "results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
    zsForFA <- t(as.matrix(zs))
    initialConds <- estimateKFInitialCondFA(z=zsForFA, nFactors=nFactors)

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
    Gamma0 <- 1e-3*diag(rep(1, ncol(A0)))
    SRSigmaW0 <- chol(x=Gamma0)
    C0 <- initialConds$C
    # C0 <- C
    B0 <- B
    D0 <- D
    Sigma0 <- diag(initialConds$sigmaDiag)
    SRSigmaV0 <- chol(x=Sigma0)
    xHat00 <- xHat0
    V00 <- V0
    SRSigmaX00 <- chol(x=V00)

    emRes = EM0(num=ncol(zs), y=t(zs), A=C, mu0=xHat00, Sigma0=V00, Phi=A0, cQ=SRSigmaW0, cR=SRSigmaV0, max.iter=maxIter, tol=tol)

    df <- data.frame(x=1:length(emRes$like), 
                     y=emRes$like)
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + xlab("Time")
    p <- p + ylab("Log Likelihood")
    p <- ggplotly(p)
    llFigFilename <- "figures//circleASTSA_LogLik.html"
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(llFigFilename)), basename(llFigFilename)))
    print(p)

    AFigFilename <- "figures//circleASTSA_A.html"
    plotTrueInitialAndEstimatedMatrices(trueM=A, initialM=A0, estimatedM=emRes$Phi, title="A", figFilename=AFigFilename)

    CFigFilename <- "figures//circleASTSA_C.html"
    plotTrueInitialAndEstimatedMatrices(trueM=C, initialM=C0, title="C", figFilename=CFigFilename)

    GammaFigFilename <- "figures//circleASTSA_Gamma.html"
    plotTrueInitialAndEstimatedMatrices(trueM=Gamma, initialM=Gamma0, estimatedM=emRes$Q, title="Gamma", figFilename=GammaFigFilename)

    SigmaFigFilename <- "figures//circleASTSA_Sigma.html"
    plotTrueInitialAndEstimatedMatrices(trueM=Sigma, initialM=Sigma0, estimatedM=emRes$R, title="Sigma", figFilename=SigmaFigFilename)

    V0FigFilename <- "figures//circleASTSA_V0.html"
    plotTrueInitialAndEstimatedMatrices(trueM=V0, initialM=V00, estimatedM=emRes$Sigma0, title="V0", figFilename=V0FigFilename)

    xHat0FigFilename <- "figures//circleASTSA_XHat0.html"
    plotTrueInitialAndEstimatedVectors(trueV=xHat0, initialV=xHat00, estimatedV=emRes$mu0, title="xHat0", figFilename=xHat0FigFilename)

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
    p <- ggplot(data, aes(x=sample, y=latent, 
                           color=factor(latentID),
                           linetype=factor(latentType)))
    p <- p + geom_line()
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + ylab("Latent Value")
    p <- p + xlab("Time")
    p <- p + theme(legend.title = element_blank()) 
    p <- ggplotly(p)
    latentsFigFilename <- "figures//circleASTSA_Latents.html"
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(latentsFigFilename)), basename(latentsFigFilename)))
    print(p)

    browser()
}

processAll()
