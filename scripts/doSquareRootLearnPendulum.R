
require(ggplot2)
require(mvtnorm)
require(ramcmc)                                                       
require(gridExtra)
source("../src/squareRootKF.R")
source("../src/smoothLDS.R")
source("../src/emEstimationSquareRootKF.R")
# source("makeSymmetricMatrices.R")

processAll <- function() {
    initialCondNoiseSD <- 0.01
    nIter <- 10
    simulationFilename <- "results/simulationPendulum.RData"
    figFilename <- "figures/srLearningPendulum.png"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
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

    A0 <- A+matrix(rnorm(n=length(A), sd=initialCondNoiseSD), ncol=ncol(A))
    # A0 <- A
    noiseMatrix <- matrix(rnorm(n=length(Gamma), sd=initialCondNoiseSD), ncol=ncol(Gamma))
    noiseMatrix <- noiseMatrix%*%t(noiseMatrix)
    Gamma0 <- Gamma + noiseMatrix
    # Gamma0 <- Gamma
    SRSigmaW0 <- chol(x=Gamma0)
    C0 <- C+matrix(rnorm(n=length(C), sd=initialCondNoiseSD), ncol=ncol(C))
    # C0 <- C
    # B0 <- C+matrix(rnorm(n=length(B), sd=initialCondNoiseSD), ncol=ncol(B))
    B0 <- B
    # D0 <- C+matrix(rnorm(n=length(D), sd=initialCondNoiseSD), ncol=ncol(D))
    D0 <- D
    # Sigma0 <- Sigma+makeSymmetricMatrix(m=matrix(rnorm(n=length(Sigma), sd=sd(Sigma)/initialCondSNR), ncol=ncol(Sigma)))
    Sigma0 <- Sigma
    SRSigmaV0 <- chol(x=Sigma0)
    # xHat00 <- xHat0+rnorm(n=length(xHat0), sd=initialCondNoiseSD)
    xHat00 <- xHat0
    # V00 <- V0+makeSymmetricMatrix(m=matrix(rnorm(n=length(V0), sd=initialCondNoiseSD), ncol=ncol(V0)))
    V00 <- V0
    SRSigmaX00 <- chol(x=V00)

    res <- squareRootKF(A=A, B=B, C=C, D=D, xHat0=xHat0, SRSigmaX0=SRSigmaX0, SRSigmaW=SRSigmaW, SRSigmaV=SRSigmaV, us=us, zs=zs)
    show(sprintf("True data log likelihood: %.4f", sum(log(res$c))))

    eRes <- emEstimationSquareRootKF(zs=zs, A0=A0, SRSigmaW0=SRSigmaW0, C0=C0, SRSigmaV0=SRSigmaV0, B0=B0, D0=D0, xHat00=xHat00, SRSigmaX0=SRSigmaX00, nIter=nIter, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE))

    df <- data.frame(x=1:length(eRes$logLikelihoods), 
                     y=eRes$logLikelihoods)
    p1 <- ggplot(df, aes(x=x, y=y))
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_point()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("Log Likelihood")

    fRes0 <- squareRootKF(A=A0, B=B0, C=C0, D=D0, xHat0=xHat00, SRSigmaX0=SRSigmaX00, SRSigmaW=SRSigmaW0, SRSigmaV=SRSigmaV0, us=us, zs=zs)
    sRes0 <- smoothLDS(A=A0, mu=fRes0$xHat, V=fRes0$SigmaXHat, P=fRes0$SigmaX[2:length(fRes0$SigmaX)])

    fRes <- squareRootKF(A=eRes$A, B=B, C=eRes$C, D=D, xHat0=eRes$xHat0, SRSigmaX0=chol(x=eRes$V0), SRSigmaW=chol(x=eRes$Gamma), SRSigmaV=chol(x=eRes$Sigma), us=us, zs=zs)
    sRes <- smoothLDS(A=eRes$A, mu=fRes$xHat, V=fRes$SigmaXHat, P=fRes$SigmaX[2:length(res$SigmaX)])

    nObs <- length(simRes$t)
    stateVariances <- matrix(NA, nrow=nrow(A), ncol=nObs)
    for(i in 1:nObs) {
        stateVariances[,i] <- diag(sRes$VHat[[i]])
    }

    dfPosTrue <- data.frame(t=simRes$t, value=simRes$z[2,], measureType=rep("position", times=length(simRes$t)), estimateType=rep("true", times=length(simRes$t)))
    dfPosInit <- data.frame(t=simRes$t, value=sRes0$muHat[2,], measureType=rep("position", times=length(simRes$t)), estimateType=rep("init", times=length(simRes$t)))
    dfAccTrue <- data.frame(t=simRes$t, value=simRes$z[1,], measureType=rep("acceleration", times=length(simRes$t)), estimateType=rep("true", times=length(simRes$t)))
    dfAccInit <- data.frame(t=simRes$t, value=sRes0$muHat[1,], measureType=rep("acceleration", times=length(simRes$t)), estimateType=rep("init", times=length(simRes$t)))
    dfPosHat <- data.frame(t=simRes$t, value=sRes$muHat[2,], measureType=rep("position", times=length(simRes$t)), estimateType=rep("hat", times=length(simRes$t)))
    dfAccHat <- data.frame(t=simRes$t, value=sRes$muHat[1,], measureType=rep("acceleration", times=length(simRes$t)), estimateType=rep("hat", times=length(simRes$t)))
    df <- rbind(dfPosTrue, dfPosInit, dfPosHat, dfAccTrue, dfAccInit, dfAccHat)

    dfAccHatBounds <- data.frame(t=simRes$t, lower=sRes$muHat[1,]-1.96*sd(stateVariances[1,]), upper=sRes$muHat[1,]+1.96*sd(stateVariances[1,]), measureType=rep("acceleration", times=length(simRes$t)))
    dfPosHatBound <- data.frame(t=simRes$t, lower=sRes$muHat[2,]-1.96*sd(stateVariances[2,]), upper=sRes$muHat[2,]+1.96*sd(stateVariances[2,]), measureType=rep("position", times=length(simRes$t)))
    dfHatBounds <- rbind(dfAccHatBounds, dfPosHatBound)

    p2 <- ggplot()
    p2 <- p2 + geom_line(data=df, aes(x=t, y=value, colour=measureType, linetype=estimateType))
    p2 <- p2 + geom_ribbon(data=dfHatBounds, aes(x=t, ymin=lower, ymax=upper, fill=measureType), alpha=.7, show.legend=FALSE)
    p2 <- p2 + xlab("Time")
    p2 <- p2 + ylab("")
    p2 <- p2 + theme(legend.title = element_blank()) 
    p2 <- ggplotly(p2)
    print(p2)
    #
    # grid.arrange(p1, p2, ncol = 1)
    # p <- arrangeGrob(p1, p2, ncol=1)
    # ggsave(filename=figFilename, plot=p)

    browser()
}

processAll()
