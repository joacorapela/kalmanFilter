
require(MASS)
require(pcaMethods)
require(ramcmc)
require(plotly)
require(mvtnorm)
require(gridExtra)
source("../src/chol_downdate_higherOrder.R")
source("../src/lsolve.R")
source("../src/squareRootKF.R")
source("../src/smoothLDS.R")
source("../src/emEstimationSquareRootKF.R")
source("../src/estimateKFInitialCondFA.R")
# source("makeSymmetricMatrices.R")

processAll <- function() {
    propTrain <- .5
    nFactors <- 2
    # maxIter <- 10
    maxIter <- 1000
    tol <- 1e-3
    gamma0Var <- 1e4
    simulationFilename <- "results/simulationCircle.RData"
    figFilename <- "figures/srLearningCircle.png"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
    nNeurons <- nrow(zs)
    zsTrain <- zs[,1:round(ncol(zs)*propTrain)]
    zsTest <- zs[,(round(ncol(zs)*propTrain)-1):ncol(zs)]
    trueLatentsTrain <- simRes$z[,1:round(ncol(zs)*propTrain)]
    trueLatentsTest <- simRes$z[,(round(ncol(zs)*propTrain)-1):ncol(simRes$z)]
    # x0 <- simRes$z[,1]
    # V0Var <- 1e-4
    # x0 <- rep(0, times=nrow(simRes$z))
    x0 <- runif(nFactors)*2-1
    initialStateType <- "init01"
    V0Var <- 1e+4
    rm(zs)
    zsForFA <- t(as.matrix(zsTrain))
    initialConds <- estimateKFInitialCondFA(z=zsForFA, nFactors=nFactors, rotation="none")

    usTrain <- matrix(0, nrow=1, ncol=ncol(zsTrain))
    usTest <- matrix(0, nrow=1, ncol=ncol(zsTest))
    Gamma0 <- gamma0Var*diag(rep(1, times=nFactors))
    V0 <- V0Var*diag(rep(1, nFactors))
    SRSigmaX0 <- chol(x=V0)

    A0 <- initialConds$A
    SRSigmaW0 <- chol(x=Gamma0)
    C0 <- initialConds$C
    B0 <- matrix(0, nrow=nFactors, ncol=1)
    D0 <- matrix(0, nrow=nNeurons, ncol=1)
    Sigma0 <- diag(initialConds$sigmaDiag)
    SRSigmaV0 <- chol(x=Sigma0)
    x00 <- x0
    V00 <- V0
    SRSigmaX00 <- chol(x=V00)

    logMsgPattern <- "loglik(iter=%04d)=%.4f"
    eRes <- emEstimationSquareRootKF(zs=zsTrain, 
                                     A0=A0, 
                                     SRSigmaW0=SRSigmaW0, 
                                     C0=C0, 
                                     SRSigmaV0=SRSigmaV0, 
                                     B0=B0, 
                                     D0=D0, 
                                     x00=x00, 
                                     initialStateType=initialStateType,
                                     SRSigmaX0=SRSigmaX00, 
                                     maxIter=maxIter, tol=tol, 
                                     varsToEstimate=list(initialStateMean=TRUE, 
                                                         initialStateCovariance=TRUE, 
                                                         transitionMatrix=TRUE, 
                                                         transitionCovariance=TRUE, 
                                                         observationMatrix=TRUE, 
                                                         observationCovariance=TRUE), 
#                                      covsConstraints=list(SigmaX0="diagonal", 
#                                                           SigmaW="diagonal", 
#                                                           SigmaV="diagonal"), 
                                     covsConstraints=list(SigmaX0="unconstrained", 
                                                          SigmaW="unconstrained", 
                                                          SigmaV="unconstrained"), 
                                     logMsgPattern=logMsgPattern)

    df <- data.frame(x=1:length(eRes$logLikelihoods),
                     y=eRes$logLikelihoods)
    p1 <- ggplot(df, aes(x=x, y=y))
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_point()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("Log Likelihood")
    p1 <- ggplotly(p1)
    print(p1)

    fRes0 <- squareRootKF(A=A0, B=B0, C=C0, D=D0, x0=x00, initialStateType=initialStateType, SRSigmaX0=chol(x=V00), SRSigmaW=chol(x=Gamma0), SRSigmaV=chol(x=Sigma0), us=usTrain, zs=zsTrain)
    sRes0 <- smoothLDS(A=A0, mu=fRes0$x, V=fRes0$SigmaXHat, P=fRes0$SigmaX[2:length(fRes0$SigmaX)])

    fRes <- squareRootKF(A=eRes$A, B=B0, C=eRes$C, D=D0, x0=eRes$x0, initialStateType=initialStateType, SRSigmaX0=chol(x=eRes$V0), SRSigmaW=chol(x=eRes$Gamma), SRSigmaV=chol(x=eRes$Sigma), us=usTrain, zs=zsTrain)
    sRes <- smoothLDS(A=eRes$A, mu=fRes$x, V=fRes$SigmaXHat, P=fRes$SigmaX[2:length(fRes$SigmaX)])

    fResCrossValidated <- squareRootKF(A=eRes$A, B=B0, C=eRes$C, D=D0, x0=fRes$x[,ncol(fRes$x)], initialStateType=initialStateType, SRSigmaX0=chol(x=fRes$SigmaXHat[[length(fRes$SigmaXHat)]]), SRSigmaW=chol(x=eRes$Gamma), SRSigmaV=chol(x=eRes$Sigma), us=usTest, zs=zsTest)
    show(sprintf("cross-validated loglik=%.4f", sum(log(fResCrossValidated$c))))

    data <- data.frame()
    for(i in 1:nrow(trueLatentsTrain)) {
        dataBlock <- data.frame(sample=1:length(trueLatentsTrain[i,]),
                                latent=trueLatentsTrain[i,],
                                latentID=rep(i, length(trueLatentsTrain[i,])),
                                latentType=rep("true", length(trueLatentsTrain[i,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:nrow(sRes0$muHat)) {
        dataBlock <- data.frame(sample=1:length(sRes0$muHat[i,]),
                                latent=sRes0$muHat[i,],
                                latentID=rep(i, length(sRes0$muHat[i,])),
                                latentType=rep("initial", length(sRes0$muHat[i,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:nrow(sRes$muHat)) {
        dataBlock <- data.frame(sample=1:length(sRes$muHat[i,]),
                                latent=sRes$muHat[i,],
                                latentID=rep(i, length(sRes$muHat[i,])),
                                latentType=rep("estimated", length(sRes$muHat[i,])))
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

    p2 <- ggplotly(p2)
    print(p2)

#     data <- data.frame()
#     for(i in 1:nrow(trueLatentsTest)) {
#         dataBlock <- data.frame(sample=1:length(trueLatentsTest[i,]),
#                                 latentMean=trueLatentsTest[i,],
#                                 latentID=rep(i, length(trueLatentsTest[i,])),
#                                 latentType=rep("true", length(trueLatentsTest[i,])))
#         data <- rbind(data, dataBlock)
#     }
#     for(i in 1:nrow(fResCrossValidated$x)) {
#         dataBlock <- data.frame(sample=1:length(fResCrossValidated$x[i,]),
#                                 # wrong I should use the predicted, and not the
#                                 # smoothed mean
#                                 # latentMean=fResCrossValidated$x[i,],
#                                 # latentMean=fResCrossValidated$xtt1[i,],
#                                 latentID=rep(i, length(fResCrossValidated$x[i,])),
#                                 latentType=rep("predicted", length(fResCrossValidated$x[i,])))
#         data <- rbind(data, dataBlock)
#     }
#     nTestObs <- length(fResCrossValidated$SigmaXHat)
#     latentsSDs <- matrix(NA, nrow=nFactors, ncol=nTestObs)
#     for(n in 1:nTestObs) {
#         latentsSDs[,n] <- sqrt(diag(fResCrossValidated$SigmaXHat[[n]]))
#     }
#     dataErrorBars <- c()
#     for(i in 1:nrow(fResCrossValidated$x)) {
#         dataBlock <- data.frame(sample=1:length(fResCrossValidated$x[i,]),
#                                 latentMean=fResCrossValidated$x[i,],
#                                 latentSD=latentsSDs[i,],
#                                 latentID=rep(i, length(fResCrossValidated$x[i,])))
#         dataErrorBars <- rbind(dataErrorBars, dataBlock)
#     }
#     p3 <- ggplot()
#     p3 <- p3 + geom_line(data=data, mapping=aes(x=sample, y=latentMean, color=factor(latentID), linetype=factor(latentType)))
#     p3 <- p3 + scale_linetype_manual(values=c("twodash", "solid"))
#     p3 <- p3 + geom_ribbon(data=dataErrorBars, mapping=aes(x=sample, ymin=latentMean-1.96*latentSD, ymax=latentMean+1.96*latentSD, fill=factor(latentID)))
#     p3 <- p3 + geom_hline(yintercept=0)
#     p3 <- p3 + geom_vline(xintercept=0)
#     p3 <- p3 + ylab("Latent Value")
#     p3 <- p3 + xlab("Time")
#     p3 <- p3 + theme(legend.title = element_blank())
#     p3Plotly <- ggplotly(p3)
#     print(p3)
#     print(p3Plotly)

    # begin compute one-lag ahead observation predictions stats
#     ytt1 <- eRes$A%*%kfRes$xtt1+as.numeric(coef(kem)$A)
#     Wtt1 <- array(NA, dim=c(dimObs, dimObs, nObs))
#     if(kem$call$model$R=="diagonal and unequal") {
#         R <- diag(coef(kem)$R)
#     } else {
#         error("Functionality not yet implemented for R!=<diagonal and unequal>")
#     }
#     for(n in 1:nObs) {
#         Wtt1[,,n] <- Z%*%kfRes$Vtt1[,,n]%*%t(Z)+R
#     }
    # end compute one-lag ahead observation predictions stats


    browser()
}

processAll()
