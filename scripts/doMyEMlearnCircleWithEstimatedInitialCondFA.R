
require(MASS)
require(reshape2)
require(plotly)
source("../src/emEstimationKF_R.R")
source("../src/filterLDS_R.R")
source("../src/smoothLDS_R.R")
source("../src/squareRootKF.R")
source("../src/smoothLDS_B.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/estimateKFInitialCondPPCA.R")
source("../src/plotTrueInitialAndEstimatedMatrices.R")
source("../src/plotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    dimObs <- 2*1
    nFactors <- 2
    tol <- 1e-6
    V0Var <- 1e1
    gamma0Var <- 1e1
    sigma0Var <- 1e1
    maxIter <- 100
    stateType0 <- "init00"
    # simulationFilename <- "results/simulationCircle.RData"
    simulationFilenamePattern <- "results/simulationCircleDObs%02d.RData"
    resultsFilnamePattern <- "results/circleDObs%02dMyEM_estimation.RData"

    simulationFilename <- sprintf(simulationFilenamePattern, dimObs)
    resultsFilename <- sprintf(resultsFilnamePattern, dimObs)
    simRes <- get(load(simulationFilename))
    dimLat <- ncol(simRes$A)
    # dimObs <- nrow(simRes$C)
    zs <- simRes$x

    # initialConds <- estimateKFInitialCondFA(z=t(as.matrix(zs)), nFactors=nFactors)
    initialConds <- estimateKFInitialCondPPCA(z=t(as.matrix(zs)), nFactors=nFactors)

    A0 <- initialConds$A
    Gamma0 <- gamma0Var*diag(rep(1, ncol(A0)))
    C0 <- initialConds$C
    # C0 <- simRes$C
    Sigma0 <- sigma0Var*diag(rep(1, times=dimObs))
    x0 <- simRes$mu0
    V0 <- V0Var*diag(rep(1, times=nFactors))

    res_emEstimationKF_R <- emEstimationKF_R(ys=zs, B0=initialConds$A, Q0=Gamma0, Z0=initialConds$C, R0=Sigma0, x00=x0, V00=V0, maxNIter=maxIter, tol=tol, stateType0=stateType0, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE), checkSymmetry=FALSE, checkPD=FALSE)

    save(file=resultsFilename, res_emEstimationKF_R)

    logLikFigFilename <- "figures/circleMmyEM_logLik.html"
    logLik <- res_emEstimationKF_R$logLike
    df <- data.frame(x=1:length(logLik), y=logLik)
    p <- ggplot(data=df, mapping=aes(x=x, y=y))
    p <- p + geom_point() + geom_line()
    p <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(logLikFigFilename)), basename(logLikFigFilename)))
    print(p)

    AFigFilename <- "figures//circleMyEM_A.html"
    plotTrueInitialAndEstimatedMatrices(true=simRes$A, initial=initialConds$A, estimated=res_emEstimationKF_R$B, title="B", figFilename=AFigFilename)

    CFigFilename <- "figures//circleMyEM_C.html"
    plotTrueInitialAndEstimatedMatrices(true=simRes$C, initial=initialConds$C, estimated=res_emEstimationKF_R$Z, title="Z", figFilename=CFigFilename)

    GammaFigFilename <- "figures//circleMyEM_varGamma.html"
    plotTrueInitialAndEstimatedMatrices(true=simRes$Gamma, estimated=res_emEstimationKF_R$Q, title="Gamma", figFilename=GammaFigFilename)

    SigmaFigFilename <- "figures//circleMyEM_Sigma.html"
    plotTrueInitialAndEstimatedMatrices(true=simRes$Sigma, estimated=res_emEstimationKF_R$R, title="Sigma", figFilename=SigmaFigFilename)

    x0FigFilename <- "figures//circleMyEM_x0.html"
    plotTrueInitialAndEstimatedVectors(true=simRes$mu0, estimated=res_emEstimationKF_R$x0, title="x0", figFilename=x0FigFilename)

    V0FigFilename <- "figures//circleMyEM_varV0.html"
    plotTrueInitialAndEstimatedMatrices(estimated=res_emEstimationKF_R$V0, title="V0", figFilename=V0FigFilename)

    # start run KF and KS using with initial parameters
    checkSymmetry <- FALSE
    symmetryTol <- 1e-4
    checkPD <- FALSE
    pdTol <- 1e-4
    fRes0 <- filterLDS_R(B=A0, Z=C0, x0=x0, V0=V0, stateType0=stateType0, Q=Gamma0, R=Sigma0, ys=zs, checkSymmetry=checkSymmetry, symmetryTol=symmetryTol, checkPD=checkPD, pdTol=pdTol)
    sRes0 <- smoothLDS_R(B=A0, xnn=fRes0$xnn, Vnn=fRes0$Vnn, xnn1=fRes0$xnn1, Vnn1=fRes0$Vnn1, stateType0=stateType0, x00=x0, V00=V0)
    # end run KF and KS using with initial conditions

    # start run KF and KS using with estimated parameters
    fRes <- filterLDS_R(B=res_emEstimationKF_R$B, Z=res_emEstimationKF_R$Z, x0=res_emEstimationKF_R$x0, V0=res_emEstimationKF_R$V0, stateType0=stateType0, Q=res_emEstimationKF_R$Q, R=res_emEstimationKF_R$R, ys=zs, checkSymmetry=checkSymmetry, symmetryTol=symmetryTol, checkPD=checkPD, pdTol=pdTol)
    sRes <- smoothLDS_R(B=res_emEstimationKF_R$B, xnn=fRes$xnn, Vnn=fRes$Vnn, xnn1=fRes$xnn1, Vnn1=fRes$Vnn1, stateType0=stateType0, x00=res_emEstimationKF_R$x0, V00=res_emEstimationKF_R$V0)
    # end run KF and KS using with estimated parameters


    data <- data.frame()
    for(i in 1:nrow(simRes$z)) {
        dataBlock <- data.frame(sample=1:length(simRes$z[i,]),
                                latent=simRes$z[i,],
                                latentID=rep(i, length(simRes$z[i,])),
                                latentType=rep("true",
                                               length(simRes$z[i,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:dim(sRes$xnN)[1]) {
        dataBlock <- data.frame(sample=1:length(sRes$xnN[i,1,]),
                                latent=sRes$xnN[i,1,],
                                latentID=rep(i, length(sRes$xnN[i,1,])),
                                latentType=rep("estimated",
                                               length(sRes$xnN[i,1,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:dim(sRes0$xnN)[1]) {
        dataBlock <- data.frame(sample=1:length(sRes0$xnN[i,1,]),
                                latent=sRes0$xnN[i,1,],
                                latentID=rep(i, length(sRes0$xnN[i,1,])),
                                latentType=rep("initial",
                                               length(sRes0$xnN[i,1,])))
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
    pPlotly <- ggplotly(p)
    latentsFigFilename <- "figures//circleMyEM_Latents.html"
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(latentsFigFilename)), basename(latentsFigFilename)))
    print(pPlotly)

    browser()
}

processAll()
