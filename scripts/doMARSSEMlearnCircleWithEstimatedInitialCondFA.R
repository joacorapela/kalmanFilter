
require(MASS)
require(reshape2)
require(MARSS)
require(plotly)
source("../src/squareRootKF.R")
source("../src/smoothLDS_B.R")
source("../src/estimateKFInitialCondFA.R")
source("../src/plotTrueInitialAndEstimatedMatrices.R")
source("../src/plotTrueInitialAndEstimatedVectors.R")

processAll <- function() {
    nFactors <- 2
    maxIter <- 100
    simulationFilename <- "results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    dimLat <- ncol(simRes$A)
    dimObs <- nrow(simRes$C)
    zs <- simRes$x

    # begin create model
    B1  <- matrix(list("b11", "b21", "b12", "b22"), nrow=2)
    U1  <- "zero"
    # Q1  <- "diagonal and equal"
    Q1  <- "unconstrained"
    # Z1  <- factor(c(1,1,1,1,2,2,2,2))
    Z1  <- matrix(list("z11", "z21", "z31", "z41", "z51", "z61", "z71", "z81", "z12", "z22", "z32", "z42", "z52", "z62", "z72", "z82"), nrow=8)
    A1  <- "zero"
    # R1  <- "diagonal and unequal"
    R1  <- "unconstrained"
    pi1 <- "unequal"
    V01 <- "diagonal and equal"

    model.list <- list(B=B1, U=U1, Q=Q1, Z=Z1, A=A1, R=R1, x0=pi1, V0=V01)
    # end create model

    # begin set initial conditions
    zsForFA <- t(as.matrix(zs))
    initialConds <- estimateKFInitialCondFA(z=zsForFA, nFactors=nFactors)

    B0 <- matrix(as.vector(initialConds$A), ncol=1)
    Z0 <- matrix(as.vector(initialConds$C), ncol=1)
    R0 <- diag(initialConds$sigmaDiag)
    R0lowerTri <- matrix(R0[lower.tri(R0, diag=TRUE)], ncol=1)
    control <- list(maxit=maxIter, trace=1)

    inits <- list(B=B0, Z=Z0, R=R0lowerTri)
    # inits <- list(B=B0, Z=Z0)
    # end set initial conditions

    kem <- MARSS(zs, model=model.list, inits=inits, control=control, silent=2)

    if(kem$convergence>1) {
        warning(sprintf("MARSS did not converged (convergence=%d)", kem$convergence))
    }
    logLikFigFilename <- "figures/circleMARSS_logLik.html"
    logLik <- kem$iter.record$logLik
    df <- data.frame(x=1:length(logLik), y=logLik)
    p <- ggplot(data=df, mapping=aes(x=x, y=y))
    p <- p + geom_point() + geom_line()
    p <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(logLikFigFilename)), basename(logLikFigFilename)))
    print(p)

    AFigFilename <- "figures//circleMARSS_A.html"
    plotTrueInitialAndEstimatedMatrices(true=simRes$A, initial=initialConds$A, estimated=matrix(coef(kem)$B, nrow=dimLat), title="A", figFilename=AFigFilename)

    CFigFilename <- "figures//circleMARSS_C.html"
    plotTrueInitialAndEstimatedMatrices(true=simRes$C, initial=initialConds$C, estimated=matrix(coef(kem)$Z, nrow=dimObs), title="C", figFilename=CFigFilename)

    GammaFigFilename <- "figures//circleMARSS_varGamma.html"
    plotTrueInitialAndEstimatedVectors(true=simRes$Gamma[1,1], estimated=coef(kem)$Q, title="Gamma Variance", figFilename=GammaFigFilename)

    SigmaFigFilename <- "figures//circleMARSS_Sigma.html"
    plotTrueInitialAndEstimatedVectors(true=diag(simRes$Sigma), initial=initialConds$sigmaDiag, estimated=coef(kem)$R, title="Sigma Diagonal", figFilename=SigmaFigFilename)

    x0FigFilename <- "figures//circleMARSS_x0.html"
    plotTrueInitialAndEstimatedVectors(true=simRes$mu0, estimated=coef(kem)$x0, title="x0", figFilename=x0FigFilename)

    V0FigFilename <- "figures//circleMARSS_varV0.html"
    plotTrueInitialAndEstimatedVectors(true=simRes$V0[1,1], estimated=coef(kem)$V0, title="V0 Variance", figFilename=V0FigFilename)

    kfRes <- MARSSkf(kem)
    kem0 <- kem
    kem0$par <- kem0$start
    kfRes0 <- MARSSkf(kem0)

    data <- data.frame()
    for(i in 1:nrow(simRes$z)) {
        dataBlock <- data.frame(sample=1:length(simRes$z[i,]), latent=simRes$z[i,], latentID=rep(i, length(simRes$z[i,])), latentType=rep("true", length(simRes$z[i,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:nrow(kfRes$xtT)) {
        dataBlock <- data.frame(sample=1:length(kfRes$xtT[i,]),
                                latent=kfRes$xtT[i,],
                                latentID=rep(i, length(kfRes$xtT[i,])),
                                latentType=rep("estimated",
                                               length(kfRes$xtT[i,])))
        data <- rbind(data, dataBlock)
    }
    for(i in 1:nrow(kfRes0$xtT)) {
        dataBlock <- data.frame(sample=1:length(kfRes0$xtT[i,]),
                                latent=kfRes0$xtT[i,],
                                latentID=rep(i, length(kfRes0$xtT[i,])),
                                latentType=rep("initial",
                                               length(kfRes0$xtT[i,])))
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
    latentsFigFilename <- "figures//circleMARSS_Latents.html"
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(latentsFigFilename)), basename(latentsFigFilename)))
    print(p)

    browser()
}

processAll()
