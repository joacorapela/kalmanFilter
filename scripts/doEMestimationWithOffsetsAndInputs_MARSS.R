
# require(MASS)
require(MARSS)
# require(reshape2)
require(ini)
source("../src/estimateKFInitialCondFA.R")
source("../src/estimateKFInitialCondPPCA.R")

processAll <- function() {
    estConfigNumber <- 1
    simResNumber <- 58388369
    simConfigFilenamePattern <- "data/%08d_simulation_metaData.ini"
    simResMetaDataFilenamePattern <- "results/%08d_simulation.ini"
    simResFilenamePattern <- "results/%s_simulation.RData"

    estResFilenamePattern <- "results/%08d_estimation.RData"
    estConfigFilenamePattern <- "data/%08d_estimation_metaData.ini"
    estResMetaDataFilenamePattern <- "results/%08d_estimation_metaData.ini"
    estFilenamePattern <- "results/%08d_estimation.RData"

    simResMetaDataFilename <- sprintf(simResMetaDataFilenamePattern, simResNumber)
    simResMetaData <- read.ini(simResMetaDataFilename)
    simConfigNumber <- as.numeric(simResMetaData$simulation_info$simConfigNumber)
    simConfigFilename <- sprintf(simConfigFilenamePattern, simConfigNumber)
    simConfig <- read.ini(simConfigFilename)

    sRate <- as.double(simConfig$control_variables$sRate)
    dt <- 1/sRate
    c <- as.matrix(read.table(simConfig$state_variables$cFilename))
    c <- c*dt
    dim(c) <- c(1, 1, length(c))
    # d <- as.matrix(read.table(simConfig$measurements_variables$dFilename))
    # dim(d) <- c(1, 1, length(d))

    simResFilename <- sprintf(simResFilenamePattern, simResNumber)

    exit <- FALSE
    while(!exit) {
        estResNumber <- sample(1e8, 1)
        estResFilename <- sprintf(estResFilenamePattern, estResNumber)
        if(!file.exists(estResFilename)) {
            exit <- TRUE
        }
    }
    estResMetaDataFilename <- sprintf(estResMetaDataFilenamePattern, estResNumber)
    show(sprintf("Saving estimation results in: %s", estResFilename))
    show(sprintf("Saving estimation meta data in: %s", estResMetaDataFilename))

    estConfigFilename <- sprintf(estConfigFilenamePattern, estConfigNumber)
    estConfig <- read.ini(estConfigFilename)

    # EM convergence tolerance
    tol <- as.double(estConfig$control_variables$tol)
    maxIter <- as.numeric(estConfig$control_variables$maxIter)

    V0 <- eval(parse(text=estConfig$initial_values$V0))
    M <- nrow(V0)
    u0 <- eval(parse(text=estConfig$initial_values$u0))
    C0 <- eval(parse(text=estConfig$initial_values$C0))
    Q0 <- eval(parse(text=estConfig$initial_values$Q0))
    a0 <- eval(parse(text=estConfig$initial_values$a0))
    # D0 <- eval(parse(text=estConfig$initial_values$D0))
    R0 <- eval(parse(text=estConfig$initial_values$R0))
    maxIter <- 100

    simRes <- get(load(simResFilename))
    y <- simRes$y
    P <- nrow(y)
    if(tolower(estConfig$initial_values$m0)=="simulated") {
        m0 <- simRes$m0
    } else {
        if(tolower(estConfig$initial_values$m0)=="randomuniform") {
            m0Min <- as.double(estConfig$initial_value$m0Min)
            m0Max <- as.double(estConfig$initial_value$m0Max)
            m0 <- runif(n=M, min=m0Min, max=m0Max)
        } else {
            m0 <- eval(parse(text=estConfig$initial_values$m0))
        }
    }

    # begin create model
    B1  <- matrix(list("b11", "b21", "b12", "b22"), nrow=2)
    U1  <- "unconstrained"
    C1  <- "unconstrained"
    Q1  <- "unconstrained"
    Z1  <- matrix(list("z11", "z21", "z12", "z22"), nrow=2)
    A1  <- "unconstrained"
    # D1  <- "unconstrained"
    R1  <- "unconstrained"
    pi1 <- "unequal"
    V01 <- "unconstrained"

    # model.list <- list(B=B1, U=U1, C=C1, Q=Q1, Z=Z1, A=A1, D=D1, d=d, R=R1, x0=pi1, V0=V01)
    model.list <- list(B=B1, U=U1, C=C1, Q=Q1, Z=Z1, A=A1, R=R1, x0=pi1, V0=V01)
    # end create model

    V0 <- eval(parse(text=estConfig$initial_values$V0))
    u0 <- eval(parse(text=estConfig$initial_values$u0))
    C0 <- eval(parse(text=estConfig$initial_values$C0))
    Q0 <- eval(parse(text=estConfig$initial_values$Q0))
    a0 <- eval(parse(text=estConfig$initial_values$a0))
    # D0 <- eval(parse(text=estConfig$initial_values$D0))
    R0 <- eval(parse(text=estConfig$initial_values$R0))

    # initialConds <- estimateKFInitialCondFA(z=t(as.matrix(y)), nFactors=M)
    initialConds <- estimateKFInitialCondPPCA(z=t(as.matrix(y)), nFactors=M)
    # initialConds <- c(initialConds, list(u=u0, C=C0, Q=Q0, a=a0, D=D0, R=R0, m0=m0, V0=V0))
    initialConds <- c(initialConds, list(u=u0, C=C0, Q=Q0, a=a0, R=R0, m0=m0, V0=V0))

    B0 <- matrix(as.vector(initialConds$B), ncol=1)
    u0 <- matrix(as.vector(initialConds$u), ncol=1)
    C0 <- matrix(as.vector(initialConds$C), ncol=1)
    Q0 <- initialConds$Q
    Q0lowerTri <- matrix(Q0[lower.tri(Q0, diag=TRUE)], ncol=1)
    x0 <- matrix(as.vector(initialConds$m0), ncol=1)
    V0 <- initialConds$V0
    V0lowerTri <- matrix(V0[lower.tri(V0, diag=TRUE)], ncol=1)
    Z0 <- matrix(as.vector(initialConds$Z), ncol=1)
    a0 <- matrix(as.vector(initialConds$a), ncol=1)
    # D0 <- matrix(as.vector(initialConds$D), ncol=1)
    R0 <- initialConds$R
    R0lowerTri <- matrix(R0[lower.tri(R0, diag=TRUE)], ncol=1)
    inits <- list(B=B0, U=u0, C=C0, Q=Q0lowerTri, x0=x0, V0=V0lowerTri, Z=Z0, A=a0, R=R0lowerTri)
    # inits <- list(B=B0, Z=Z0)

    control <- list(maxit=maxIter, trace=1, safe=TRUE)

    kem <- MARSS(y, model=model.list, inits=inits, control=control, silent=2, fun.kf="MARSSkfss")
    # kem <- MARSS(y, model=model.list, control=control, silent=2, fun.kf="MARSSkfss")

    estRes <- c(kem, list(initialConds=initialConds))
    save(file=estResFilename, estRes)

    metaData <- list()
    metaData[["simulation_info"]] <- list(simResNumber=simResNumber)
    metaData[["estimation_config_info"]] <- list(estConfigNumber=estConfigNumber)
    write.ini(x=metaData, filepath=estResMetaDataFilename)

    browser()
}

processAll()
