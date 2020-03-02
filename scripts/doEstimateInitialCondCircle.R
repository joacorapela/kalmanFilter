
source("~/dev/learning/kalmanFilter/code/src/estimateKFInitialCond.R")

processAll <- function() {
    nFactors <- 2
    simulationFilename <- "results/simulationCircle.RData"

    simRes <- get(load(simulationFilename))
    zs <- simRes$x
    zsForFA <- t(as.matrix(zs))
    initialCond <- estimateKFInitialCond(z=zsForFA, nFactors=nFactors)

    browser()
}

processAll()

