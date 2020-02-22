
source("estimateKFInitialCond.R")

processAll <- function() {
    nFactors <- 4
    simulationFilename <- "results/simulation2DTrajectory.RData"
    simRes <- get(load(simulationFilename))
    res <- estimateKFInitialCond(z=t(simRes$x), nFactors=nFactors)

    browser()
}

processAll()
