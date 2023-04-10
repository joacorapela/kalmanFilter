
source("../src/computeAIC.R")

processAll <- function() {
    estResNumber <- 39522788
    estResFilenamePattern <- "results/%08d_estimation.RData"

    estFilename <- 
    estResFilename <- sprintf(estResFilenamePattern, estResNumber)
    estRes <- get(load(estResFilename))
    AIC <- computeAIC(dsSSM=estRes$dsSSM)
    estRes <- c(estRes, list(AIC=AIC))

    stateDim <- nrow(estRes$dsSSM$B)
    show(sprintf("stateDim=%d, crossValidatedLogLike=%f, AIC=%f", stateDim, estRes$crossValidatedLogLike, AIC))
    save(file=estResFilename, estRes)
    browser()
}

processAll()
rm(processAll)
