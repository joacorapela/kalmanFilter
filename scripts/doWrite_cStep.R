
processAll <- function() {
    durationSamples <- 2000
    baselineValue <- 0.0
    stepValue <- 1.0
    stepStart <- 750
    stepEnd <- 1250
    filenamePattern <- "data/c_step_durationSamples%d_baselineValue%.02f_stepValue%.02f_stepStart%d_stepEnd%d.csv"

    filename <- sprintf(filenamePattern, durationSamples, baselineValue, stepValue, stepStart, stepEnd)
    c <- matrix(rep(baselineValue, times=durationSamples), nrow=1)
    c[1,stepStart:stepEnd] <- stepValue
    write.table(c, file=filename, row.names=FALSE, quote=FALSE, col.names=FALSE)
}

processAll()
rm(processAll)

