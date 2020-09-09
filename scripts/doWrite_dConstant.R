
processAll <- function() {
    durationSamples <- 2000
    constantValue <- 0.0
    filenamePattern <- "data/d_constant_durationSamples%d_constantValue%.02f.csv"

    filename <- sprintf(filenamePattern, durationSamples, constantValue)
    d <- matrix(rep(constantValue, times=durationSamples), nrow=1)
    write.table(d, file=filename, row.names=FALSE, quote=FALSE, col.names=FALSE)
}

processAll()
rm(processAll)

