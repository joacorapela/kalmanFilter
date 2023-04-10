
require(R.matlab)
require(plotly)
require(htmlwidgets)

processAll <- function() {
    mEstNumer <- 74925126
    rEstNumer <- 84389417

    mResultsFilenamePattern <- "/nfs/ghome/live/rapela/dev/research/gatsby-swc/gatsby/PLDS/matlabCode/results/%08d_PLDS_estimation.mat"
    rResultsFilenamePattern <- "results/%08d_estimation.RData"
    figFilenamePattern <- "figures/r%08d_m%08d_percentageOfExaplainedVariance.%s"

    mResultsFilename <- sprintf(mResultsFilenamePattern, mEstNumer)
    rResultsFilename <- sprintf(rResultsFilenamePattern, rEstNumer)
    rRes <- get(load(rResultsFilename))
    neuronIndices <- 1:length(rRes$pevs)
    rPEVs <- rRes$pevs

    mRes <- readMat(mResultsFilename)
    mPEVs <- mRes$pevs

    title <- sprintf("Mean Percentage of Explained Variance: GLDS=%.02f, PLDS=%.02f", mean(rPEVs), mean(mPEVs))

    fig <- plot_ly(type="scatter", mode="lines+markers")
    fig <- fig %>% add_trace(x=neuronIndices, y=rPEVs, name="GLDS")
    fig <- fig %>% add_trace(x=neuronIndices, y=mPEVs, name="PLDS")
    fig <- fig %>% layout(title=title, xaxis=list(title="Neuron Index"), yaxis=list(title="Percentage of Explained Variance"))
    # pngFilename <- sprintf(figFilenamePattern, rEstNumer, mEstNumer, "png")
    htmlFilename <- sprintf(figFilenamePattern, rEstNumer, mEstNumer, "html")
    # orca(p=fig, file=simPNGFilename)
    saveWidget(widget=fig, file=file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
    browser()
}

processAll()

rm(processAll)
