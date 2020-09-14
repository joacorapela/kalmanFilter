
require(ggplot2)
require(plotly)

getPlotTrueInitialAndEstimatedVectors <- function(true=NA, initial=NA, estimated=NA, title="", xlab="Index", ylab="Value", trueLinetype="solid", initialLinetype="dot", estimatedLinetype="dash") {
    if(is.na(true[1]) && is.na(initial[1]) && is.na(estimated[1])) {
        error("At least one of true, initial or estimated must not be NA")
    }

    allData <- c()
    if(!is.na(true[1])) {
        true <- matrix(true, ncol=1)
        colnames(true) <- "true"
        allData <- cbind(allData, true)
    }
    if(!is.na(initial[1])) {
        initial <- matrix(initial, ncol=1)
        colnames(initial) <- "initial"
        allData <- cbind(allData, initial)
    }
    if(!is.na(estimated[1])) {
        estimated <- matrix(estimated, ncol=1)
        colnames(estimated) <- "estimated"
        allData <- cbind(allData, estimated)
    }
    allData <- data.frame(allData)

    fig <- plot_ly(data=allData, type='scatter', mode='lines+markers')
    if(!is.na(true[1])) {
        fig <- fig%>%add_trace(x=1:nrow(allData), y=allData[["true"]], name="true", line=list(dash=trueLinetype), type="scatter", mode="lines+markers")
    }
    if(!is.na(initial[1])) {
        fig <- fig%>%add_trace(x=1:nrow(allData), y=allData[["initial"]], name="initial", line=list(dash=initialLinetype), type="scatter", mode="lines+markers")
    }
    if(!is.na(estimated[1])) {
        fig <- fig%>%add_trace(x=1:nrow(allData), y=allData[["estimated"]], name="estimated", line=list(dash=estimatedLinetype), type="scatter", mode="lines+markers")
    }
    fig <- fig %>% layout(title=title, xaxis=list(title=xlab, tickvals=1:nrow(allData)), yaxis=list(title=ylab))
    return(fig)
}

