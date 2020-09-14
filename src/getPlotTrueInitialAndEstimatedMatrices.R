
require(plotly)
require(ggplot2)

getPlotTrueInitialAndEstimatedMatrices <- function(true=NA, initial=NA, estimated=NA, title="", xlab="Row Index", ylab="Value", trueLinetype="solid", initialLinetype="dot", estimatedLinetype="dash") {
    if(is.na(true[1]) && is.na(initial[1]) && is.na(estimated[1])) {
        error("At least one of true, initial or estimated must not be NA")
    }

    allData <- c()
    if(!is.na(true[1])) {
        trueColnames <- sprintf("true%d", 1:ncol(true))
        trueRownames <- sprintf("%d", 1:nrow(true))
        colnames(true) <- trueColnames
        rownames(true) <- trueRownames
        allData <- cbind(allData, true)
    }

    if(!is.na(initial[1])) {
        initialColnames <- sprintf("initial%d", 1:ncol(initial))
        initialRownames <- sprintf("%d", 1:nrow(initial))
        colnames(initial) <- initialColnames
        rownames(initial) <- initialRownames
        allData <- cbind(allData, initial)
    }

    if(!is.na(estimated[1])) {
        estimatedColnames <- sprintf("estimated%d", 1:ncol(estimated))
        estimatedRownames <- sprintf("%d", 1:nrow(estimated))
        colnames(estimated) <- estimatedColnames
        rownames(estimated) <- estimatedRownames
        allData <- cbind(allData, estimated)
    }
    allData <- data.frame(allData)

    # fig <- plot_ly(data=allMelted, x=~row, y=~value, linetype=~col, colour=~type, type='scatter', mode='markers')
    cols <- brewer.pal(max(3, ncol(true)), "Set1")
    fig <- plot_ly(data=allData, type='scatter', mode='lines+markers')
    if(!is.na(true[1])) {
        for(j in 1:ncol(true)) {
            fig <- fig%>%add_trace(x=1:nrow(allData), y=allData[[sprintf("true%d", j)]], name=sprintf("true[,%d]", j), line=list(color=cols[j], dash=trueLinetype), marker=list(color=cols[j]), type="scatter", mode="lines+markers")
        }
    }
    if(!is.na(initial[1])) {
        for(j in 1:ncol(true)) {
            fig <- fig%>%add_trace(x=1:nrow(allData), y=allData[[sprintf("initial%d", j)]], name=sprintf("initial[,%d]", j), line=list(color=cols[j], dash=initialLinetype), marker=list(color=cols[j]), type="scatter", mode="lines+markers")
        }
    }
    if(!is.na(estimated[1])) {
        for(j in 1:ncol(true)) {
            fig <- fig%>%add_trace(x=1:nrow(allData), y=allData[[sprintf("estimated%d", j)]], name=sprintf("estimated[,%d]", j), line=list(color=cols[j], dash=estimatedLinetype), marker=list(color=cols[j]), type="scatter", mode="lines+markers")
        }
    }
    fig <- fig %>% layout(title=title, xaxis=list(title=xlab, tickvals=1:nrow(allData)), yaxis=list(title=ylab))
    print(fig)
    browser()
    return(fig)
}

