
require(ggplot2)
require(plotly)

plotTrueInitialAndEstimatedVectors <- function(true=NA, initial=NA, estimated=NA, title="", xlab="Index", ylab="Value", figFilename=NA) {
    if(is.na(true[1]) && is.na(initial[1]) && is.na(estimated[1])) {
        error("At least one of true, initial or estimated must not be NA")
    }

    if(!is.na(estimated[1])) {
        estimatedDF <- data.frame(index=1:length(estimated), value=estimated, type=rep("Estimated", length(estimated)))
    }
    if(!is.na(true[1])) {
        trueDF <- data.frame(index=1:length(true), value=true, type=rep("True", length(true)))
    }
    if(!is.na(initial[1])) {
        initialDF <- data.frame(index=1:length(initial), value=initial, type=rep("Initial", length(initial)))
    }

    if(!is.na(true[1])) {
        if(!is.na(initial[1])) {
            if(!is.na(estimated[1])) {
                allDF <- rbind(trueDF, initialDF, estimatedDF)
            } else {
                allDF <- rbind(trueDF, initialDF)
            }
        } else {
            if(!is.na(estimated[1])) {
                allDF <- rbind(trueDF, estimatedDF)
            } else {
                DF <- rbind(trueDF)
            }
        }
    } else {
        if(!is.na(initial[1])) {
            if(!is.na(estimated[1])) {
                allDF <- rbind(initialDF, estimatedDF)
            } else {
                allDF <- rbind(initialDF)
            }
        } else {
            if(!is.na(estimated[1])) {
                allDF <- rbind(estimatedDF)
            } else {
                error("At least one of true, initiall or estimated must not be NA")
            }
        }
    }

    p <- ggplot(data=allDF, aes(x=index, y=value, colour=type))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + ggtitle(title)
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)
    p <- ggplotly(p)
    if(!is.na(figFilename)) {
        htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(figFilename)), basename(figFilename)))
    }
    # print(p)
}

