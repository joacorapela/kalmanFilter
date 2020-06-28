
require(plotly)
require(ggplot2)

plotTrueInitialAndEstimatedMatrices <- function(true=NA, initial=NA, estimated=NA, title="", xlab="Row Index", ylab="Value", figFilename=NA) {
    if(is.na(true[1]) && is.na(initial[1]) && is.na(estimated[1])) {
        error("At least one of true, initial or estimated must not be NA")
    }

    if(!is.na(true[1])) {
        trueColnames <- sprintf("%d", 1:ncol(true))
        trueRownames <- sprintf("%d", 1:nrow(true))
        colnames(true) <- trueColnames
        rownames(true) <- trueRownames
        trueMelted <- melt(true)
        trueMelted$type <- rep("True", nrow(true))
    }

    if(!is.na(initial[1])) {
        initialColnames <- sprintf("%d", 1:ncol(initial))
        initialRownames <- sprintf("%d", 1:nrow(initial))
        colnames(initial) <- initialColnames
        rownames(initial) <- initialRownames
        initialMelted <- melt(initial)
        initialMelted$type <- rep("Initial", nrow(initial))
    }

    if(!is.na(estimated[1])) {
        estimatedColnames <- sprintf("%d", 1:ncol(estimated))
        estimatedRownames <- sprintf("%d", 1:nrow(estimated))
        colnames(estimated) <- estimatedColnames
        rownames(estimated) <- estimatedRownames

        estimatedMelted <- melt(estimated)
        estimatedMelted$type <- rep("Estimated", nrow(estimatedMelted))
    }

    if(!is.na(true[1])) {
        if(!is.na(initial[1])) {
            if(!is.na(estimated[1])) {
                allMelted <- rbind(trueMelted, initialMelted, estimatedMelted)
            } else {
                allMelted <- rbind(trueMelted, initialMelted)
            }
        } else {
            if(!is.na(estimated[1])) {
                allMelted <- rbind(trueMelted, estimatedMelted)
            } else {
                allMelted <- rbind(trueMelted)
            }
        }
    } else {
        if(!is.na(initial[1])) {
            if(!is.na(estimated[1])) {
                allMelted <- rbind(initialMelted, estimatedMelted)
            } else {
                allMelted <- rbind(initialMelted)
            }
        } else {
            if(!is.na(estimated[1])) {
                allMelted <- rbind(estimatedMelted)
            } else {
                error("At least one of true, initiall or estimated must not be NA")
            }
        }
    }

    colnames(allMelted)[1] <- "row"
    colnames(allMelted)[2] <- "col"
    allMelted[,"row"] <- allMelted[,"row"]
    allMelted[,"col"] <- as.character(allMelted[,"col"])
    p <- ggplot(data=allMelted, aes(x=row, y=value, linetype=col, colour=type))
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

