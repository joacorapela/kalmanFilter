
library(mvtnorm)
source("getNormalEllipse.R")

processAll <- function() {

    nDataPoints <- 100
    nEllipsePoints <- 100
    criticalValue <- .99

    covar <- matrix(c(5, -3, -3, 5), ncol=2)
 
    mu <- c(5, 5)

    data <- rmvnorm(n=nDataPoints, mean=mu, sigma=covar)
    ellipse <- getNormalEllipse(mu=mu, covar=covar, nPoints=nEllipsePoints, 
                                criticalValue=criticalValue)

    xlim <- range(c(data[,1], ellipse[,1]))
    ylim <- range(c(data[,2], ellipse[,2]))
    plot(data[,1], data[,2], type='p', col="black", xlim=xlim, ylim=ylim)
    lines(ellipse[,1], ellipse[,2], col="red")

    browser()
}

processAll()
