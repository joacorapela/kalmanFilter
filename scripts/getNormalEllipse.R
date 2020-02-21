getNormalEllipse <- function(mu, covar, nPoints=100, criticalValue=.95) {
    # from https://www.r-bloggers.com/drawing-a-95-confidence-interval-in-r/
    # see also https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

    eigenRes <- eigen(covar)
    evals <- eigenRes$values
    evecs <- eigenRes$vectors
    largestEval <- evals[1]
    largestEvec <- evecs[,1]
    smallestEval <- evals[2]
    smallestEvec <- evecs[,2]

    circleGrid <- seq(0, 2*pi, len=nPoints)

    # get major and minor elipse axis
    chisquareVal <- sqrt(qchisq(criticalValue, 2))
    a <- chisquareVal*sqrt(largestEval)
    b <- chisquareVal*sqrt(smallestEval)
    ellipseX <- a*cos(circleGrid)
    ellipseY <- b*sin(circleGrid)

    # rotate elipse to largest eigenvector
    rotationAngle <- atan2(largestEvec[2], largestEvec[1])
    rotationMatrix <- matrix(c(cos(rotationAngle), -sin(rotationAngle), 
                               sin(rotationAngle),  cos(rotationAngle)),
                             nrow=2)
    rEllipse <- cbind(ellipseX, ellipseY)%*%rotationMatrix

    # translate the ellipse to mean
    answer <- cbind(rEllipse[,1]+mu[1], rEllipse[,2]+mu[2])

    return(answer)
}
