
library(mvtnorm) # References rmvnorm()
library(ellipse) # References ellipse()

processAll <- function() {

    # Set the covariance matrix
    sigma2 <- matrix(c(5, 2, 2, 5), ncol=2)
 
    # Set the means
    mu <- c(5, 5)
 
    # Get the correlation matrix
    P <- cov2cor(sigma2)
    # P <- sigma2
 
    # Generate the data
    p <- rmvnorm(n=50, mean=mu, sigma=sqrt(sigma2))
 
    # Plot the data
    plot(p)
 
    # Plot the ellipse
    lines(ellipse(P, centre = c(5,5)) , col='red')

    browser()
}

processAll()
