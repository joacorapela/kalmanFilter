require(ggplot2)
require(ggforce)

processAll <- function() {
    p <- ggplot() + 
         geom_ellipse(aes(x0=c(0,10), y0=c(0,10), a=c(10, 3), b=c(3, 3), angle=c(pi/4, 3*pi/4), m1=c(1, 3))) + 
         coord_fixed()

    print(p)

    browser()
}

processAll()
