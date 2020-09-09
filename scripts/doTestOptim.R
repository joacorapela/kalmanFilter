myFunc <- function(a, b, x) {
    y <- exp(a*x)+b*x
    return(y)
}

generateData <- function(a, b, N, noiseSD=1e-4) {
    x <- rnorm(N)
    noise <- rnorm(N)*noiseSD
    y <- myFunc(a=a, b=b, x=x)
    answer <- cbind(y, x)
    return(answer)
}

processAll <- function() {
    aTrue <- 1
    sdNoise_a <- 1e-0
    bTrue <- .5
    sdNoise_b <- 1e+1
    N <- 1000

    data <- generateData(a=aTrue, b=bTrue, N=N)
    y <- data[,1]
    x <- data[,2]
    paramsVec0 <- c(aTrue+rnorm(1)*sdNoise_a, bTrue+rnorm(1)*sdNoise_b)
    f <- function(paramsVec) {
        aHat <- paramsVec[1]
        bHat <- paramsVec[2]

        yHat <- myFunc(a=aHat, b=bHat, x=x)
        error <- mean((y-yHat)^2)
        return(error)
    }
    est <- optim(par=paramsVec0, fn=f, gr=NULL, method='BFGS', hessian=TRUE, control=list(trace=1, REPORT=1))

    show("True params")
    show(c(aTrue, bTrue))
    show("Initial params")
    show(paramsVec0)
    show("Estimated params")
    show(est$par)

    browser()
}

processAll()

