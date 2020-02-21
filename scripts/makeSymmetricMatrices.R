
makeSymmetricMatrices <- function(ms) {
    sMs <- list()
    for(i in 1:length(ms)) {
        sMs <- append(sMs, list(makeSymmetricMatrix(m=ms[[i]])))
    }
    return(sMs)
}

makeSymmetricMatrix <- function(m) {
    for(i in 1:(nrow(m)-1)) {
        for(j in (i+1):ncol(m)) {
            m[i,j] = m[j,i] = (m[i,j]+m[j,i])/2
        }
    }
    return(m)
}
