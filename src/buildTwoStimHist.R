buildTwoStimHist <- function(stim, memory) {
    buildInputsBlock <- function(stim, memory) {
        inputs <- c()
        N <- length(stim)
        for(i in 0:(memory)) {
            inputs <- rbind(inputs, c(rep(0, times=i), stim[1:(N-i)]))
        }
        return(inputs)
    }

    stim1Hist <- buildInputsBlock(stim=stim[1,], memory=memory)
    stim2Hist <- buildInputsBlock(stim=stim[2,], memory=memory)
    stimHist <- rbind(stim1Hist, stim2Hist)
    return(stimHist)
}
