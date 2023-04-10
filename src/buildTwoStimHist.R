buildStimHist <- function(stim, memory) {
    buildInputsBlock <- function(stim, memory) {
        inputs <- matrix(stim, nrow=1)
        if(memory>0) {
            N <- length(stim)
            for(i in 1:(memory)) {
                inputs <- rbind(inputs, c(rep(0, times=i), stim[1:(N-i)]))
            }
        }
        return(inputs)
    }

    stimHist <- buildInputsBlock(stim=stim, memory=memory)
    return(stimHist)
}