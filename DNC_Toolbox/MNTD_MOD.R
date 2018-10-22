# First value compares NN to N, don't need the second value
MNTD_MOD <- function (samp, dis, abundance.weighted = FALSE){
    N <- dim(samp)[1]
    MNTD_MOD <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        sppOutSample <- names(samp[i, samp[i, ] < 1])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppOutSample]
            if (abundance.weighted) {
                mntds <- apply(sample.dis, 2, min, na.rm = TRUE)
                sample.weights <- samp[i, sppInSample]
                MNTD_MOD[i] <- weighted.mean(mntds, sample.weights)
            }
            else {
                MNTD_MOD[i] <- mean(apply(sample.dis, 2, min, na.rm = TRUE))
            }
        }
        else {
            MNTD_MOD[i] <- NA
        }
    }
    MNTD_MOD
}