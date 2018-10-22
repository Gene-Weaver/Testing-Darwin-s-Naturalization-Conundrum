# Modified 12-29-17 -> apply(X = mntd.rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
# The 1 looks along the invasives ROW to find min.
# The 2 (wrong) looks along columns, the value from N-NN still
SES_MNTD_MOD <- function (samp, dis, null.model = c("taxa.labels", "richness", 
                                    "frequency", "sample.pool", "phylogeny.pool", "independentswap", 
                                    "trialswap"), abundance.weighted = FALSE, runs = 999, iterations = 1000) 
{
    dis <- as.matrix(dis)
    mntd.obs <- MNTD_MOD(samp, dis, abundance.weighted)
    null.model <- match.arg(null.model)
    mntd.rand <- switch(null.model, taxa.labels = t(replicate(runs,
                                                              MNTD_MOD(samp, taxaShuffle(dis), abundance.weighted))), richness = t(replicate(runs,
                                                                                                                                             MNTD_MOD(randomizeMatrix(samp, null.model = "richness"), 
                                                                dis, abundance.weighted))),
                                                                frequency = t(replicate(runs, 
                                                                                        MNTD_MOD(randomizeMatrix(samp, null.model = "frequency"), 
                                                                dis, abundance.weighted))), sample.pool = t(replicate(runs, 
                                                                                                                      MNTD_MOD(randomizeMatrix(samp, null.model = "richness"), 
                                                                dis, abundance.weighted))), phylogeny.pool = t(replicate(runs, 
                                                                                                                         MNTD_MOD(randomizeMatrix(samp, null.model = "richness"), 
                                                                taxaShuffle(dis), abundance.weighted))), independentswap = t(replicate(runs, 
                                                                                                                                       MNTD_MOD(randomizeMatrix(samp, null.model = "independentswap", 
                                                                iterations), dis, abundance.weighted))), trialswap = t(replicate(runs, 
                                                                                                                                 MNTD_MOD(randomizeMatrix(samp, null.model = "trialswap", 
                                                                iterations), dis, abundance.weighted))))
    mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 2, FUN = mean, 
                            na.rm = TRUE)
    mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 2, FUN = sd, 
                          na.rm = TRUE)
    mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
    mntd.obs.rank <- apply(X = rbind(mntd.obs, mntd.rand), MARGIN = 2, 
                           FUN = rank)[1, ]
    mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)
    data.frame(ntaxa = specnumber(samp), mntd.obs, mntd.rand.mean, 
               mntd.rand.sd, mntd.obs.rank, mntd.obs.z, mntd.obs.p = mntd.obs.rank/(runs + 
                                                                                        1), runs = runs, row.names = row.names(samp))
}