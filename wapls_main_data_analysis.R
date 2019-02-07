#---------------------------------------------------------------------#
#                                                                     #
# Script performing data analysis for JoPL paper in Rick SI           #
#                                                                     #
# Gavin L. Simpson                                                    #
# 13 April 2012                                                       #
#                                                                     #
#---------------------------------------------------------------------#

## load functions for this paper
source("paper_fun.R")

## load the data
source("load_data.R")

## need analogue but load before rioja
library('analogue')

## WA-PLS
library('rioja')

## do the internal CV
N <- 10000
set.seed(2)
WAPLS.err <- t(sapply(seq_len(N),
                      function(x, spp, env, ...) fitWAPLS(spp, env, ...),
                      spp = swap, env = ph))

saveRDS(WAPLS.err, file = "wapls_internal.rds")

## predict pH at the UK sites ---------------------------------------##

## load the UK data now as we need to join
source("load_uk_data.R")

## select an optimisation set from the uk data that can be used in all
## methods that require it. Remaining samples for a 50-sample test set
##
## There are 73 UK sites, so select 50 for test and 23 for optimisation
##
## object to hold results
ukPred.wapls.rmsep <- ukPred.wapls.max.bias <- ukPred.wapls.ncomp <-
    numeric(length = N)
## seed
set.seed(25)
## loop
for (i in seq_len(N)) {
    ukOptiTake <- splitSample(ukph, chunk = 10, take = 23, fill = "random")
    ukTestSet <- uk[-ukOptiTake, ]
    ukOptiSet <- uk[ukOptiTake, ]
    ukTestSetEnv <- ukph[-ukOptiTake]
    ukOptiSetEnv <- ukph[ukOptiTake]
    ## Fit the full model
    mod.wapls <- WAPLS(swap, ph, npls = 6)
    ## choose ncomp from the optimisation set
    predOpti <- predict(mod.wapls, newdata = ukOptiSet)
    ncomp <- sqrt(colMeans((ukOptiSetEnv - predOpti$fit))^2)
    ncomp <- which.min(ncomp)
    ## do the predictions
    ukPred.wapls <- predict(mod.wapls, newdata = ukTestSet)
    tmp.pred <- ukPred.wapls$fit[, ncomp]
    ukPred.wapls.rmsep[i] <- sqrt(mean((ukTestSetEnv - tmp.pred)^2))
    ukPred.wapls.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - tmp.pred, ukTestSetEnv, n = 5)
    ukPred.wapls.ncomp[i] <- ncomp
}
##-------------------------------------------------------------------##

saveRDS(data.frame(rmsep = ukPred.wapls.rmsep,
                   max.bias = ukPred.wapls.max.bias,
                   ncomp = ukPred.wapls.ncomp),
        file = "ukPred_wapls.rds")
