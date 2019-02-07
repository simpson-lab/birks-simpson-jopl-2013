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

## Modern Analogue Technique
library('analogue')

## do the internal CV
N <- 10000
set.seed(2)
MAT.err <- t(sapply(seq_len(N),
                    function(x, spp, env, ...) fitMAT(spp, env, ...),
                    spp = swap, env = ph))

## save it out
saveRDS(MAT.err, file = "mat_internal.rds")

## predict pH at the UK sites ---------------------------------------##

## load the UK data now as we need to join
source("load_uk_data.R")

## select an optimisation set from the uk data that can be used in all
## methods that require it. Remaining samples for a 50-sample test set
##
## There are 73 UK sites, so select 50 for test and 23 for optimisation
##
## object to hold results
ukPred.mat.rmsep <- ukPred.mat.max.bias <- ukPred.mat.k <- numeric(length = N)

tmpDat <- join(swap, uk, split = TRUE, verbose = TRUE)
swap <- tmpDat[["swap"]]
uk <- tmpDat[["uk"]]
rm(tmpDat)
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
    mod.mat <- mat(swap, ph, method = "SQchord", kmax = 30)
    ## choose k from the optimisation set
    predOpti <- predict(mod.mat, newdata = ukOptiSet)
    k.rmsep <- sqrt(colMeans((ukOptiSetEnv -
                              t(predOpti$predictions$model$predicted))^2))
    K <- which.min(k.rmsep)
    ## do the predictions
    ukPred.mat <- predict(mod.mat, newdata = ukTestSet)
    tmp.pred <- ukPred.mat$predictions$model$predicted[K,]
    ukPred.mat.rmsep[i] <- sqrt(mean((ukTestSetEnv - tmp.pred)^2))
    ukPred.mat.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - tmp.pred, ukTestSetEnv, n = 5)
    ukPred.mat.k[i] <- K
}
## ------------------------------------------------------------------##

saveRDS(data.frame(rmsep = ukPred.mat.rmsep,
                   max.bias = ukPred.mat.max.bias,
                   k = ukPred.mat.k),
        file = "ukPred_mat.rds")
