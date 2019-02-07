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

## Weighted Averaging
library('analogue')

## Inverse deshrinking

## Fit the full model
mod.wainv <- wa(swap, ph, deshrink = "inverse")
mod.wainv.tol <- wa(swap, ph, deshrink = "inverse", tol.dw = TRUE,
                    min.tol = 0.1, small.tol = "fraction")
mod.wacla <- wa(swap, ph, deshrink = "classical")
mod.wacla.tol <- wa(swap, ph, deshrink = "classical", tol.dw = TRUE,
                    min.tol = 0.1, small.tol = "fraction")
mod.wamon <- wa(swap, ph, deshrink = "monotonic")
mod.wamon.tol <- wa(swap, ph, deshrink = "monotonic", tol.dw = TRUE,
                    min.tol = 0.1, small.tol = "fraction")

## Do the internal CV
N <- 10000
set.seed(2)
WA.inv.err <- t(sapply(seq_len(N),
                       function(x, spp, env, ...) fitWA(spp, env, ...),
                       spp = swap, env = ph))
set.seed(2)
WA.cla.err <- t(sapply(seq_len(N),
                       function(x, spp, env, ...) fitWA(spp, env, ...),
                       spp = swap, env = ph, deshrink = "classical"))
set.seed(2)
WA.mon.err <- t(sapply(seq_len(N),
                       function(x, spp, env, ...) fitWA(spp, env, ...),
                       spp = swap, env = ph, deshrink = "monotonic"))
set.seed(2)
WA.inv.tol.err <- t(sapply(seq_len(N),
                           function(x, spp, env, ...) fitWA(spp, env, ...),
                           spp = swap, env = ph, tol.dw = TRUE,
                           min.tol = 0.1, small.tol = "fraction"))
set.seed(2)
WA.cla.tol.err <- t(sapply(seq_len(N),
                           function(x, spp, env, ...) fitWA(spp, env, ...),
                           spp = swap, env = ph,
                           deshrink = "classical", tol.dw = TRUE,
                           min.tol = 0.1, small.tol = "fraction"))
set.seed(2)
WA.mon.tol.err <- t(sapply(seq_len(N),
                           function(x, spp, env, ...) fitWA(spp, env, ...),
                           spp = swap, env = ph,
                           deshrink = "monotonic", tol.dw = TRUE,
                           min.tol = 0.1, small.tol = "fraction"))

waErrors <- list(Inverse = WA.inv.err,
                 Classical = WA.cla.err,
                 Monotonic = WA.mon.err,
                 InverseTol = WA.inv.tol.err,
                 ClassicalTol = WA.cla.tol.err,
                 MonotonicTol = WA.mon.tol.err)

saveRDS(waErrors, file = "wa_internal.rds")

## predict pH at the UK sites ---------------------------------------##

## load the UK data and process
source("load_uk_data.R")

## select an optimisation set from the uk data that can be used in all
## methods that require it. Remaining samples for a 50-sample test set
##
## There are 73 UK sites, so select 50 for test and 23 for optimisation
##
## object to hold results
ukPred.inv.rmsep <- ukPred.inv.tol.rmsep <- ukPred.cla.rmsep <-
    ukPred.cla.tol.rmsep <- ukPred.mon.rmsep <- ukPred.mon.tol.rmsep <-
    numeric(length = N)
ukPred.inv.max.bias <- ukPred.inv.tol.max.bias <- ukPred.cla.max.bias <-
    ukPred.cla.tol.max.bias <- ukPred.mon.max.bias <- ukPred.mon.tol.max.bias <-
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
    ##
    ## do the predictions
    ## WA Inv
    ukPred.inv <- predict(mod.wainv, newdata = ukTestSet)
    ukPred.inv.rmsep[i] <- sqrt(mean((ukTestSetEnv - ukPred.inv$pred$pred)^2))
    ukPred.inv.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.inv$pred$pred,
                           ukTestSetEnv, n = 5)
    ## WA Inv TolDW
    ukPred.inv.tol <- predict(mod.wainv.tol, newdata = ukTestSet)
    ukPred.inv.tol.rmsep[i] <- sqrt(mean((ukTestSetEnv -
                                          ukPred.inv.tol$pred$pred)^2))
    ukPred.inv.tol.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.inv.tol$pred$pred,
                           ukTestSetEnv, n = 5)
    ## WA Cla
    ukPred.cla <- predict(mod.wacla, newdata = ukTestSet)
    ukPred.cla.rmsep[i] <- sqrt(mean((ukTestSetEnv - ukPred.cla$pred$pred)^2))
    ukPred.cla.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.cla$pred$pred,
                           ukTestSetEnv, n = 5)
    ## WA Cla TolDW
    ukPred.cla.tol <- predict(mod.wacla.tol, newdata = ukTestSet)
    ukPred.cla.tol.rmsep[i] <- sqrt(mean((ukTestSetEnv -
                                          ukPred.cla.tol$pred$pred)^2))
    ukPred.cla.tol.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.cla.tol$pred$pred,
                           ukTestSetEnv, n = 5)
    ## WA Monotonic
    ukPred.mon <- predict(mod.wamon, newdata = ukTestSet)
    ukPred.mon.rmsep[i] <- sqrt(mean((ukTestSetEnv -
                                      ukPred.mon$pred$pred)^2))
    ukPred.mon.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.mon$pred$pred,
                           ukTestSetEnv, n = 5)
    ## WA Monotonic TolDW
    ukPred.mon.tol <- predict(mod.wamon.tol, newdata = ukTestSet)
    ukPred.mon.tol.rmsep[i] <- sqrt(mean((ukTestSetEnv -
                                          ukPred.mon.tol$pred$pred)^2))
    ukPred.mon.tol.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.mon.tol$pred$pred,
                           ukTestSetEnv, n = 5)
}
ukErrors <- list(Inverse = data.frame(rmsep = ukPred.inv.rmsep,
                 max.bias = ukPred.inv.max.bias),
                 Classical = data.frame(rmsep = ukPred.cla.rmsep,
                 max.bias = ukPred.cla.max.bias),
                 Monotonic = data.frame(rmsep = ukPred.mon.rmsep,
                 max.bias = ukPred.mon.max.bias),
                 InverseTol = data.frame(rmsep = ukPred.inv.tol.rmsep,
                 max.bias = ukPred.inv.tol.max.bias),
                 ClassicalTol = data.frame(rmsep = ukPred.cla.tol.rmsep,
                 max.bias = ukPred.cla.tol.max.bias),
                 MonotonicTol = data.frame(rmsep = ukPred.mon.tol.rmsep,
                 max.bias = ukPred.mon.tol.max.bias))
## ------------------------------------------------------------------##

saveRDS(ukErrors, file = "ukPred_wa.rds")
