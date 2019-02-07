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

## MLRC
library('rioja')

## do the internal CV
system.time({
N <- 10000
set.seed(2)
GLR.err <- t(sapply(seq_len(N),
                    function(x, spp, env, ...) fitGLR(spp, env, ...),
                    spp = swap / 100, env = ph,
                    n.cut = 5 ## makes no difference if use.glm is envoked
                    ))
})

saveRDS(GLR.err, file = "glr_internal.rds")

## predict pH at the UK sites ---------------------------------------##

## load the UK data now as we need to join
source("load_uk_data.R")

## select an optimisation set from the uk data that can be used in all
## methods that require it. Remaining samples for a 50-sample test set
##
## There are 73 UK sites, so select 50 for test and 23 for optimisation
##
## object to hold results
ukPred.glr.rmsep <- ukPred.glr.max.bias <- numeric(length = N)
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
    mod.glr <- MLRC(swap / 100, ph, use.glm = TRUE, max.iter = 100, n.cut = 5)
    ## do the predictions
    ukPred.glr <- predict(mod.glr, newdata = ukTestSet / 100)
    ukPred.glr.rmsep[i] <- sqrt(mean((ukTestSetEnv - ukPred.glr$fit)^2))
    ukPred.glr.max.bias[i] <-
        analogue:::maxBias(ukTestSetEnv - ukPred.glr$fit, ukTestSetEnv, n = 5)
}
## ------------------------------------------------------------------##

saveRDS(data.frame(rmsep = ukPred.glr.rmsep,
                   max.bias = ukPred.glr.max.bias),
        file = "ukPred_glr.rds")
