#---------------------------------------------------------------------#
#                                                                     #
# Script containing functions for JoPL paper in Rick SI               #
#                                                                     #
# Gavin L. Simpson                                                    #
# 13 April 2012                                                       #
#                                                                     #
#---------------------------------------------------------------------#

fitWA <- function(spp, env,
                  chunk = 10, nTest = 37, nOpti = 20,
                  deshrink = "inverse", tol.dw = FALSE,
                  fill = "random", maxit = 200,
                  ...) {
    ## function to fit a WA model generating an RMSEP from
    ## an independently sampled test set

    ## Split sample
    ## opti : optimisation test - not needed for all models, but
    ##        still need to extract it so training sets are correct
    ##        size
    ## test : test set
    opti <- splitSample(env, chunk = chunk, take = nOpti, maxit = maxit)
    test <- splitSample(env[-opti], chunk = chunk, take = nTest,
                        fill = fill, maxit = maxit)

    ## select the samples
    train <- spp[-c(opti, test), ]
    trainEnv <- env[-c(opti, test)]
    testSet <- spp[test, ]
    testSetEnv <- env[test]

    ## which species does this get rid of for the training set
    spp.drop <- which(colSums(train) <= 0)

    ## drop these species
    if(length(spp.drop) > 0) {
        train <- train[, -spp.drop]
    }

    ## fit the model
    mod <- wa(train, trainEnv, deshrink = deshrink, tol.dw = tol.dw,
              ...)

    ## predict on the test
    pred <- predict(mod, newdata = testSet)$pred$pred
    rmsep <- sqrt(mean((testSetEnv - pred)^2))
    max.bias <- analogue:::maxBias(testSetEnv - pred, testSetEnv)
    out <- c(rmsep, max.bias)
    names(out) <- c("rmsep", "max.bias")
    out
}

fitMAT <- function(spp, env,
                   chunk = 10, nTest = 37, nOpti = 20,
                   method = "SQchord", ##kmax = 30,
                   fill = "random",
                   ...) {
    ## function to fit a MAT model generating an RMSEP from
    ## an independently sampled test set

    ## Split sample
    ## opti : optimisation test - not needed for all models, but
    ##        still need to extract it so training sets are correct
    ##        size
    ## test : test set
    opti <- splitSample(env, chunk = chunk, take = nOpti)
    test <- splitSample(env[-opti], chunk = chunk, take = nTest,
                        fill = fill)

    ## select the samples
    train <- spp[-c(opti, test), ]
    trainEnv <- env[-c(opti, test)]
    testSet <- spp[test, ]
    testSetEnv <- env[test]
    optiSet <- spp[opti, ]
    optiSetEnv <- env[opti]

    ## fit the model
    mod <- mat(train, trainEnv, method = method, ##kmax = kmax,
               ...)

    ## choose k from the optimisation set
    predOpti <- predict(mod, newdata = optiSet)
    k.rmsep <- sqrt(colMeans((optiSetEnv -
                              t(predOpti$predictions$model$predicted))^2))
    K <- which.min(k.rmsep)

    ## predict on the test
    pred <- predict(mod, newdata = testSet,
                    k = K)$predictions$model$predicted[K,]
    rmsep <- sqrt(mean((testSetEnv - pred)^2))
    max.bias <- analogue:::maxBias(testSetEnv - pred, testSetEnv)
    out <- c(rmsep, max.bias, K)
    names(out) <- c("rmsep", "max.bias", "k")
    out
}

## wrapper for WAPLS
fitWAPLS <- function(spp, env,
                     chunk = 10, nTest = 37, nOpti = 20,
                     npls = 6, fill = "random",
                     ...) {
    ## function to fit a WAPLS model generating an RMSEP from
    ## an independently sampled test set

    ## Split sample
    ## opti : optimisation test - not needed for all models, but
    ##        still need to extract it so training sets are correct
    ##        size
    ## test : test set
    opti <- splitSample(env, chunk = chunk, take = nOpti)
    test <- splitSample(env[-opti], chunk = chunk, take = nTest,
                        fill = fill)

    ## select the samples
    train <- spp[-c(opti, test), ]
    trainEnv <- env[-c(opti, test)]
    testSet <- spp[test, ]
    testSetEnv <- env[test]
    optiSet <- spp[opti, ]
    optiSetEnv <- env[opti]

    ## drop species not in the train
    spp.drop <- which(colSums(train) <= 0)
    if(length(spp.drop) > 0) {
        train <- train[, -spp.drop]
    }

    ## fit the model
    mod <- WAPLS(train, trainEnv, npls = npls, ...)

    ## choose number of WAPLS components from the optimisation set
    predOpti <- predict(mod, newdata = optiSet)
    ncomp <- sqrt(colMeans((optiSetEnv - predOpti$fit))^2)
    ncomp <- which.min(ncomp)

    ## predict on the test
    pred <- predict(mod, newdata = testSet)$fit[, ncomp]
    max.bias <- analogue:::maxBias(testSetEnv - pred, testSetEnv)
    rmsep <- sqrt(mean((testSetEnv - pred)^2))
    out <- c(rmsep, max.bias, ncomp)
    names(out) <- c("rmsep", "max.bias", "ncomp")
    out
}

fitGLR <- function(spp, env,
                   chunk = 10, nTest = 37, nOpti = 20,
                   max.iter = 50, use.glm = TRUE, n.cut = 5,
                   fill = "random",
                   ...) {
    ## function to fit a WA model generating an RMSEP from
    ## an independently sampled test set

    ## Split sample
    ## opti : optimisation test - not needed for all models, but
    ##        still need to extract it so training sets are correct
    ##        size
    ## test : test set
    opti <- splitSample(env, chunk = chunk, take = nOpti)
    test <- splitSample(env[-opti], chunk = chunk, take = nTest,
                        fill = fill)

    ## select the samples
    train <- spp[-c(opti, test), ]
    trainEnv <- env[-c(opti, test)]
    testSet <- spp[test, ]
    testSetEnv <- env[test]
    optiSet <- spp[opti, ]
    optiSetEnv <- env[opti]

    ## drop species not in the train
    spp.drop <- which(colSums(train) <= 0)
    if(length(spp.drop) > 0) {
        train <- train[, -spp.drop]
    }

    ## fit the model
    mod <- MLRC(train, trainEnv, max.iter = max.iter,
                use.glm = use.glm, n.cut = n.cut)

    ## predict on the test
    pred <- drop(predict(mod, newdata = testSet)$fit)
    rmsep <- sqrt(mean((testSetEnv - pred)^2))
    max.bias <- analogue:::maxBias(testSetEnv - pred, testSetEnv, n = 5)
    out <- c(rmsep, max.bias)
    names(out) <- c("rmsep", "max.bias")
    out
}
