#---------------------------------------------------------------------#
#                                                                     #
# Script performing data analysis for JoPL paper in Rick SI           #
#                                                                     #
# Gavin L. Simpson                                                    #
# 13 April 2012                                                       #
#                                                                     #
#---------------------------------------------------------------------#

##--- load packages required to analyse the data --------------------##
require(ggplot2)
require(lme4)
require(multcomp)
##-------------------------------------------------------------------##

##--- Defaults ------------------------------------------------------##
OPTS <- opts(axis.text.x = theme_text(angle = 45, colour = "grey50",
             hjust = 1, vjust = 1))
##-------------------------------------------------------------------##

##--- Internal Test Set ---------------------------------------------##
## load the result vectors
WA <- readRDS("wa_internal.rds")
MAT <- readRDS("mat_internal.rds")
WAPLS <- readRDS("wapls_internal.rds")
GLR <- readRDS("glr_internal.rds")

## RMSEP: bind together
tfResults <- data.frame(WA = do.call("cbind.data.frame",
                        lapply(WA, `[`, i = , j = 1)),
                        WAPLS = WAPLS[,"rmsep"],
                        MAT = MAT[, "rmsep"],
                        GLR = GLR[, "rmsep"])

## stack the results
stackRes <- stack(tfResults)
names(stackRes) <- c("Error","Method")
stackRes <- stackRes[, c(2,1)]
stackRes <- within(stackRes,
                   Method <- factor(Method, levels = names(tfResults)))
## add a factor for "experiment" whic we use for a ranodm effect
stackRes <- within(stackRes,
                   Run <- rep(seq_len(nrow(tfResults)),
                              times = ncol(tfResults)))

## boxplot of the internal error
bw.internal <- ggplot(stackRes, aes(x = Method, y = Error)) +
    geom_boxplot() +
    ylab(expression(RMSEP[test])) +
    OPTS
bw.internal

## Maximum Bias: bind together
tfMaxBias <- data.frame(WA = do.call("cbind.data.frame",
                        lapply(WA, `[`, i = , j = 2)),
                        WAPLS = WAPLS[,"max.bias"],
                        MAT = MAT[, "max.bias"],
                        GLR = GLR[, "max.bias"])

## stack the results
stackBias <- stack(tfMaxBias)
names(stackBias) <- c("MaxBias","Method")
stackBias <- stackBias[, c(2,1)]
stackBias <- within(stackBias,
                    Method <- factor(Method, levels = names(tfMaxBias)))
## add a factor for "experiment" whic we use for a ranodm effect
stackBias <- within(stackBias,
                    Run <- rep(seq_len(nrow(tfMaxBias)),
                               times = ncol(tfMaxBias)))

## boxplot of the internal error
bw.internal.bias <- ggplot(stackBias, aes(x = Method, y = MaxBias)) +
    geom_boxplot() +
    ylab(expression(Maximum~Bias[test])) +
    OPTS
bw.internal.bias

##-------------------------------------------------------------------##

##--- UK Test Set ---------------------------------------------------##
## load the result vectors
uk.WA <- readRDS("ukPred_wa.rds")
uk.MAT <- readRDS("ukPred_mat.rds")
uk.WAPLS <- readRDS("ukPred_wapls.rds")
uk.GLR <- readRDS("ukPred_glr.rds")

## RMSEP: bind together
ukResults <- data.frame(WA = do.call("cbind.data.frame",
                        lapply(uk.WA, `[`, i = , j = 1)),
                        WAPLS = uk.WAPLS[,"rmsep"],
                        MAT = uk.MAT[,"rmsep"],
                        GLR = uk.GLR[,"rmsep"])

## stack the results
uk.stackRes <- stack(ukResults)
names(uk.stackRes) <- c("Error","Method")
uk.stackRes <- uk.stackRes[, c(2,1)]
uk.stackRes <- within(uk.stackRes,
                      Method <- factor(Method, levels = names(ukResults)))
## add a factor for "experiment" whic we use for a ranodm effect
uk.stackRes <- within(uk.stackRes,
                      Run <- rep(seq_len(nrow(ukResults)),
                                 times = ncol(ukResults)))

## boxplot of the internal error
bw.uk <- ggplot(uk.stackRes, aes(x = Method, y = Error)) +
    geom_boxplot() +
    ylab(expression(RMSEP[test])) +
    OPTS
bw.uk

## Max Bias: bind together
ukMaxBias <- data.frame(WA = do.call("cbind.data.frame",
                        lapply(uk.WA, `[`, i = , j = 2)),
                        WAPLS = uk.WAPLS[,"max.bias"],
                        MAT = uk.MAT[,"max.bias"],
                        GLR = uk.GLR[,"max.bias"])
## stack the results
uk.stackBias <- stack(ukResults)
names(uk.stackBias) <- c("MaxBias","Method")
uk.stackBias <- uk.stackBias[, c(2,1)]
uk.stackBias <- within(uk.stackBias,
                       Method <- factor(Method, levels = names(ukResults)))
## add a factor for "experiment" whic we use for a ranodm effect
uk.stackBias <- within(uk.stackBias,
                       Run <- rep(seq_len(nrow(ukResults)),
                                  times = ncol(ukResults)))

## boxplot of the internal error
bw.uk.bias <- ggplot(uk.stackBias, aes(x = Method, y = MaxBias)) +
    geom_boxplot() +
    ylab(expression(Maximum~Bias[test])) +
    OPTS
bw.uk.bias
##-------------------------------------------------------------------##

##--- Save the Plots ------------------------------------------------##
ggsave("bw_plot_rmsep_internal.pdf", plot = bw.internal)
ggsave("bw_plot_rmsep_uk_test.pdf", plot = bw.uk)
ggsave("bw_plot_maxbias_internal.pdf", plot = bw.internal.bias)
ggsave("bw_plot_maxbias_uk_test.pdf", plot = bw.uk.bias)
ggsave("figure1.png", plot = bw.internal)
ggsave("figure2.png", plot = bw.uk)
ggsave("figure3.png", plot = bw.internal.bias)
ggsave("figure4.png", plot = bw.uk.bias)
ggsave("figure1.eps", plot = bw.internal)
ggsave("figure2.eps", plot = bw.uk)
ggsave("figure3.eps", plot = bw.internal.bias)
ggsave("figure4.eps", plot = bw.uk.bias)
##-------------------------------------------------------------------##

##--- Test if the methods have different performances ---------------##
##
## RMSEP
## Internal
mod.int0 <- lmer(Error ~ 1 + (1 | Run), data = stackRes)
mod.int <- lmer(Error ~ Method + (1 | Run), data = stackRes)
anova(mod.int0, mod.int)
summary(mod.int)
tukey.int <- glht(mod.int, linfct = mcp(Method = "Tukey"))
summary(tukey.int)

## UK test
mod.uk0 <- lmer(Error ~ 1 + (1 | Run), data = uk.stackRes)
mod.uk <- lmer(Error ~ Method + (1 | Run), data = uk.stackRes)
anova(mod.uk0, mod.uk)
summary(mod.uk)
tukey.uk <- glht(mod.uk, linfct = mcp(Method = "Tukey"))
summary(tukey.uk)

## MaxBias
## Internal
mod.int0.b <- lmer(MaxBias ~ 1 + (1 | Run), data = stackBias)
mod.int.b <- lmer(MaxBias ~ Method + (1 | Run), data = stackBias)
anova(mod.int0.b, mod.int.b)
summary(mod.int.b)
tukey.int.b <- glht(mod.int.b, linfct = mcp(Method = "Tukey"))
summary(tukey.int.b)

## UK test
mod.uk0.b <- lmer(MaxBias ~ 1 + (1 | Run), data = uk.stackBias)
mod.uk.b <- lmer(MaxBias ~ Method + (1 | Run), data = uk.stackBias)
anova(mod.uk0.b, mod.uk.b)
summary(mod.uk.b)
tukey.uk.b <- glht(mod.uk.b, linfct = mcp(Method = "Tukey"))
summary(tukey.uk.b)
##-------------------------------------------------------------------##

##--- Pairwise Comparison Plots -------------------------------------##
pdf("tukey_internal_rmsep.pdf", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, version = "1.4")
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.int)
par(op)
dev.off()
pdf("tukey_uk_rmsep.pdf", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, version = "1.4")
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.uk)
par(op)
dev.off()
pdf("tukey_internal_maxbias.pdf", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, version = "1.4")
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.int.b)
par(op)
dev.off()
pdf("tukey_uk_maxbias.pdf", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, version = "1.4")
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.uk.b)
par(op)
dev.off()
##-------------------------------------------------------------------##

##--- Pairwise Comparison Plots - EPS -------------------------------##
postscript("figure5.eps", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, paper = "special", horizontal = FALSE)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.int)
par(op)
dev.off()
postscript("figure6.eps", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, paper = "special", horizontal = FALSE)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.uk)
par(op)
dev.off()
postscript("figure7.eps", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, paper = "special", horizontal = FALSE)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.int.b)
par(op)
dev.off()
postscript("figure8.eps", width = 10, height = 10, pointsize = 14,
    onefile = FALSE, paper = "special", horizontal = FALSE)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.uk.b)
par(op)
dev.off()
##-------------------------------------------------------------------##

##--- Pairwise Comparison Plots - PNGs ------------------------------##
png("figure5.png", width = 10, height = 10, pointsize = 14,
    units = "in", res = 75)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.int)
par(op)
dev.off()
png("figure6.png", width = 10, height = 10, pointsize = 14,
    units = "in", res = 75)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.uk)
par(op)
dev.off()
png("figure7.png", width = 10, height = 10, pointsize = 14,
    units = "in", res = 75)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.int.b)
par(op)
dev.off()
png("figure8.png", width = 10, height = 10, pointsize = 14,
    units = "in", res = 75)
op <- par(mar = c(5,15,4,2) + 0.1)
plot(tukey.uk.b)
par(op)
dev.off()
##-------------------------------------------------------------------##

##--- Write out tables of results -----------------------------------##
tab.res <- data.frame(Mean.RMSEP      = colMeans(tfResults),
                      SD.RMSEP        = sapply(tfResults, sd),
                      Mean.MaxBias      = colMeans(tfMaxBias),
                      SD.MaxBias        = sapply(tfMaxBias, sd),
                      UK.Mean.RMSEP = colMeans(ukResults),
                      UK.SD.RMSEP     = sapply(ukResults, sd),
                      UK.Mean.MaxBias = colMeans(ukMaxBias),
                      UK.SD.MaxBias     = sapply(ukMaxBias, sd))
write.csv(round(tab.res, 3), file = "tabulated_results.csv")
##-------------------------------------------------------------------##

##--- Number of analogues used --------------------------------------##
ana <- data.frame(SWAP = MAT[,3], UK = uk.MAT[,3])
stackAna <- stack(ana)
names(stackAna) <- c("k", "TestSet")
## build plot
ana.plt <- ggplot(stackAna, aes(x = k)) +
    geom_histogram(binwidth = 1) + facet_wrap( ~ TestSet) +
    ylab("Number of Runs") +
    xlab(expression(italic(k) ~~ "Number of analogues")) +
    scale_x_continuous(limits = c(1, 40))
ana.plt
## save
ggsave("number_of_analogues_hist.pdf", ana.plt,
       width = 7, height = 4.7, units = "in")
##-------------------------------------------------------------------##

##--- Number of components used -------------------------------------##
comp <- data.frame(SWAP = WAPLS[,3], UK = uk.WAPLS[,3])
stackComp <- stack(comp)
names(stackComp) <- c("ncomp", "TestSet")
## build plot
comp.plt <- ggplot(stackComp, aes(x = ncomp)) +
    geom_histogram(binwidth = 1) + facet_wrap( ~ TestSet) +
    ylab("Number of Runs") +
    xlab("Number of components") +
    scale_x_continuous(limits = c(1, 6))
comp.plt
## save
ggsave("number_of_components_hist.pdf", comp.plt,
       width = 7, height = 4.7, units = "in")
##-------------------------------------------------------------------##
