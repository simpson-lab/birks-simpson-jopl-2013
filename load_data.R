#---------------------------------------------------------------------#
#                                                                     #
# Script loads and processes the SWAP 167 diatom pH data set          #
#                                                                     #
# Gavin L. Simpson                                                    #
# 13 April 2012                                                       #
#                                                                     #
#---------------------------------------------------------------------#

## load vegan for read.cep
library('vegan')

## data are in ./data
fpath <- "./data/"
swap <- read.cep(paste(fpath, "SWAP167.cep", sep = ""), force = TRUE)
ph <- read.cep(paste(fpath, "SWAP167.env", sep = ""), force = TRUE)[,1]

## fix up the rownames on swap and add names to ph
rnams <- rownames(swap)
rnams <- sub("^X", "", rnams)
rownames(swap) <- names(ph) <- rnams

## clean up
rm(rnams)
rm(fpath)
