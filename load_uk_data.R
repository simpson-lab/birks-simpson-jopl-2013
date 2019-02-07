## Packages required
library('analogue')
library('vegan')

## load UK data
uk <- read.cep("./data/uk-diatom-data.cep", force = TRUE)
ukph <- read.cep("./data/uk-ph-data.cep", force = TRUE)
nams <- rownames(ukph)
ukph <- ukph[,1]
names(ukph) <- nams
rm(nams)

## which samples are alread in SWAP?
ukSites <- rownames(uk)
swapSites <- rownames(swap)
ukTake <- ukSites %in% swapSites
ukSites[!ukTake]

## subset the ones we want
uk <- uk[!ukTake, ]
ukph <- ukph[!ukTake]

## drop species now not included
spp.drop <- which(colSums(uk) <= 0)
uk <- uk[, -spp.drop]

## now get rid of taxa that don't hit 1% abundance
uk <- chooseTaxa(uk, max.abun = 1, type = "AND")
