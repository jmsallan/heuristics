setwd("~/Dropbox (UPC)/quanti_docs/TSP")

source("code/functionsTSP.R")
source("code/functionsTSPedges.R")

#-----two small instances to test the code----

d1 <- matrix(c(0, 3, 4, 2, 7, 3, 0, 4, 6, 3, 4, 4, 0, 5, 8, 2, 6, 5, 0, 6, 7, 3, 8, 6, 0), 5, 5)
d2 <- matrix(c(0, 2, 3, 4, 5, 1, 0, 6, 8, 3, 2, 9, 0, 1, 2, 3, 2, 1, 0, 4, 4, 2, 9, 6, 0), 5, 5)

bandb.d1 <- BandBTSPdist(d1, TRUE)
dist.d1 <- HeuristicDistances(d1)
nn.d1 <- NearestNeighbour(d1)

bandb.d2 <- BandBTSPdist(d2, TRUE)
dist.d2 <- HeuristicDistances(d2)
nn.d2 <- NearestNeighbour(d2)

InstanceTSP10 <- SampleTSP(10)

bandb.TSP10 <- BandBTSPdist(InstanceTSP10, TRUE)
dist.TSP10 <- HeuristicDistances(InstanceTSP10)
nn.TSP10 <- NearestNeighbour(InstanceTSP10)

save.image("results/resultsdemoTSPedges.RData")
