source("code/functionsTSP2.R")

sample <- SampleTSP(30)
D <- sample$distances

nn.sol <- NearestNeighbour(D)

ts.sol <- TSTSP2opt(D, nn.sol$sol)

ga01 <- GATSP(D, npop=30, pmut=1, iter=100, alpha = 0.1, verbose = TRUE)
ga02 <- GATSP(D, npop=30, pmut=1, iter=10, alpha = 0.1, memetic=TRUE, verbose = TRUE)

plotTSP(sample$coordinates, nn.sol$sol)
plotTSP(sample$coordinates, ga01$sol)
plotTSP(sample$coordinates, ga02$sol)
