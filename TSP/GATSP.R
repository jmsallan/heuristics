
source("code/functionsTSP2.R")

Instance <- SampleTSP(30)
TestInstance <- Instance$distances

NN.TestInstance <- NearestNeighbour(TestInstance)

TS.TestInstance <- TSTSP2opt(TestInstance, NN.TestInstance$sol)

GA01.TestInstance <- GATSP(TestInstance, npop=100, iter = 500, pmut = 1)

GA02.TestInstance <- GATSP(TestInstance, npop=100, crOX=FALSE, iter = 500, pmut = 1)

GA03.TestInstance <- GATSP(TestInstance, npop=100, crOX=FALSE, iter = 100, pmut = 1, memetic = TRUE, verbose=TRUE, alpha=0.1)

