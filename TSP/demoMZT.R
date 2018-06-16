#set wd!

source("code/functionsTSP2.R")

test <- SampleTSP(25)

sol <- MZT(test$distances, verbose = TRUE)
plotTSP(test$coordinates, sol$sol)

save.image("results/MZT.RData")
