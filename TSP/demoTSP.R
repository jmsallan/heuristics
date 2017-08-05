setwd("~/Dropbox (UPC)/quanti_docs/TSP")
source("functionsTSP.R")

#-----Creating an instance of size 20 for testing----

instance <- SampleTSP(20, seed=1313)

#----Testing different heuristics----

nn.instance <- NearestNeighbour(instance, 1)

hc.instance <- HillClimbing(1:20, instance)

set.seed(1414)
sa.instance <- SA(1:20, instance, 19000, 100, TRUE)

ts.instance <- TS(1:20, instance, 100, report.evol = TRUE)

set.seed(1212)
ga.instance <- GA(instance, npop=100, iterations=100, pmut=0.2, report.evol = TRUE)

#----Plotting the results----

pdf("instanceTSP01.pdf", width=10, height=15)
par(mfrow=c(3,1))

plot(1:19000, sa.instance$evol, type="n", xlab="iter", ylab="dist", ylim = c(350, 1400), main="evolution of simulated annealing")
lines(1:19000, sa.instance$evol, col="blue")

plot(1:100, ts.instance$evol, type="n", xlab="iter", ylab="dist", ylim = c(350, 1000), main="evolution of tabu search")
lines(1:100, ts.instance$evol, col="blue")

plot(1:100, ga.instance$evol, type="n", xlab="iter", ylab="dist", ylim = c(350, 1000), main="evolution of genetic algorithm")
lines(1:100, ga.instance$evol, col="blue")

dev.off()

#----Loading matrix of att48 instance----

data48 <- read.table(file="att48data.tsp", header=FALSE)
colnames(data48) <- c("node", "xcoord", "ycoord")

att48 <- DistanceMatrixAtt(data48)

nn.att48 <- NearestNeighbour(att48)

set.seed(1212)
sa.att48 <- SA(1:48, att48, 112800, 300)

ts.att48 <- TS(1:48, att48, 150)

set.seed(1414)
ga.att48 <- GA(att48, npop=200, iterations=1000, pmut=0.2)
