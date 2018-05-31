#set wd!

source("code/functionsTSP2.R")

load("results/TSPcities.Rdata")

#----loading country instances-----

library(netgen)

#----att48. optimal tour has length 10628----

att48 <- importFromTSPlibFormat("instances/att48.tsp")
Datt48 <- att48$distance.matrix

#nearest neighbour
NNatt48 <- NearestNeighbour(Datt48)

#tabu search
TS01att48 <- TSTSP2opt(Datt48, NNatt48$sol, iter=40, tabu.size = 5, eval = TRUE)

pdf("results/tsatt48.pdf", height = 5, width = 5)
plot(1:40, TS01att48$evalbest, type="l", col="blue", xlab="", ylab="")
lines(1:40, TS01att48$evalfit, col="red")
dev.off()


#simulated annealing
set.seed(1111)
SA01att48 <- SATSP2opt(Datt48, NNatt48$sol, Tmax=10000, mu=1, eval = TRUE)

pdf("results/saatt48.pdf", height = 5, width = 5)
plot(1:10000, SA01att48$evaltest, type="l", col="green", xlab="", ylab="")
lines(1:10000, SA01att48$evalfit, col="red")
lines(1:10000, SA01att48$evalbest, col="blue")
dev.off()


#ILS
set.seed(1313)
ILS01att48 <- ILSTSTSP2opt(Datt48, NNatt48$sol, rounds=50, iter=50)
ILS02att48 <- ILSTSTSP2opt(Datt48, 1:48, rounds=50, iter=50)

#comparing ILS
pdf("results/ilsatt48.pdf", width = 5, height = 5)
plot(1:50, ILS02att48$evalbest, type = "l", col="red", xlab="", ylab="")
lines(1:50, ILS01att48$evalbest, col="blue")
legend("topright", c("random", "nn"), lty=c(1,1), col=c("red", "blue"))
dev.off()

#GRASP
set.seed(2222)
GRASP01att48 <- GRASP2opt(Datt48, rcl.size = 2, tries=50, iter=50)

pdf("results/pathatt48.pdf", width=10, height=5)
par(mfrow=c(1,2))
plotTSP(att48$coordinates, NNatt48$sol)
plotTSP(att48$coordinates, ILS01att48$sol)
dev.off()

#GAs
set.seed(1111)
GA01att48 <- GATSP(Datt48, npop=100, iter=200, pmut=1, alpha=0)
GA02att48 <- GATSP(Datt48, npop=100, iter=200, memetic=TRUE, alpha=0.1)

#----a280. optimum 2579----

a280 <- importFromTSPlibFormat("instances/a280.tsp")
Da280 <- a280$distance.matrix

#nearest neighbour
NNa280 <- NearestNeighbour(Da280)

#tabu search
TS01a280 <- TSTSP2opt(Da280, NNa280$sol, iter=40, tabu.size = 5, eval = TRUE)

plot(1:40, TS01a280$evalbest, type="l", col="blue", xlab="", ylab="")
lines(1:40, TS01a280$evalfit, col="red")


set.seed(1313)
SA01a280 <- SATSP2opt(Da280, NNa280$sol, Tmax=310000, mu=1, eval = FALSE)

#ILS
set.seed(1313)
ILS01a280 <- ILSTSTSP2opt(Da280, NNa280$sol, rounds=30, iter=40)
ILS02a280 <- ILSTSTSP2opt(Da280, 1:280, rounds=30, iter=40)

plot(1:30, ILS02a280$evalbest, type = "l", col="red", xlab="", ylab="")
lines(1:30, ILS01a280$evalbest, col="blue")

set.seed(1111)
GRASP01a280 <- GRASP2opt(Da280, rcl.size = 2, tries=30, iter=40)

opta280 <- getOptTour("instances/a280.opt.tour")

pdf("results/patha280.pdf", width=10, height=5)
par(mfrow=c(1,2))
plotTSP(a280$coordinates, ILS01a280$sol)
plotTSP(a280$coordinates, opta280)
dev.off()

#GAs
set.seed(1111)
GA01a280 <- GATSP(Da280, npop=100, iter=200, pmut=1, alpha=0)
GA02a280 <- GATSP(Da280, npop=100, iter=200, memetic=TRUE, alpha=0.1)

#-----qatar: optimal 9352 ----

qa <- importFromTSPlibFormat("instances/qa194.tsp")
Dqa <- qa$distance.matrix

#nearest neighbour
NNqa <- NearestNeighbour(Dqa)

#nearest neighbour
NNa280 <- NearestNeighbour(Da280)

#tabu search
TS01qa <- TSTSP2opt(Dqa, NNqa$sol, iter=40, tabu.size = 5, eval = TRUE)

#simulated annealing
set.seed(1313)
SA01qa <- SATSP2opt(Dqa, NNqa$sol, Tmax=310000, mu=1, eval = FALSE)

#ILS
set.seed(1313)
ILS01qa <- ILSTSTSP2opt(Dqa, NNqa$sol, rounds=30, iter=40)
ILS02qa <- ILSTSTSP2opt(Dqa, 1:194, rounds=30, iter=40)

set.seed(1111)
GRASP01qa <- GRASP2opt(Dqa, rcl.size = 2, tries=30, iter=40)

set.seed(1111)
GA01qa <- GATSP(Dqa, npop=100, iter=200, pmut=1, alpha=0)
GA02qa <- GATSP(Dqa, npop=100, iter=200, memetic=TRUE, alpha=0.1)


save.image("results/TSPcities.Rdata")
