setwd("~/Documents/heuristics/FlowShop") #replace with the path of your computer

source("code/FSfunctions.R")
load("instances/TaillardFS.RData")

#---- testing a small random instance -----

set.seed(2020)
instance <- matrix(sample(10:90, 100, replace = TRUE), 5, 20)

makespan(instance, PalmerTrapezes(instance)$pal) #1273
makespan(instance, PalmerTrapezes(instance)$tra) #1294

#simulated annealing
set.seed(2020)
SA01 <- SAFS(instance, PalmerTrapezes(instance)$pal, mu=1000, eval = TRUE)
SA01$fit #1214
SA02 <- SAFS(instance, PalmerTrapezes(instance)$pal, mu=1000, op="insertion", eval = TRUE)
SA02$fit #1211

#tabu search
TS01 <- TSFS(instance, PalmerTrapezes(instance)$pal)
TS01$obj #1230
TS01b <- TSFS(instance, PalmerTrapezes(instance)$pal, early = TRUE, iter=20)
TS01b$obj #1230

TS02 <- TSFS(instance, PalmerTrapezes(instance)$pal, op="insertion")
TS02$obj #1211

TS02b <- TSFS(instance, PalmerTrapezes(instance)$pal, op="insertion", tabu.size=7, eval=TRUE)
plotTS(TS02b)


# GRASP algorithm
GRASP01 <- GRASPFS(instance) #1209

set.seed(2020)
GRASP02 <- GRASPFS(instance, opt="SA", Tmax=10000, mu=1000) #1209

GRASP03 <- GRASPFS(instance, opt="TS") #1209

ILS01 <- ILSFS(instance) #1209

set.seed(2020)
ILS02 <- ILSFS(instance, opt="SA", Tmax=10000, mu=1000) #1209

ILS03 <- ILSFS(instance, opt="TS") #1209

# Genetic algorithm
GA01 <- GAFS(instance, npop=10, iter=1000, pmut=1, verbose = TRUE) #1209

#---- Taillard Instances -----

set.seed(1313)

pa.tai20.5 <- sapply(tai20.5, function(x) makespan(x$tij, PalmerTrapezes(x$tij)$pal))
tr.tai20.5 <- sapply(tai20.5, function(x) makespan(x$tij, PalmerTrapezes(x$tij)$tra))

hcs.tai20.5 <- sapply(tai20.5, function(x) HCFS(x$tij, 1:20)$obj)
hci.tai20.5 <- sapply(tai20.5, function(x) HCFS(x$tij, inisol=PalmerTrapezes(x$tij)$tra, op="insertion")$obj)

sas.tai20.5 <- sapply(tai20.5, function(x) SAFS(x$tij, inisol=PalmerTrapezes(x$tij)$tra, T=10000, mu=1000)$obj)
sai.tai20.5 <- sapply(tai20.5, function(x) SAFS(x$tij, inisol=PalmerTrapezes(x$tij)$tra, T=10000, mu=1000, op="insertion")$obj)

tss.tai20.5 <- sapply(tai20.5, function(x) TSFS(x$tij, inisol=PalmerTrapezes(x$tij)$tra)$obj)
tsi.tai20.5 <- sapply(tai20.5, function(x) TSFS(x$tij, inisol=PalmerTrapezes(x$tij)$tra, op="insertion")$obj)


r.tai20.5 <- data.frame(pa=pa.tai20.5, tr=tr.tai20.5, hcs=hcs.tai20.5, hci=hci.tai20.5, sas=sas.tai20.5, sai=sai.tai20.5, tss=tss.tai20.5, tsi=tsi.tai20.5)

