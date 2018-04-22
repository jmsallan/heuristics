setwd("~/Dropbox (UPC)/quanti_docs/FlowShop") #replace this with the path to the file in your computer

load("instances/TaillardFS.RData") #load instances
source("code/FSfunctions.R")       #load code with functions

test <- tai20.5[[1]]

makespan(test$tij, 1:20)
HillClimbing(test$tij, 1:20)

SA(test$tij, 1:20, 500, 10000)  #simulated annealing
TS(test$tij, 1:20, 1000)        #tabu search
GA(test$tij, 100, 100, 0.8)
