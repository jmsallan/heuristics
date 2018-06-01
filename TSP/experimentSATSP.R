

source("code/functionsTSP2.R")

sample <- SampleTSP(50)

NN <- NearestNeighbour(sample$distances)

TS <- TSTSP2opt(sample$distances, iter=50, NN$sol, eval = TRUE)

plot(1:50, TS$evalfit, type="l", col="red")
lines(1:50, TS$evalbest, col="blue")

#---- experimenting with size of tabu list

sizes <- c(5, 10, 15)

TS.sizes <- lapply(sizes, function(x) TSTSP2opt(sample$distances, iter=100, tabu.size=x, NN$sol))

#simulated annealing: if n=50, each move of TS has 50*(50-3)/2=1175 evaluations of objective. So 50 runs of TS equal to 50*1175=58750 iterations of SA.

SA <- SATSP2opt(sample$distances, NN$sol, Tmax=58750, mu=1, eval=TRUE)

ExperimentSA <- function(instance, inisol, label, Tmax, muValues, runs){
  
  totalruns <- length(muValues)*runs
  
  results <- as.data.frame(matrix(numeric(totalruns*3), totalruns, 3))
  names(results) <- c("initialSol", "mu", "distance")
  sols <- list()
  count <- 1
  
  for(r in 1:runs){
    for(mu in muValues){
      test <- SATSP2opt(instance, inisol = inisol, Tmax=Tmax, mu=mu)
      results[count, 1] <- label
      results[count, 2] <- mu
      results[count, 3] <- test$fit
      sols[[count]] <- test$sol
      count <- count+1
    }
  }
  
  results$distance <- as.numeric(results$distance)

  return(list(results=results, sols=sols))
}

set.seed(1313)
resultsSANN <- ExperimentSA(sample$distances, NN$sol, label="NN", Tmax=58750, muValues = c(1, 10, 50), runs=10)

set.seed(1111)
solution <- sample(1:50, 50)
resultsSArandom <- ExperimentSA(sample$distances, solution, label="random", Tmax=58750, muValues = c(1, 10, 50), runs=10)

allresults <- rbind(resultsSANN$results, resultsSArandom$results)

library(dplyr)

dataSA <- tbl_df(allresults)

dataSA %>% group_by(initialSol) %>% summarise(max=max(distance), mean=mean(distance), min=min(distance))
dataSA %>% group_by(mu) %>% summarise(max=max(distance), mean=mean(distance), min=min(distance))

library(ggplot2)

ggplot(dataSA, aes(initialSol, distance)) + geom_boxplot()
ggplot(dataSA, aes(factor(mu), distance)) + geom_boxplot()

ggplot(dataSA, aes(factor(mu), distance)) + geom_boxplot() + facet_grid(. ~ initialSol)

