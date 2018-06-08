

source("code/functionsTSP2.R")

sample <- SampleTSP(50)

NN <- NearestNeighbour(sample$distances)

t1 <- Sys.time()
TS <- TSTSP2opt(sample$distances, iter=50, NN$sol, eval = TRUE)
t2 <- Sys.time()

time.ts1 <- t2 - t1

plot(1:50, TS$evalfit, type="l", col="red")
lines(1:50, TS$evalbest, col="blue")

#---- experimenting with size of tabu list

sizes <- c(5, 10, 15)

TS.sizes <- lapply(sizes, function(x) TSTSP2opt(sample$distances, iter=100, tabu.size=x, NN$sol))

#simulated annealing: if n=50, each move of TS has 50*(50-3)/2=1175 evaluations of objective. So 50 runs of TS equal to 50*1175=58750 iterations of SA.

t1 <- Sys.time()
set.seed(1111)
SA <- SATSP2opt(sample$distances, NN$sol, Tmax=50000, mu=1, eval=TRUE)
t2 <- Sys.time()

time.sa1 <- t2 - t1

#instance
#inisol, with a label
# Tmax equal for all runs
# muValues is a vector of possible values of mu
# runs is the number of runs of each value of mu
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
resultsSANN <- ExperimentSA(sample$distances, NN$sol, label="NN", Tmax=50000, muValues = c(1, 10, 50), runs=10)

set.seed(1111)
solution <- sample(1:50, 50)
resultsSArandom <- ExperimentSA(sample$distances, solution, label="random", Tmax=50000, muValues = c(1, 10, 50), runs=10)

#binding data from NN and random starting solutions
allresults <- rbind(resultsSANN$results, resultsSArandom$results)

# https://dplyr.tidyverse.org/
# https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html
library(dplyr)

#turning results into a tibble (not necessary)
dataSA <- tbl_df(allresults)


dataSA %>% group_by(initialSol) %>% summarise(max=max(distance), mean=mean(distance), min=min(distance))
dataSA %>% group_by(mu) %>% summarise(max=max(distance), mean=mean(distance), min=min(distance))

# http://ggplot2.tidyverse.org/
# https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf
library(ggplot2)

ggplot(dataSA, aes(initialSol, distance)) + geom_boxplot()
ggplot(dataSA, aes(factor(mu), distance)) + geom_boxplot()

ggplot(dataSA, aes(factor(mu), distance)) + geom_boxplot() + facet_grid(. ~ initialSol)

dataSA %>% filter(mu==1) %>% ggplot(aes(distance)) + geom_density()

save.image("results/experimentSATSP.RData")
