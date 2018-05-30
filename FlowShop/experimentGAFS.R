#set wd!

source("code/functionsFS2.R")

#---- experimental design of a GA for the Flowshop ----

#values of npop: 10, 50, 100
#values of pmut: 0.4, 0.8, 1
#two crossover operators: OX, 2point crossover
#every combination is performed 5 runs

experimentFS <- function(Instance, npopValues, pmutValues, crOXValues, runs){
  
  totalruns <- length(npopValues)*length(pmutValues)*length(crOXValues)*runs
  
  results <- as.data.frame(matrix(numeric(totalruns*4), totalruns, 4))
  names(results) <- c("npop", "pmut", "OX cross", "makespan")
  sols <- list()
  count <- 1
  
  for( r in 1:runs){
    for(npop in npopValues){
      for(pmut in pmutValues){
        for(crOX in crOXValues){
            test <- GAFS(problem = Instance, npop = npop, iter=100, pmut = pmut, crOX = crOX, alpha=0.1)
            results[count, ] <- c(npop, pmut, crOX, test$fit)
            sols[[count]] <- test$sol
            count <- count+1
        }
      }
    }
  }
  
  return(list(results=results, sols=sols))
  
}

load("instances/TaillardFS.RData")

results <- experimentFS(tai100.5[[1]]$tij, npopValues = c(10, 20, 50), pmutValues = c(0.4, 0.8, 1), crOXValues = c(TRUE, FALSE), runs = 5)

#saving the results on RDS

saveRDS(results, "results/GAtai100_5_1.rds")

#---- extracting results ----
