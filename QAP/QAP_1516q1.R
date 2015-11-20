#-------running QAP functions------ 
setwd("~/Dropbox/quanti_docs/QAP")
source("QAPfunctions.R")

#------reading instances------

url.base <- "http://www.opt.math.tu-graz.ac.at/qaplib/data.d/"

link <- function(x) paste0(url.base, x, ".dat")

sample.med <- c("rou12", "rou15", "rou20", "scr12", "scr15", "scr20")
link.med <- lapply(sample.med, link)
data.med <- lapply(link.med, read.QAPLIB)
names(data.med) <- sample.med

sample.tfda <- c("had12", "had16", "had18", "had20")
link.tfda <- lapply(sample.tfda, link)
data.tfda <- lapply(link.tfda, function(x) read.QAPLIB(x, swap="TRUE"))
names(data.tfda) <- sample.tfda

data <- c(data.med, data.tfda)
sample <- c(sample.med, sample.tfda)

rm(data.med, data.tfda, sample.med, sample.tfda, link.med, link.tfda)

#----reading optimal solutions------

url.sol <- "http://www.opt.math.tu-graz.ac.at/qaplib/soln.d/"

link.sol <- function(x) paste0(url.sol, x, ".sln")

sol.link <- lapply(sample, link.sol)
sol.opt <- lapply(sol.link, readsol.QAPLIB)
names(sol.opt) <- sample

#-----solving by constructive heuristic-----

solve.by.heuristic <- function(x){
  flow <- x$flow
  distance <- x$distance
  sol.heuristic <- Heuristic.QAP(flow, distance)
  return(list(obj = sol.heuristic$obj, sol=sol.heuristic$sol, time=0))
}

sol.heuristic <- lapply(data, solve.by.heuristic)
names(sol.heuristic) <- sample

#-----solve by simulated annealing-----

solve.by.sa <- function(x){
  flow <- x$flow
  distance <- x$distance
  
  starting <- Heuristic.QAP(flow, distance)
  
  set.seed(1313)
  time1 <- Sys.time()
  sol.heuristic <- SimulatedAnnealing.QAP(flow, distance, 30000, starting$sol)
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(obj = sol.heuristic$obj, sol=sol.heuristic$sol, time=time))
  
  return(sol.heuristic)
}

sol.sa <- lapply(data, solve.by.sa)
names(sol.sa) <- sample

#-----solve by tabu search----

solve.by.ts <- function(x){
  flow <- x$flow
  distance <- x$distance
  
  starting <- Heuristic.QAP(flow, distance)
  
  set.seed(1313)
  time1 <- Sys.time()
  sol.heuristic <- TabuSearch.QAP(flow, distance, iter=30000, sol=starting$sol)
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(obj = sol.heuristic$obj, sol=sol.heuristic$sol, time=time))
  
  return(sol.heuristic)
}

sol.ts <- lapply(data, solve.by.ts)
names(sol.ts) <- sample

#------solve by genetic algorithm-----

solve.by.ga <- function(x){
  flow <- x$flow
  distance <- x$distance
  
  set.seed(1313)
  time1 <- Sys.time()
  sol.heuristic <- GeneticAlgorithm.QAP(flow, distance)
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(obj = sol.heuristic$obj, sol=sol.heuristic$sol, time=time))
  
  return(sol.heuristic)
}

sol.ga <- lapply(data, solve.by.ga) 
names(sol.ga) <- sample

#-----turning solution (permutation) into string-----

SolString <- function(sol){
  text <- character(0)
  for (i in 1:length(sol)) text <- paste(text, sol[i])
  text <- substring(text, 2)
  return(text)
}

for(i in 1:10){
  sol.opt[[i]]$sols <- SolString(sol.opt[[i]]$sol)
  sol.heuristic[[i]]$sols <- SolString(sol.heuristic[[i]]$sol)
  sol.sa[[i]]$sols <- SolString(sol.sa[[i]]$sol)  
  sol.ts[[i]]$sols <- SolString(sol.ts[[i]]$sol)
  sol.ga[[i]]$sols <- SolString(sol.ga[[i]]$sol) 
  
  sol.opt[[i]]$eff <- 0
  sol.heuristic[[i]]$eff <- (sol.heuristic[[i]]$obj - sol.opt[[i]]$obj)/sol.opt[[i]]$obj
  sol.sa[[i]]$eff <- (sol.sa[[i]]$obj - sol.opt[[i]]$obj)/sol.opt[[i]]$obj
  sol.ts[[i]]$eff <- (sol.ts[[i]]$obj - sol.opt[[i]]$obj)/sol.opt[[i]]$obj
  sol.ga[[i]]$eff <- (sol.ga[[i]]$obj - sol.opt[[i]]$obj)/sol.opt[[i]]$obj
  
}

results <- list(0)

n <- length(sample)

for(i in 1:n){
  results[[i]] <- data.frame(matrix(nrow=5, ncol=0))
  
  results[[i]]$obj <- c(sol.heuristic[[i]]$obj, sol.sa[[i]]$obj, sol.ts[[i]]$obj, sol.ga[[i]]$obj, sol.opt[[i]]$obj)
  
  results[[i]]$sol <- c(SolString(sol.heuristic[[i]]$sol), SolString(sol.sa[[i]]$sol), SolString(sol.ts[[i]]$sol), SolString(sol.ga[[i]]$sol), SolString(sol.opt[[i]]$sol))
  
  results[[i]]$time <- c(sol.heuristic[[i]]$time, sol.sa[[i]]$time, sol.ts[[i]]$time, sol.ga[[i]]$time, sol.opt[[i]]$time)
  
  results[[i]]$time <- round(results[[i]]$time, 2)
  
  for(j in 1:5){
    results[[i]]$eff[j] <- (results[[i]]$obj[j] - results[[i]]$obj[5])/results[[i]]$obj[5]
    results[[i]]$eff[j] <- round(results[[i]]$eff[j], 2)
  }
  rownames(results[[i]]) <- c(paste("Heuristic", sample[[i]]), paste("Sim. Annealing", sample[[i]]), paste("Tabu", sample[[i]]), paste("Genetic", sample[[i]]), paste("Optimal", sample[[i]]))
}


library(stargazer)

sink("resultsQAP.tex", append=FALSE, split=FALSE)
for(i in 1:10) stargazer(results[[i]], summary=FALSE, dep.var.caption=paste("Results for", sample[i]))
sink()
save.image("QAPresults.RData")
load("QAPresults.Rdata")

