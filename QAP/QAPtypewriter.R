#-------running QAP functions------ 
setwd("~/Dropbox/quanti_docs/QAP")
source("QAPfunctions.R")

#------reading instances------

url.base <- "http://www.opt.math.tu-graz.ac.at/qaplib/data.d/"

link <- function(x) paste0(url.base, x, ".dat")

sample <- c("bur26a", "bur26b", "bur26c", "bur26d", "bur26e", "bur26f", "bur26g", "bur26h")
link.sample <- lapply(sample, link)
data <- lapply(link.sample, function(x) read.QAPLIB(x, swap="TRUE"))

#----reading optimal solutions------

url.sol <- "http://www.opt.math.tu-graz.ac.at/qaplib/soln.d/"

link.sol <- function(x) paste0(url.sol, x, ".sln")

sol.link <- lapply(sample, link.sol)
sol.opt <- lapply(sol.link, readsol.QAPLIB)

#-----solving by constructive heuristic-----

solve.by.heuristic <- function(x){
  flow <- x$flow
  distance <- x$distance
  sol.heuristic <- Heuristic.QAP(flow, distance)
  return(list(obj = sol.heuristic$obj, sol=sol.heuristic$sol, time=0))
}

sol.heuristic <- lapply(data, solve.by.heuristic)

#-----solve by simulated annealing-----

solve.by.sa <- function(x){
  flow <- x$flow
  distance <- x$distance
  
  starting <- Heuristic.QAP(flow, distance)
  
  set.seed(1313)
  time1 <- Sys.time()
  sol.annealing <- SA(starting$sol, flow, distance, 30000, 3000)
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(obj = sol.annealing$obj, sol=sol.annealing$sol, time=time))
}

sol.sa <- lapply(data, solve.by.sa)

#-----solve by tabu search----

solve.by.ts <- function(x){
  flow <- x$flow
  distance <- x$distance
  
  starting <- Heuristic.QAP(flow, distance)
  
  set.seed(1313)
  time1 <- Sys.time()
  sol.tabu <- TS(starting$sol, flow, distance, iter=300)
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(obj = sol.tabu$obj, sol=sol.tabu$sol, time=time))
}

sol.ts <- lapply(data, solve.by.ts)

#------solve by genetic algorithm-----

solve.by.ga <- function(x){
  flow <- x$flow
  distance <- x$distance
  
  set.seed(1313)
  time1 <- Sys.time()
  sol.genetic <- GA(flow, distance, 300, 200, 0.2)
  time2 <- Sys.time()
  time <- time2 - time1
  return(list(obj = sol.genetic$obj, sol=sol.genetic$sol, time=time))
}

sol.ga <- lapply(data, solve.by.ga) 


#-----turning solution (permutation) into string-----

SolString <- function(sol){
  text <- character(0)
  for (i in 1:length(sol)) text <- paste(text, sol[i])
  text <- substring(text, 2)
  return(text)
}

#----presenting the results---- 

n <- length(sample)

results <- list(0)

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
for(i in 1:length(sample)) stargazer(results[[i]], type="text", summary=FALSE)

sink("QAPtypewriter.txt")
for(i in 1:length(sample)) stargazer(results[[i]], summary=FALSE)
sink()

save.image("QAPresultsTypewriters.RData")
load("QAPresultsTypewriters.Rdata")

