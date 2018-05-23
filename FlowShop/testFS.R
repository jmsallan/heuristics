# set wd!
  
source("code/functionsFS2.R")

#----testing makespan function----

mInstance <- matrix(c(8,3,6,9,
                      4,5,3,2,
                      1,2,7,4), 3, 4, byrow = TRUE)

solution <- c(4, 3, 1, 2)

makespan(mInstance, solution)

set.seed(1313)
Instance <- matrix(sample(10:50, 100, replace = TRUE), 5, 20)

makespans <- sapply(1:1000, function(x) makespan(Instance, sample(1:20, 20)))

makespan(Instance, 1:20)

set.seed(1111)
perms <- lapply(1:1000, function(x) sample(1:20, 20))

sols <- lapply(perms, function(x) makespan(Instance, x))
solsb <- sapply(perms, function(x) makespan(Instance, x))

hist(solsb, xlab="", ylab="", main="makespan")

#----crossover and mutation operators-----


sol1 <- sample(1:10, 10)
sol2 <- sample(1:10, 10)
crossover2pointI(sol1, sol2)

#shift change

sol <- 1:10
shift(sol)

#----testing the GA----

GA01 <- GAFS(Instance, npop=10, iter=1000, pmut=1, verbose = TRUE)

load("instances/TaillardFS.RData")

GA02 <- GAFS(tai20.5[[1]]$tij, npop=10, iter=1000, pmut=1, verbose = TRUE)


#----NEH heuristic----

NEH <- function(Instance){
  
  m <- dim(Instance)[1]
  n <- dim(Instance)[2]
  
  tasks <- order(apply(Instance, 2, sum))
  
  sol <- tasks[1]
  
  for(i in 2:(n-1)){
    
    best <- Inf
    pos <- -1
    
    for(j in 0:(i-1)){
      
      test <- c(sol[0:j], tasks[i], sol[(j+1):i])
      if (makespan(Instance, test) < best){
        pos <- j
        best <- makespan(Instance, test)
      }
      
      if(makespan(Instance, c(sol, task[i])) < best){
        pos <- i
        best <- makespan(Instance, c(sol, task[i])) 
      }
      
      if(pos < i)
        sol <- c(sol[0:pos], task[i], sol[(pos+1):i])
      else
        sol <- c(sol, task[i])
    }
  }
  return(sol)
}



#----shift move----

test <- 1:5
for(i in 1:5)
  for(j in 1:5)
    if (i!=j) {print(c(i,j))
      print(shiftmove(test, i, j))}
