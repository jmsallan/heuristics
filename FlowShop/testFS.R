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

Instance.NEH <- matrix(c(
  5, 9, 9, 4,
  9, 3, 4, 8,
  8, 10, 5, 8,
  10, 1, 8, 7,
  1, 8, 6, 2), 5, 4, byrow=TRUE)

NEH(Instance.NEH)

NEH(tai20.5[[1]]$tij)

#----shift move----

test <- 1:5
for(i in 1:5)
  for(j in 1:5)
    if (i!=j) {print(c(i,j))
      print(shiftmove(test, i, j))}
