#setwd!

source("code/functionsTSP.R")

#----- a 2-opt move for the TSP ----

move2opt <- function(sol, i, k){
  n <- length(sol)
  move <- c(sol[1:(i-1)], rev(sol[i:k]),  sol[(k+1):n])
  
  if(i==n){
    
    move <- move[-1]
    move <- c(move[k:n], move[1:(k-1)])
  }
  
  return(move)
}

vec <- 1:7

move2opt(vec, 3, 5)
move2opt(vec, 2, 6)
move2opt(vec, 1, 3)

# all 2opt moves:

n <- 7
k <- 0
for(i in 2:(n-1)){
  for(j in (i+2):(n+1)){
    
    print(c(i,j))
    # k <- k+1
    # print(k)
}
}
#---- a sample instance of size 10 ----

Instance10 <- SampleTSP(10)

sol110 <- TSP(Instance10, 1:10)
solnn <- NearestNeighbour(Instance10)

HillClimbing2opt <- function(sini, D){
  
  
  
  
}
