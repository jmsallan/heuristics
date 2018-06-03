#---- reading SCP instances ----

# set-covering instances from OR-library
# http://people.brunel.ac.uk/~mastjjb/jeb/orlib/scpinfo.html

# The format of all of these 80 data files is:
# number of rows (m), number of columns (n)
# the cost of each column c(j),j=1,...,n
# for each row i (i=1,...,m): the number of columns which cover
# row i followed by a list of the columns which cover row i

readCSP <- function(url){
    a <- readLines(url)
    b <- lapply(a, function(x) unlist(strsplit(x, " ")))
    c <- lapply(b, function(x) as.numeric(x[-1]))
    d <- sapply(c, function(x) length(x))
    m <- c[[1]][1]
    n <- c[[1]][2]

    counter <- 2
    
    elements <- 0
    nextcounter <- counter
    while(elements < n){
      elements <- elements + d[nextcounter]
      nextcounter <- nextcounter + 1
    }
    
    costs <- numeric(0)
    for(i in counter:(nextcounter-1)) costs <- c(costs, c[[i]])
    
    counter <- nextcounter
    
    M <- matrix(rep(0, m*n), m, n)
    
    for(i in 1:m){
    
      elements <- 0
      nextcounter <- nextcounter+1
      while(elements < c[[counter]]){    
        elements <- elements + d[nextcounter]
        nextcounter <- nextcounter + 1
      }
      
      for(j in (counter+1):(nextcounter-1)){
        M[i , c[[j]]] <- 1
      }
      counter <- nextcounter
    }

return(list(M=M, costs=costs))
}

#---- solving LP with linear programming ----

solveLP_SCP <- function(instance){
  
  obj <- instance$costs
  mat <- instance$M
  
  m <- dim(mat)[1]
  n <- length(obj)
  
  dir <- rep(">=", m)
  rhs <- rep(1, m)
  types <- rep("B", n)
  
  solution <- Rglpk_solve_LP(obj=obj, mat=mat, dir=dir, rhs=rhs, types=types, max=FALSE)
  
  optcost <- solution$optimum
  
  optsolution <- which(solution$solution==1)
  
  return(list(optcost=optcost, optsolution=optsolution))
  
}

