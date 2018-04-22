#-------defining objective function------- 

makespan <- function(mat, p){
  m <- dim(mat)[1]
  n <- dim(mat)[2]
  
  if(length(p)!=n) return(success=FALSE)
  
  #reorder columns of m accordin to sequence
  mat <- mat[, p]
  
  #create matrix
  mak <- matrix(numeric(m*n), m, n)
  
  mak[1,1] <- mat[1,1]
  for(i in 2:m) mak[i,1] <- mak[i-1,1] + mat[i,1]
  for(j in 2:n) mak[1,j] <- mak[1,j-1] + mat[1,j]
  for(i in 2:m){
    for(j in 2:n) mak[i,j] <- max(mak[i-1,j], mak[i, j-1]) + mat[i,j] 
  }
  
  return(mak[m,n])
  
}

#----swap of two elements (neighbourhood definiton)----

swap <- function(v, i , j){
  
    aux <- v[i]
    v[i] <- v[j]
    v[j] <- aux
    return(v)
}

#-----hill climbing----

HillClimbing <- function(mat, sini){
  
  n <- length(sini)
  sol <- sini
  obj <- makespan(mat, sol)
  k <- TRUE
  
  while(k){
    #obtaining the solution of minimum value of o.f. from neighbourhood
    soltest <- numeric(n)
    objtest <- Inf
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        test <- makespan(mat, swap(sol, i, j))
        if(test < objtest){
          objtest <- test
          soltest <- swap(sol, i, j)
        }
      }
    }
    
    #comparing the obtained solution with the obtained in previous iteration
    #loop breaks if there is no improvement
    if(objtest < obj){
      obj <- objtest
      sol <- soltest
    }
    else
      k <- FALSE
  }
  
  return(list(sol=sol, obj=obj))
}

#-----simulated annealing-----

SA <- function(mat, sini, Tmax, mu, report.evol=FALSE){
  
  n <- length(sini)
  sol <- sini
  obj <- makespan(mat, sol)
  
  solbest <- sol
  objbest <- obj
  
  evol <- numeric(Tmax)
  
  T <- Tmax
  
  while(T>0){
    move <- sample(1:n, 2)
    soltest <- swap(sol, move[1], move[2])
    objtest <- makespan(mat, soltest)
    
    print(obj)
    print(objtest)
    print(T)
    if(exp(-mu*(objtest-obj)/T) > runif(1)){
      sol <- soltest
      obj <- objtest
    }
    
    if(obj < objbest){
      solbest <- sol
      objbest <- obj
    }
    
    evol[Tmax-T] <- obj
    T <- T-1
  }
  
  if(report.evol)
    return(list(sol=solbest, obj=objbest, evol=evol))
  else
    return(list(sol=solbest, obj=objbest))
}

#----Tabu search----

TS <- function(mat, sini, iter, tabu.size=7, report.evol=FALSE){
  
  n <- length(sini)
  sol <- sini
  obj <- makespan(mat, sol)
  evol <- numeric(iter)
  
  solbest <- sol
  objbest <- obj
  
  k <- 1
  
  tabu.list <- matrix(numeric(tabu.size*2), tabu.size, 2)
  
  
  while(k<=iter){
    
    #we need to find the best solution of the iteration
    #if the aspiration condition is reached, tabu list will not be updated
    move <- numeric(2)
    objiter <- Inf
    soliter <- numeric(n)
    update.tabu <- FALSE
    
    #examining all moves to find best of iteration
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        #compute move solution and objective function
        soltest <- swap(sol, i, j)
        objtest <- makespan(mat, soltest)
        
        #see if move is tabu
        is.tabu <- FALSE
        for(s in 1:tabu.size) if(tabu.list[s,1]==i & tabu.list[s,2]==j) is.tabu <- TRUE
        
        if(is.tabu){
          #move is tabu
          #aspiration condition: replace if better than best
          if(objtest < objbest){
            soliter <- soltest
            objiter <- objtest
            update.tabu <- FALSE
          }
        }
        else{
          #move is not tabu
          #replace iter solution if better than iter
          #tabu list has to be updated, save move
          if(objtest < objiter){
            soliter <- soltest
            objiter <- objtest
            update.tabu <- TRUE
            move <- c(i,j)
          }
        }
      }
    }
    
    #update solution to explore
    sol <- soliter
    obj <- objiter
    
    #update best solution found
    if(obj < objbest){
      objbest <- obj
      solbest <- sol
    }
    
    #update tabu list
    if(update.tabu){
      for(s in tabu.size:2) tabu.list[s, ] <- tabu.list[s-1, ]
      tabu.list[1, ] <- move
    }
    
    evol[k] <- obj
    k <- k+1
    
  }
  
  if(report.evol)
    return(list(sol=solbest, obj=objbest, evol=evol))
  else
    return(list(sol=solbest, obj=objbest))
}

#-----Genetic Algorithm for the TSP------

#function that calculates the value of objective function for each element of the population, and also returns the probability distribution for selection
obj <- function(mat, population){
  npop <- length(population)
  fitness <- sapply(population, function(x) makespan(mat, x))
  
  max.fitness <- max(fitness)
  min.fitness <- min(fitness)
  
  if(max.fitness==min.fitness){
    prob <- rep(1/npop, npop)
  }
  else{
    prob <- (max.fitness - fitness + min.fitness)^2
    sum.prob <- sum(prob)
    prob <- prob/sum.prob  
  }
  
  return(list(fitness=fitness, prob=prob))
}


#function that performs the 1-point crossover operator
crossover1point <- function(sol1, sol2){
  n <- length(sol1)
  son <- numeric(n)
  point <- sample(1:n, 1)
  
  son <- c(sol1[1:point], sol2[!(sol2 %in% sol1[1:point])])
  
  return(son)
  
}

#function that performs a 2-points swap mutation operator (same as SwapNodes)
mutation2nodes <- function(sol){
  
  n <- length(sol)
  points <- sort(sample(1:n, 2))
  i <- points[1]
  j <- points[2]
  
  swap <- sol
  swap[i] <- sol[j]
  swap[j] <- sol[i]
  
  if(i==1 & j!=1){
    swap <- swap[c(j:n, 1:(j-1))]
  }
  
  return(swap)
}

GA <- function(mat, npop, iterations, pmut, report.evol=FALSE, verbose=FALSE){
  
  #problem size  
  n <- dim(mat)[2]
  
  #evolution of best fit
  evol <- numeric(iterations)
  
  #initial population
  population.parent <- replicate(npop, sample(1:n), simplify=FALSE)
  
  #initializing next generation
  population.son <- replicate(npop, numeric(n), simplify=FALSE)
  
  best.obj <- Inf
  best.sol <- numeric(n)
  
  count <- 1
  
  while(count <= iterations){
    
    #assess fitness of population parent
    fitness <- obj(mat, population.parent)
    iter.obj <- min(fitness$fitness)
    pos.min <- which(fitness$fitness == min(fitness$fitness))[1]
    iter.sol <- population.parent[[pos.min]]
    #print(iter.sol)
    
    #elitist strategy: replace worst solution by best solution so far
    pos.max <- which(fitness$fitness == max(fitness$fitness))[1]
    population.son[[pos.max]] <- best.sol
    
    #update best solution
    if(iter.obj < best.obj){
      best.sol <- iter.sol
      best.obj <- iter.obj
    }
    
    #creating son population
    for(i in 1:npop){
      #crossover
      parents <- sample(1:npop, 2, prob = fitness$prob)
      v1 <- population.parent[[parents[1]]]
      v2 <- population.parent[[parents[2]]]
      population.son[[i]] <- crossover1point(v1, v2)
      #mutation
      if(pmut > runif(1)) population.son[[i]] <- mutation2nodes(population.son[[i]])
    }
    
    
    #population son becomes parent
    population.parent <- population.son
    
    evol[count] <- best.obj
    count <- count + 1
  }
  
  if(verbose) print(c(pmut, npop))
  
  if(report.evol)
    return(list(sol=best.sol, obj=best.obj, evol=evol))
  else
    return(list(sol=best.sol, obj=best.obj))
  
}


