#-----matrix distance from a set of points-----

Distancia <- function(df, euc=TRUE){
  
  names(df) <- c("x", "y")
  n <- nrow(df)
  D <- matrix(numeric(n*n), n, n)
  
  if(euc){#euc=TRUE
    for(i in 1:n){
      for(j in 1:n){
        D[i,j] <- sqrt((df$x[i]-df$x[j])^2 + (df$y[i]-df$y[j])^2) 
      }
    }  
    return(D)
  }else{#euc=FALSE
    for(i in 1:n){
      for(j in 1:n){
        D[i,j] <- abs(df$x[i]-df$x[j]) + abs(df$y[i]-df$y[j])
      } 
    }
    return(D)
  }
}

#------Function that obtains distance matrix for Att TSP instances----

DistanceMatrixAtt <- function(data){
  n <- nrow(data)
  d <- matrix(numeric(), n, n)
  
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      xd <- data$xcoord[i] - data$xcoord[j]
      yd <- data$ycoord[i] - data$ycoord[j]
      
      rij <- sqrt((xd*xd + yd*yd) /10)
      
      tij <- round(rij) #caveat: rounding to nearest event number
      
      if(tij < rij){
        d[i, j] <- tij + 1}
      else{ 
        d[i, j] <- tij}
      
      d[j, i] <- d[i, j]
    }
  }
  
  for(i in 1:n)
    d[i, i] <- Inf
  
  return(d)
}

#----TSP instance generation from a set of points----

SampleTSP <- function(size, seed=1111, euc=TRUE){
  set.seed(seed)
  points <- matrix(sample(0:100, size*2, replace = TRUE), size, 2)
  points <- as.data.frame(points)
  d <- Distancia(points, euc)
  return(d)
}

#----Objective function----

TSP <- function(D, sol){
  
  if((dim(D)[1]==dim(D)[2]) & (dim(D)[1] == length(sol))){
    n <- length(sol)
    d <- 0
    for(i in 1:(n-1)) d <- d + D[sol[i],sol[i+1]]
    d <- d + D[sol[n],sol[1]] 
    return(list(d=d,SUCCESS=TRUE))
  }else{
    return(SUCCESS=FALSE)
  }
} 

#-----NearestNeighbour function for the TSP----
#obtaining solution by best neighbour (starting from node start) of tsp for a distance matrix
#returns a list with path obtained and total distance

NearestNeighbour <- function(d, start=1){
  
  n <- dim(d)[1]
  
  #test if d is a square matrix
  #test if start <= n
  
  for(i in 1:n)
    d[i, i] <- Inf
  
  #flag checks if node has been included in solution
  #sol is the solution starting from start
  
  flag <- logical(n)
  sol <- integer(n)
  
  sol[1] <- start
  flag[start] <- TRUE
  
  for(i in 1:(n-1)){
    
    min <- Inf
    k <- sol[i]
    for(j in 1:n){
      if(d[k, j] < min & flag[j]==FALSE){
        min <- d[k, j]
        sol[i+1] <- j
      }
    }
    flag[sol[i+1]] <- TRUE
  }
  
  distance <- TSP(d, sol)
  return(list(sol=sol, obj=distance))
  
}

#----swap of two elements (neighbourhood definiton)----

swap <- function(v, i , j){
  
  n <- length(v)
  if(i > j){
    aux <- i
    i <- j
    j <- aux
  }
  
  if((n >= j)&(i!=j)) {
    aux <- v[i]
    v[i] <- v[j]
    v[j] <- aux
    if(i==1) v <- v[c(j:n, 1:(j-1))]
    return(list(v=v, success=TRUE))
  }    
  else{
    if(i==j & j <=n) return(list(v=v, success=TRUE))
    if(n < j) return(success=FALSE)
  }  
}

#-----hill climbing----

HillClimbing <- function(sini, D){
  
  n <- length(sini)
  sol <- sini
  dist <- TSP(D, sol)$d
  k <- TRUE
  
  while(k){
    #obtaining the solution of minimum value of o.f. from neighbourhood
    soltest <- numeric(n)
    disttest <- Inf
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        test <- TSP(D, swap(sol, i, j)$v)$d
        if(test < disttest){
          disttest <- test
          soltest <- swap(sol, i, j)$v
        }
      }
    }
    
    #comparing the obtained solution with the obtained in previous iteration
    #loop breaks if there is no improvement
    if(disttest < dist){
      dist <- disttest
      sol <- soltest
    }
    else
      k <- FALSE
  }
  
  return(list(sol=sol, dist=dist))
}

#-----simulated annealing-----

SA <- function(sini, D, Tmax, mu, report.evol=FALSE){
    
    n <- length(sini)
    sol <- sini
    dist <- TSP(D, sol)$d
    
    solbest <- sol
    distbest <- dist
    
    evol <- numeric(Tmax)
    
    T <- Tmax
    
    while(T>=0){
        move <- sample(1:n, 2)
        soltest <- swap(sol, move[1], move[2])$v
        disttest <- TSP(D, soltest)$d
        
        #print(soltest)
        if(exp(-mu*(disttest-dist)/T) > runif(1)){
            sol <- soltest
            dist <- disttest
        }
        
        if(dist < distbest){
            solbest <- sol
            distbest <- dist
        }
        
        evol[Tmax-T] <- dist
        T <- T-1
    }
    
    if(report.evol)
        return(list(sol=solbest, dist=distbest, evol=evol))
    else
        return(list(sol=solbest, dist=distbest))
}

#----Tabu search----

TS <- function(sini, D, iter, tabu.size=7, report.evol=FALSE){
  
  n <- length(sini)
  sol <- sini
  dist <- TSP(D, sol)$d
  
  evol <- numeric(iter)
  
  solbest <- sol
  distbest <- dist
  k <- 1
  
  tabu.list <- matrix(numeric(tabu.size*2), tabu.size, 2)
  
  
  while(k<=iter){
    
    #we need to find the best solution of the iteration
    #if the aspiration condition is reached, tabu list will not be updated
    move <- numeric(2)
    distiter <- Inf
    soliter <- numeric(n)
    update.tabu <- FALSE
    
    #examining all moves to find best of iteration
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        #compute move solution and objective function
        soltest <- swap(sol, i, j)$v
        disttest <- TSP(D, soltest)$d
        
        #see if move is tabu
        is.tabu <- FALSE
        for(s in 1:tabu.size) if(tabu.list[s,1]==i & tabu.list[s,2]==j) is.tabu <- TRUE
        
        if(is.tabu){
          #move is tabu
          #aspiration condition: replace if better than best
          if(disttest < distbest){
            soliter <- soltest
            distiter <- disttest
            update.tabu <- FALSE
          }
        }
        else{
          #move is not tabu
          #replace iter solution if better than iter
          #tabu list has to be updated, save move
          if(disttest < distiter){
            soliter <- soltest
            distiter <- disttest
            update.tabu <- TRUE
            move <- c(i,j)
          }
        }
      }
    }
    
    #update solution to explore
    sol <- soliter
    dist <- distiter
    
    #update best solution found
    if(dist < distbest){
      distbest <- dist
      solbest <- sol
    }
    
    #update tabu list
    if(update.tabu){
      for(s in tabu.size:2) tabu.list[s, ] <- tabu.list[s-1, ]
      tabu.list[1, ] <- move
    }
    
    evol[k] <- dist
    k <- k+1
    
  }
  
  if(report.evol)
      return(list(sol=solbest, obj=distbest, evol=evol))
  else
      return(list(sol=solbest, obj=distbest))
}

#----A GRASP constructive heuristic for the TSP based on nearest neighbour heuristic----

GRASP <- function(d, s, start=1){
  
  n <- dim(d)[1]
  
  #test if d is a square matrix
  #test if start <= n
  
  for(i in 1:n)
    d[i, i] <- Inf
  
  #flag checks if node has been included in solution
  #sol is the solution starting from start
  
  flag <- logical(n)
  sol <- integer(n)
  
  sol[1] <- start
  flag[start] <- TRUE
  
  for(i in 1:(n-1)){
    
    row <- d[sol[i], ]
    pos <- which(flag==FALSE)
    list <- pos[order(row[pos], decreasing = FALSE)]
    sizelist <- min(s, n-sum(flag))
    
    rank <- sample(1:sizelist, 1)
    j <- list[rank]
    sol[i+1] <- j
    flag[j] <- TRUE
    
  }
  
  distance <- TSP(d, sol)
  return(list(sol=sol, obj=distance))
  
}

#---- combining GRASP constructive heuristic with TS ----

#-----Genetic Algorithm for the TSP------

#function that calculates the value of objective function for each element of the population, and also returns the probability distribution for selection
obj <- function(d, population){
  npop <- length(population)
  fitness <- sapply(population, function(x) TSP(d, x)$d)
  
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
crossoverTSP1point <- function(sol1, sol2){
  n <- length(sol1)
  son <- numeric(n)
  point <- sample(1:n, 1)
  
  son <- c(sol1[1:point], sol2[!(sol2 %in% sol1[1:point])])
  
  return(son)
  
}

#function that performs a 2-points swap mutation operator (same as SwapNodes)
mutationTSP2nodes <- function(sol){
  
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

GA <- function(d, npop, iterations, pmut, seed=FALSE, report.evol=FALSE, verbose=FALSE){
  
  #problem size  
  n <- dim(d)[1]
  
  #evolution of best fit
  evol <- numeric(iterations)
  
  #initial population
  population.parent <- replicate(npop, c(1, sample(2:n,size=(n-1))), simplify=FALSE)
  
  #initializing next generation
  population.son <- replicate(npop, numeric(n), simplify=FALSE)
  
  #include solution with nearest neighbour if seed==TRUE
  if(seed == TRUE){
    population.parent[[1]] <- NearestNeighbour(d, 1)$sol
  }
  
  best.obj <- Inf
  best.sol <- numeric(n)
  
  count <- 1
  
  while(count <= iterations){
    
    #assess fitness of population parent
    fitness <- obj(d, population.parent)
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
      population.son[[i]] <- crossoverTSP1point(v1, v2)
      #mutation
      if(pmut > runif(1)) population.son[[i]] <- mutationTSP2nodes(population.son[[i]])
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
