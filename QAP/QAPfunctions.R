#-------defining objective function------- 
#the components of the solution design the location of each facility

perMatrix <- function(vec){
  n <- length(vec)
  m <- matrix(0, n , n)
  for(i in 1:n){
    m[i, vec[i]] <- 1
  }
  return(m)
}

QAP <- function(flow, distance, permutation){
  pmatrix <- perMatrix(permutation)
  flow_perm <- pmatrix %*% flow %*% t(pmatrix)
  return(sum(flow_perm * distance))
}


#-----reading data from QAPLIB file-----
#by default first matrix is flow, second matrix distance. In some QAPLIB instances is the opposite, then swap parameter must be turned to TRUE 

read.QAPLIB <- function(link, swap=FALSE){
  
  file <- url(link, method="libcurl")
  data <- scan(file)
  close(file)
  
  n <- data[1]
  
  flow <- matrix(data[2:(n*n + 1)], n, n, byrow=TRUE)
  
  distance <- matrix(data[(n*n + 2):(2*n*n + 1)], n , n)
  
  if(swap==TRUE){
    aux <- flow
    flow <- distance
    distance <- aux
  }
  
  return(list(flow=flow, distance=distance))
  
}

#----reading solution from QAPLIB-----

readsol.QAPLIB <- function(link, perm=TRUE){
  file <- url(link, method = "libcurl")
  data <- scan(file)
  close(file)
  
  n <- data[1]
  obj <- data[2]
  sol <- data[3:(n+2)]
  
  if(perm){
      aux <- numeric(length(sol))
      for(i in 1:length(sol)) aux[sol[i]] <- i
      sol <- aux
  }
  
  return(list(obj=obj, sol=sol, time=0))
  
}

#--------defining constructive heuristic-------
#starts assigning the highest flow value to shortest distance
#then for remaining elements:
#maximal flow from selected facility (Fa) to not selected (Fb) is detected
#facility Fb is placed as close as possible to Fa 

Heuristic.QAP <- function(flow, distance0){
  
  distance <- distance0
  
  #setting flow and distance matrices as data frames
  
  n <- dim(flow)[1]
  for(i in 1:n) distance[i, i] <- Inf
  
  row <- rep(1:n, n)
  column <- as.vector(sapply(1:n, function(x) rep(x, n)))
  
  distance.df <- data.frame(as.vector(distance), row, column)
  colnames(distance.df) <- c("distance", "row", "column")
  flow.df <- data.frame(as.vector(flow), row, column)
  colnames(flow.df) <- c("flow", "row", "column")
  
  #declaring and initializing values
  
  sol <- rep(0, n)
  
  facility <- rep(FALSE, n)
  location <- rep(FALSE, n)
  
  #assign two initial locations
  
  flow.pos <- which(flow.df$flow == max(flow.df$flow))[1]
  distance.pos <- which(distance.df$distance == min(distance.df$distance))[1]
  
  sol[flow.df[flow.pos, 2]] <- distance.df[distance.pos, 2]
  sol[flow.df[flow.pos, 3]] <- distance.df[distance.pos, 3]
  
  facility[flow.df[flow.pos, 2]] <- TRUE
  facility[flow.df[flow.pos, 3]] <- TRUE
  
  location[distance.df[distance.pos, 2]] <- TRUE
  location[distance.df[distance.pos, 3]] <- TRUE
  
  #assigning rest of locations
  
  for(i in 1:(n-2)){
    
    #elements of flow: row yet chosen, column not chosen
    flow.exam <- flow.df[which((facility[flow.df$row] == TRUE) & (facility[flow.df$column] == FALSE)), ]
    flow.pos <- which(flow.exam$flow == max(flow.exam$flow))[1]
    
    #elements of distance to explore: origin in location of facility departing flow, destiny not assigned
    
    origin <- sol[flow.exam$row[flow.pos]]
    
    distance.exam <- distance.df[which((distance.df$row == origin) & location[distance.df$column] == FALSE), ]
    distance.pos <- which(distance.exam$distance == min(distance.exam$distance))[1]
    
    sol[flow.exam[flow.pos, 3]] <- distance.exam[distance.pos, 3]
    facility[flow.exam[flow.pos, 3]] <- TRUE
    location[distance.exam[distance.pos, 3]] <- TRUE
    
  }
  
  for(i in 1:n) distance[i, i] <- 0
  
  obj <- QAP(flow, distance0, sol)
  
  return(list(obj = obj, sol = sol))
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
    return(list(v=v, success=TRUE))
  }    
  else{
    if(i==j & j <=n) return(list(v=v, success=TRUE))
    if(n < j) return(success=FALSE)
  }  
}

#-----hill climbing----

HillClimbing <- function(sini, flow, distance){
  
  n <- length(sini)
  sol <- sini
  obj <- QAP(flow, distance, sol)
  k <- TRUE
  
  while(k){
    #obtaining the solution of minimum value of o.f. from neighbourhood
    soltest <- numeric(n)
    objtest <- Inf
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        test <- QAP(flow, distance, swap(sol, i, j)$v)
        if(test < objtest){
          objtest <- test
          soltest <- swap(sol, i, j)$v
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

SA <- function(sini, flow, distance, Tmax, mu, report.evol=FALSE){
  
  n <- length(sini)
  sol <- sini
  obj <- QAP(flow, distance, sol)
  
  solbest <- sol
  objbest <- obj
  
  evol <- numeric(Tmax)
  
  T <- Tmax
  
  while(T>=0){
    move <- sample(1:n, 2)
    soltest <- swap(sol, move[1], move[2])$v
    objtest <- QAP(flow, distance, soltest)
    
    #print(soltest)
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

TS <- function(sini, flow, distance, iter, tabu.size=7, report.evol=FALSE){
  
  n <- length(sini)
  sol <- sini
  obj <- QAP(flow, distance, sol)
  
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
        soltest <- swap(sol, i, j)$v
        objtest <- QAP(flow, distance, soltest)
        
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
obj <- function(flow, distance, population){
  npop <- length(population)
  fitness <- sapply(population, function(x) QAP(flow, distance, x))
  
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

GA <- function(flow, distance, npop, iterations, pmut, report.evol=FALSE, verbose=FALSE){
  
  #problem size  
  n <- dim(distance)[1]
  
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
    fitness <- obj(flow, distance, population.parent)
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


