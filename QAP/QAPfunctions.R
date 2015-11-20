#-------defining objective function------- 
#the components of the solution design the facility that goes in each location

perMatrix <- function(vec){
  n <- length(vec)
  m <- matrix(0, n , n)
  for(i in 1:n){
    m[i, vec[i]] <- 1
  }
  return(m)
}

objective.QAP <- function(flow, distance, permutation){
  pmatrix <- perMatrix(permutation)
  flow_perm <- pmatrix %*% flow %*% t(pmatrix)
  return(sum(flow_perm * distance))
}

#--------defining constructive heuristic-------
#starts assigning the highest flow value to shortest distance
#then for remaining elements:
#maximal flow from selected facility (Fa) to not selected (Fb) is detected
#facility Fb is placed as close as possible to Fa 

Heuristic.QAP <- function(flow, distance){
  
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
  
  obj <- objective.QAP(flow, distance, sol)
  
  return(list(obj = obj, sol = sol, time = 0))
}

#-----reading data from QAPLIB file-----
#by default first matrix is flow, second matrix distance. In some QAPLIB instances is the opposite, then swap parameter must be turned to TRUE 

read.QAPLIB <- function(link, swap=FALSE){
  
  file <- url(link)
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

readsol.QAPLIB <- function(link){
  file <- url(link)
  data <- scan(file)
  close(file)
  
  n <- data[1]
  obj <- data[2]
  sol <- data[3:(n+2)]
  
  return(list(obj=obj, sol=sol, time=0))
  
}




#------Swap nodes operator (sane as in QAP)

SwapNeighbourNodes <- function(sol, k){
  n <- length(sol)
  swap <- sol
  if(k == n){
    swap[1] <- sol[n]
    swap[n] <- sol[1]
  }
  else{
    swap[k] <- sol[k+1]
    swap[k+1] <- sol[k]
  }
  return(swap)
}

#--------Simmulated annealing heuristic-------

SimulatedAnnealing.QAP <- function(flow, distance, tmax, sol){
  
  #finding problem size
  n <- dim(flow)[1]
  
  sol.best <- sol
  obj <- objective.QAP(flow, distance, sol)
  obj.best <- obj
  
  #tmax number of iterations
  count <- tmax
  
  repeat{
    #select a random move
    pos <- sample(1:n, 1)
    
    #finding new solution and total distance
    sol.test <- SwapNeighbourNodes(sol, pos)
    obj.test <- objective.QAP(flow, distance, sol)
    
    if(exp(-(obj.test - obj)/count) > runif(1)){
      sol <- sol.test
      obj <- obj.test
    }
    
    if(obj < obj.best){
      obj.best <- obj
      sol.best <- sol
    }
    
    count <- count - 1
    if(count <= 0) break
  }
  
  return(list(obj = obj.best, sol = sol.best))
  
}

#----tabu search heuristic---- 

TabuSearch.QAP <- function(flow, distance, iter=100, tabu.size = 7, sol){
  
  #finding problem size
  n <- dim(flow)[1]
  
  sol.best <- sol
  obj <- objective.QAP(flow, distance, sol)
  obj.best <- obj
  
  #tmax number of iterations
  count <- iter
  
  tabu.list <- numeric(tabu.size)
  
  repeat{
    #testing all random nodes
    
    move <- numeric()
    obj.iter <- Inf
    sol.iter <- numeric(n)
    
    for(i in 1:n){
      if(sum(i == tabu.list) == 0){
        sol.test <- SwapNeighbourNodes(sol, i)
        obj.test <- objective.QAP(flow, distance, sol.test)
      }
      if(obj.test < obj.iter){
        sol.iter <- sol.test
        obj.iter <- obj.test
        move <- i
      }
    }    
    
    #update solution and tabu list
    
    sol <- sol.iter
    obj <- obj.iter
    
    for(i in tabu.size:2)
      tabu.list[i] <- tabu.list[i-1]
    
    tabu.list[1] <- move
    
    if(obj < obj.best){
      obj.best <- obj
      sol.best <- sol
    }        
    
    count <- count - 1
    if(count <= 0) break
  }
  
  return(list(obj = obj.best, sol = sol.best))
  
}

#------Genetic algorithm heuristic------

GeneticAlgorithm.QAP <- function(flow, distance, npop=100, iterations=1000, pmut=0.8){
  
  n <- dim(flow)[1]
  
  #initial population
  
  population.father <- replicate(npop, sample(1:n,size=n), simplify=FALSE)
  
  population.mother <- replicate(npop, numeric(n), simplify=FALSE)
  
  sol.best <- numeric(n)
  obj.best <- Inf
  
  fitness <- c(numeric(), npop)
  fitness.acum <- c(numeric(), npop)
  
  iter <- 0
  
  repeat{
    
    #computing relative fitness for selection
    
    for(i in 1:npop)
      fitness[i] <- objective.QAP(flow, distance, population.father[[i]])
    
    if(min(fitness) < obj.best){
      obj.best <- min(fitness)
      sol.best <- population.father[[which.min(fitness)]]
      
    }
    
    max.fitness <- max(fitness)
    
    fitness <- max.fitness - fitness
    
    fitness.acum[i] <- fitness[i]
    
    for(i in 2:npop)
      fitness.acum[i] <- fitness.acum[i-1] + fitness[i]
    
    fitness.acum <- fitness.acum/max(fitness.acum)
    
    for(k in 1:npop){
      #selecting father and mother
      r <- runif(1)
      father <- 0
      
      repeat{
        father <- father + 1
        if(fitness.acum[father] > r | father >= n) break
      }
      
      repeat{
        s <- runif(1)
        mother <-0
        repeat{
          mother <- mother + 1
          if(fitness.acum[mother] > s | mother >= n) break
        }
        if(father!=mother) break
      }
      
      #crossover operator for father and son
      
      crossover.point <- sample(1:n, 1)
      
      c1 <- population.father[[father]][1:crossover.point]
      
      flag <- logical(n)
      
      for(i in 1:n){
        for(j in 1:crossover.point)
          if(population.father[[mother]][i] == c1[j])
            flag[i] = TRUE
      }
      
      c2 <- population.father[[mother]][which(!flag)]
      
      population.mother[[k]] <- c(c1, c2)
      
      #mutation operator
      
      if(runif(1) < pmut){
        swap <- sample(1:n, 2)
        aux <- population.mother[[k]][swap[1]]
        population.mother[[k]][swap[1]] <- population.mother[[k]][swap[2]]
        population.mother[[k]][swap[2]] <- aux
        
      }
    }
    
    population.father <- population.mother
    
    iter <- iter + 1
    if(iter >= iterations) break
  }
  
  for(i in 1:npop)
    fitness[i] <- objective.QAP(flow, distance, population.father[[i]])
  
  if(min(fitness) < obj.best){
    obj.best <- min(fitness)
    sol.best <- population.father[[which.min(fitness)]]
  }
  
  return(list(sol = sol.best, obj = obj.best))
}
