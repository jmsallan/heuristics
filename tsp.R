
#------Function that swaps neighbor nodes from a permutative solution 

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

#-----defining a random matrix distance of size n-----

RandomDMatrix <- function(n){
  m <- matrix( sample(1:10, n*n, replace=TRUE), n, n)
  for(i in 1:n)
    m[i, i] <- Inf
  
  return(m)
}

#------DistanceTSP function for the tsp------
#given a d matrix and a solution, returns the value of total distance

DistanceTSP <- function(d, sol){
  n <- dim(d)[1]
  
  distance <- 0
  
  for(i in 1:(n-1))
    distance <- distance + d[sol[i], sol[i+1]]
  
  distance <- distance + d[sol[n], sol[1]]    
}

#-----NearestNeighbour function for the TSP----
#obtaining solution by best neighbour (starting from node start) of tsp for a distance matrix
#returns a list wiht path obtained and total distance

NearestNeighbour <- function(d, start){
    
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
    
    distance <- DistanceTSP(d, sol)
    return(list(sol=sol, obj=distance))
    
}


#-----simulated annealing for the TSP-----

SimulatedAnnealingTSP <- function(d, tmax, initial.random=TRUE){
    
    #determining problem size
    n <- dim(d)[1]
    
    #sarting solution
    if(initial.random==TRUE){
        sol <- sample(1:n, n, replace=FALSE)
    }else{
        sol <- NearestNeighbour(d, 1)$sol
    }
    
    sol.best <- sol
    obj <- DistanceTSP(d, sol)
    obj.best <- obj
    
    #tmax number of iterations
    count <- tmax
    
    repeat{
        #select a random move
        pos <- sample(1:n, 1)
        
        #finding new solution and total distance
        sol.test <- SwapNeighbourNodes(sol, pos)
        obj.test <- DistanceTSP(d, sol.test)
        
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

#------A genetic algorithm for the TSP------ 

GeneticAlgorithmTSP <- function(d, npop, iterations, pmut, seed=FALSE){
    
    n <- dim(d)[1]
    
    #initial population
    
    population.father <- replicate(npop, sample(1:n,size=n), simplify=FALSE)
    
    population.son <- replicate(npop, numeric(n), simplify=FALSE)
    
    if(seed == TRUE){
        population.father[[1]] <- NearestNeighbour(d, 1)$sol
    }
    
    sol.best <- numeric(n)
    obj.best <- Inf
    
    fitness <- c(numeric(), npop)
    fitness.acum <- c(numeric(), npop)
    
    iter <- 0
    
    repeat{
        
        #computing relative fitness for selection
        
        for(i in 1:npop)
            fitness[i] <- DistanceTSP(d, population.father[[i]])
        
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
            #selecting father and son
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
            
            population.son[[k]] <- c(c1, c2)
            
            #mutation operator
            
            if(runif(1) < pmut){
                swap <- sample(1:n, 2)
                aux <- population.son[[k]][swap[1]]
                population.son[[k]][swap[1]] <- population.son[[k]][swap[2]]
                population.son[[k]][swap[2]] <- aux
                
            }
        }
                
        population.father <- population.son
        
        iter <- iter + 1
        if(iter >= iterations) break
    }
    
    for(i in 1:npop)
        fitness[i] <- DistanceTSP(d, population.father[[i]])
    
    if(min(fitness) < obj.best){
        obj.best <- min(fitness)
        sol.best <- population.father[[which.min(fitness)]]
    }
    
    return(list(sol = sol.best, obj = obj.best))
}


#-----A tabu search algorithm for the TSP-----  

TabuSearchTSP <- function(d, iter, tabu.size = 7, initial.random = TRUE){
    
    #determining problem size
    n <- dim(d)[1]
    
    #sarting solution
    if(initial.random==TRUE){
        sol <- sample(1:n, n, replace=FALSE)
    }else{
        sol <- NearestNeighbour(d, 1)$sol
    }
    
    sol.best <- sol
    obj <- DistanceTSP(d, sol)
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
                obj.test <- DistanceTSP(d, sol.test)
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


