
#----TSP instance generation from a set of points----

SampleTSP <- function(size, seed=1111, euc=TRUE){
  set.seed(seed)
  points <- matrix(sample(0:100, size*2, replace = TRUE), size, 2)
  d <- dist(points, upper = TRUE)
  d <- as.matrix(d)
  return(d)
}

#----Objective function----

TSP <- function(D, sol){
  
    n <- length(sol)
    d <- 0
    for(i in 1:(n-1)) d <- d + D[sol[i],sol[i+1]]
    d <- d + D[sol[n],sol[1]] 
    return(d)
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
  return(list(sol=sol, fit=distance))
  
}

#----- a 2-opt move for the TSP ----

move2opt <- function(sol, i, k){
  n <- length(sol)
  move <- c(sol[1:(i-1)], rev(sol[i:k]),  sol[(k+1):n])
  
  if(i==1){
    move <- c(rev(sol[i:k]), sol[(k+1):n])
    move <- c(move[k:n], move[1:(k-1)])
  }
  return(move)
}


#----- a hill climbing heuristic based on 2-opt ----

HillClimbing2opt <- function(D, sini, verbose=FALSE){
  
  #compute objective function
  sol <- sini
  obj <- TSP(D, sol)
  n <- dim(D)[1]
  
  change <- TRUE
  
  while(change){
    
    #selecting best 2opt move
    gain <- 0
    igain <- 0
    kgain <- 0
    for(i in 1:(n-2)){
      for(k in (i+1):min(n+i-3,n-1)){
        if(i==1)
          testgain <- D[sol[n], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[n] ,sol[i]] - D[sol[k], sol[k+1]]
        else
          testgain <- D[sol[i-1], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[i-1] ,sol[i]] - D[sol[k], sol[k+1]]
        
        if(testgain <= gain){
          gain <- testgain
          igain <- i
          kgain <- k
        }
      }
    }
    
    if(verbose) print(gain)
    
    #if gain is negative, update
    if(gain < 0){
      sol <- move2opt(sol, igain, kgain)
      obj <- obj  + gain
    }else change <- FALSE
    
  }
  return(list(sol=sol, obj=obj))
}


#---- tabu search based on 2opt -----

TSTSP2opt <- function(D, inisol, iter=100, tabu.size=5, eval=FALSE, asp=TRUE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(iter)
    evalbest <- numeric(iter)
    evalgain <- numeric(iter)
  }
  
  #initialization
  sol <- inisol
  bestsol <- inisol
  fit <- TSP(D, sol)
  bestfit <- fit
  T <- 1
  Ttabu <- 1
  flag.tabu <- FALSE
  tabu.list <- matrix(numeric(2*tabu.size), tabu.size, 2)
  
  n <- dim(D)[1]
  
  while (T<=iter){
    
    #find the best move
    gain <- Inf
    bestmove <- c(0,0)
    
    for(i in 1:(n-2)){
      for(k in (i+1):min(n+i-3,n-1)){
        if(i==1)
          testgain <- D[sol[n], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[n] ,sol[i]] - D[sol[k], sol[k+1]]
        else
          testgain <- D[sol[i-1], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[i-1] ,sol[i]] - D[sol[k], sol[k+1]]
      
        if(testgain < gain & sum(c(i, k) %in% tabu.list)==0){
          gain <- testgain
          bestmove <- c(i, k)
          flag.tabu <- TRUE
        }
        
        if(fit + testgain <= bestfit & sum(c(i, k) %in% tabu.list)==2 & asp){
          gain <- testgain
          bestmove <- c(i, k)
          flag.tabu <- FALSE
        }
      }
    }
    
    #update sol and bestsol
    sol <- move2opt(sol, bestmove[1], bestmove[2])
    fit <- fit + gain
    if(fit <= bestfit){
      bestsol <- sol
      bestfit <- fit
    }
    
    #update tabu list
    
    if(flag.tabu){
      Ttabu <- Ttabu + 1
      tabu.list[Ttabu%%tabu.size+1, ] <- bestmove
    }
    
    if(eval){
      evalfit[T] <- fit
      evalbest[T] <- bestfit
      evalgain[T] <- gain
    }
    
    T <- T+1
  }
  
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest, evalgain=evalgain))
  else
    return(list(sol=bestsol, fit=bestfit))
}

#---- a silly search strategy to illustrate the need of tabu list -----

SillyTSP2opt <- function(D, inisol, iter=100){
  
  #tracking evaluation

  evalfit <- numeric(iter)
  evalbest <- numeric(iter)
  evalgain <- numeric(iter)

  #initialization
  sol <- inisol
  bestsol <- inisol
  fit <- TSP(D, sol)
  bestfit <- fit
  T <- 1

  n <- dim(D)[1]
  
  while (T<=iter){
    
    #find the best move
    gain <- Inf
    bestmove <- c(0,0)
    
    for(i in 1:(n-2)){
      for(k in (i+1):min(n+i-3,n-1)){
        if(i==1)
          testgain <- D[sol[n], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[n] ,sol[i]] - D[sol[k], sol[k+1]]
        else
          testgain <- D[sol[i-1], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[i-1] ,sol[i]] - D[sol[k], sol[k+1]]
        
        if(testgain < gain){
          gain <- testgain
          bestmove <- c(i, k)
        }
      }
    }
    
    #update sol and bestsol
    sol <- move2opt(sol, bestmove[1], bestmove[2])
    fit <- fit + gain
    if(fit <= bestfit){
      bestsol <- sol
      bestfit <- fit
    }
    
    evalfit[T] <- fit
    evalbest[T] <- bestfit
    evalgain[T] <- gain

    T <- T+1
  }
  
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest, evalgain=evalgain))
}

#---- simulated annealing for the TSP based on 2opt----

SATSP2opt <- function(D, inisol, Tmax=1000, mu=1, eval=FALSE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(Tmax)
    evalbest <- numeric(Tmax)
    evaltest <- numeric(Tmax)
  }
  
  #initialization
  sol <- inisol
  bestsol <- inisol
  fit <- TSP(D, sol)
  bestfit <- fit
  T <- Tmax
  n <- dim(D)[1]
  num.moves <- n*(n-3)/2
  
  #generating possible moves
  moves <- matrix(numeric(num.moves*2), num.moves, 2)
  
  count <- 1
  for(i in 1:(n-2)){
    for(k in (i+1):min(n+i-3,n-1)){
      moves[count, ] <- c(i, k)
      count <- count+1
    }
  }
  
  while (T > 0) {
    
    test <- sample(1:num.moves, 1)
    i <- moves[test, 1]
    k <- moves[test, 2]
    
    if(i==1)
      testgain <- D[sol[n], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[n] ,sol[i]] - D[sol[k], sol[k+1]]
    else
      testgain <- D[sol[i-1], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[i-1] ,sol[i]] - D[sol[k], sol[k+1]]
    
    testfit <- fit + testgain
    
    if(testfit <= fit){
      
      sol <- move2opt(sol, i, k)
      fit <- testfit
      
      if(testfit <= bestfit){
        bestsol <- sol
        bestfit <- testfit
      } 
      
    }else{
      if(exp(-mu*(testfit-fit)) > runif(1)){
        sol <- move2opt(sol, i, k)
        fit <- testfit
      }
    }
    
    if(eval){
      evalfit[Tmax - T + 1] <- fit
      evalbest[Tmax - T + 1] <- bestfit
      evaltest[Tmax - T + 1] <- testfit
    }
    
    T <- T - 1
  }
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest, evaltest=evaltest))
  else
    return(list(sol=bestsol, fit=bestfit))
}
 
#---- A GRASP heuristic for the TSP ----

GRASP <- function(D, rcl.size, start=1){
  
  n <- dim(D)[1]
  
  #test if d is a square matrix
  #test if start <= n
  
  for(i in 1:n)
    D[i, i] <- Inf
  
  #flag checks if node has been included in solution
  #sol is the solution starting from start
  
  flag <- logical(n)
  sol <- integer(n)
  
  sol[1] <- start
  flag[start] <- TRUE
  
  for(i in 1:(n-1)){
    
    row <- D[sol[i], ]
    pos <- which(flag==FALSE)
    rcl <- pos[order(row[pos], decreasing = FALSE)]
    lengthrcl <- min(rcl.size, n-sum(flag))
    
    rank <- sample(1:lengthrcl, 1)
    j <- rcl[rank]
    sol[i+1] <- j
    flag[j] <- TRUE
    
  }
  
  distance <- TSP(D, sol)
  return(list(sol=sol, fit=distance))
  
}

GRASP2opt <- function(D, rcl.size, tries, start=1, ls="TS"){
  
  n <- dim(D)[1]
  
  bestsol <- numeric(n)
  bestfit <- Inf
  
  report <- matrix(numeric(2*tries), tries, 2)
  
  for(count in 1:tries){
    
    initial.solution <- GRASP(D, rcl.size)
    
    if(ls=="TS")
        refined.solution <- TSTSP2opt(D, initial.solution$sol, iter=10*n)
      else
        refined.solution <- SATSP2opt(D, initial.solution$sol, Tmax=n*100)
    
    if(refined.solution$fit <= bestfit){
      bestfit <- refined.solution$fit
      bestsol <- refined.solution$sol
    }  

    report[count, 1] <- initial.solution$fit     
    report[count, 2] <- refined.solution$fit
  }
  
  report <- as.data.frame(report)
  names(report) <- c("initial", "refined")
  
  return(list(sol=bestsol, fit=bestfit, report=report))
}
 