
#----TSP instance generation from a set of points----

SampleTSP <- function(size, seed=1111, euc=TRUE){
  set.seed(seed)
  points <- matrix(sample(0:100, size*2, replace = TRUE), size, 2)
  d <- dist(points, upper = TRUE)
  d <- as.matrix(d)
  return(list(distances=d, coordinates=points))
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
      if(exp(-mu*(testfit-fit)/T) > runif(1)){
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

GRASP2opt <- function(D, rcl.size, tries, start=1, iter=50, ls="TS"){
  
  n <- dim(D)[1]
  
  bestsol <- numeric(n)
  bestfit <- Inf
  
  report <- matrix(numeric(2*tries), tries, 2)
  
  for(count in 1:tries){
    
    initial.solution <- GRASP(D, rcl.size)
    
    if(ls=="TS")
        refined.solution <- TSTSP2opt(D, initial.solution$sol, iter=iter)
      else
        refined.solution <- SATSP2opt(D, initial.solution$sol, Tmax=n*(n-3)*50)
    
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

#---- iterated local search ----

ILSTSTSP2opt <- function(D, inisol, rounds=10, iter=100, tabu.size=5, eval=TRUE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(rounds)
    evalbest <- numeric(rounds)
  }
  
  
  n <- length(inisol)
  sol <- inisol
  
  tabu <- TSTSP2opt(D, inisol, iter=iter, tabu.size = tabu.size)
  
  sol <- tabu$sol
  fit <- tabu$fit
  
  bestsol <- tabu$sol
  bestfit <- tabu$fit
  
  
  for(count.rounds in 1:rounds){
    
    #perturbated solution: 4-opt move
    
    #select cutoff points
    points.notright <- TRUE
    while(points.notright){
      cp <- sort(sample(2:n, 4))
      diffs <- numeric(3)
      for(i in 1:3) diffs[i] <- cp[i+1] - cp[i]
      points.notright <- sum(c(0,1) %in% diffs)
    }
    
    sol.pert <- sol <- c(sol[1:(cp[1]-1)], sol[cp[3]:(cp[4]-1)], sol[cp[2]:(cp[3]-1)], sol[cp[1]:(cp[2]-1)], sol[cp[4]:n])
    
    tabu <- TSTSP2opt(D, sol.pert, iter=iter, tabu.size = tabu.size)
    
    #acceptance criterion
    if(tabu$fit <= fit){
      sol <- tabu$sol
      fit <- tabu$fit
    }
    
    if(tabu$fit <= bestfit){
      bestsol <- tabu$sol
      bestfit <- tabu$fit
    }
    
    #tracking evaluation
    if(eval){
      evalfit[count.rounds] <- fit
      evalbest[count.rounds] <- bestfit
    }

  }
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest))
  else
  return(list(sol=bestsol, fit=bestfit)) 
}  

#----genetic algorithm for the TSP

GATSP <- function(Instance, npop=10, iter=100, crOX=TRUE, pmut=0.8, memetic=FALSE, elitist=TRUE, alpha=0, verbose=FALSE){
  
  n <- dim(Instance)[1]
  
  #setting moves for 2-opt
  num.moves <- n*(n-3)/2
  moves <- matrix(numeric(num.moves*2), num.moves, 2)
  count <- 1
  for(i in 1:(n-2)){
    for(k in (i+1):min(n+i-3,n-1)){
      moves[count, ] <- c(i, k)
      count <- count+1
    }
  }
  
  
  bestfit <- Inf
  bestsol <- numeric(n)
  
  #initializing generation G and sons S
  G <- matrix(numeric(n*npop), npop, n)
  S <- matrix(numeric(n*npop), npop, n)
  for(i in 1:npop) G[i, ] <- sample(1:n, n)
  
  #compute objective function of G
  fit <- apply(G, 1, function(x) TSP(Instance, x))
  
  #update best solution
  testfit <- min(fit)
  if(testfit <= bestfit){
    bestfit <- testfit
    bestsol <- G[which(fit==testfit)[1], ]
  }
  
  T <- 1
  
  while(T <= iter){
    
    #compute probabilities of selection
    maxfit <- max(fit)
    probs <- ((1+alpha)*maxfit - fit)^2
    probs <- probs/sum(probs)
    
    if(verbose) print(probs)
    
    for(i in 1: npop){
      
      select <- sample(1:npop, 2, prob=probs)
      
      if(crOX)
        S[i, ] <- crossoverOX(G[select[1], ], G[select[2], ])
      else
        S[i, ] <- crossover2pointI(G[select[1], ], G[select[2], ])
      
      if(memetic)
        S[i, ] <- HillClimbing2opt(Instance, S[i, ])$sol
      else{
        if(pmut > runif(1)){
          mv <- sample(1:num.moves, 1)
          S[i, ] <- move2opt(S[i, ], moves[mv, 1], moves[mv, 2])
        }
      }
    }
    
    if (elitist) S[sample(1:npop,1), ] <- bestsol
    
    G <- S
    
    #compute distances of G
    fit <- apply(G, 1, function(x) TSP(Instance, x))
    testfit <- min(fit)
    
    #update best solution
    if(testfit <= bestfit){
      bestfit <- testfit
      bestsol <- G[which(fit==testfit)[1], ]
    }
    
    if(verbose) print(bestfit)
    
    T <- T+1
    
  }
  
  return(list(sol=bestsol, fit=bestfit))
}

#2-point crossover version I (Murata)
crossover2pointI <- function(sol1, sol2){
  
  n <- length(sol1)
  son <- numeric(n)
  points <- sort(sample(1:n,2))
  
  if(points[2]==points[1]+1) return(sol1)
  
  son[1:points[1]] <- sol1[1:points[1]]
  son[points[2]:n] <- sol1[points[2]:n]
  son[(points[1]+1):(points[2]-1)] <- sol2[-which(sol2 %in% sol1[1:points[1]] | sol2 %in% sol1[points[2]:n])]
  
  return(son)
}

#swap mutation operator
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

#order crossover (OX)
crossoverOX <- function(sol1, sol2){
  
  n <- length(sol1)
  sol <- numeric(n)
  
  cut <- sample(1:(n-1), 1)
  
  sol[1:cut] <- sol1[1:cut]
  sol[(cut+1):n] <- sol2[-which(sol2  %in% sol1[1:cut])]
  
  return(sol)
}


#----plotting a solution of TSP----

plotTSP <- function(coord, sol){
  
  #defining coordinates of segments
  
  n <- length(sol)
  
  lines <- data.frame(matrix(numeric(4*n), n, 4))
  colnames(lines) <- c("x1", "y1", "x2", "y2")
  
  for(i in 1:(n-1))
    lines[i, ] <- c(coord[sol[i], ], coord[sol[i+1], ])

  lines[n, ] <- c(coord[sol[n], ], coord[sol[1], ])
  
  colnames(lines) <- c("x1", "y1", "x2", "y2")
  
  plot(coord[, 1] , coord[, 2], xlab="", ylab="", xaxt="n", yaxt="n", pch=19)
  segments(lines$x1, lines$y1, lines$x2, lines$y2)  

}

#----gettin optimal tour from TSPLIB instances----

getOptTour <- function(file){
  
  opttour <- readLines(file)
  
  k <- which()
  
  n <- as.numeric(strsplit(opttour[4], " ")[[1]][3])
  
  tour <- numeric(n)
  
  for(i in 1:n) tour[i] <- opttour[i+5]
  
  tour <- as.numeric(tour)
  
  return(tour)
}

#---- solving TSP with LP using the MZT formulation ----

MZT <- function(D, verbose=FALSE){
  
  n <- dim(D)[1]
  
  xIndex <- rep(1:n, each=n-1)
  aux <- 1:n
  yIndex <- unlist(lapply(1:n, function (x) aux[which(aux!=x)]))
  xIndex <- c(xIndex, rep(0, n-1))
  yIndex <- c(yIndex, rep(0, n-1))
  uIndex <- c(rep(0,n*(n-1)), 2:n)
  
  ncols <- n*(n-1) + (n-1)
  uNrows <- (n-1)*(n-2)
  nrows <- 2*n + uNrows
  
  A <- matrix(numeric(ncols*nrows), nrows, ncols)
  
  for(i in 1:n){
    A[i, which(xIndex==i)] <- 1
    A[n+i, which(yIndex==i)] <- 1
  }
  
  uRows <- matrix(numeric(2*uNrows), uNrows, 2)
  
  k <- 1
  for(i in 2:n){
    for(j in 2:n){
      if(i!=j){
        uRows[k,1] <- i
        uRows[k,2] <- j
        k <- k+1
      }
    }
  }
  
  for(i in 1:uNrows){
    A[2*n+i, which(uIndex==uRows[i,1])] <- 1
    A[2*n+i, which(uIndex==uRows[i,2])] <- -1
    A[2*n+i, which(xIndex==uRows[i,1] & yIndex==uRows[i,2])] <- n
  }
  
  #rhs terms
  rhs <- c(rep(1, 2*n), rep(n-1, uNrows))
  
  #dir terms
  dir <- c(rep("==", 2*n), rep("<=", uNrows))
  
  #costs (distances)
  c <- numeric(ncols)
  
  for(i in 1:(n*(n-1)))
    c[i] <- D[xIndex[i], yIndex[i]]
  
  #variable types
  
  vars <- c(rep("B", n*(n-1)), rep("C", (n-1)))
  
  library(Rglpk)
  solLP <- Rglpk_solve_LP(obj=c, mat=A, dir=dir, rhs=rhs, types=vars, max=FALSE, verbose=verbose)
  
  # print(solLP$solution)
  
  xSol <- xIndex[which(solLP$solution[1:(n*(n-1))]==1)]
  ySol <- yIndex[which(solLP$solution[1:(n*(n-1))]==1)]
  
  sol <- matrix(c(xSol, ySol), n, 2)
  
  route <- numeric(n)
  
  route[1] <- 1
  k <- 1
  for(i in 2:n){
    s <- ySol[which(xSol==k)]
    route[i] <- s
    k <- s
  }
  
  return(list(sol=route, fit=solLP$optimum))
  
}
