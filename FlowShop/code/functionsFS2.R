#----NEH heuristic---

#----makespan of the flowshop problem----

makespan <- function(M, sol){
  
  m <- dim(M)[1]
  n <- length(sol)
  
  times <- matrix(numeric(m*n), m, n)
  
  M <- M[ , sol]
  
  times[1,1] <- M[1,1]
  
  for(j in 2:n) times[1, j] <- times[1, (j-1)] + M[1, j]
  for(i in 2:m) times[i, 1] <- times[(i-1), 1] + M[i, 1]
  
  for(i in 2:m){
    for(j in 2:n)
      times[i,j] <- max(times[i-1,j], times[i, j-1]) + M[i, j]
  }
  
  result <- times[m, n]
  return(result)
}

#---- genetic algorithm for the flowshop problem ----

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

#order crossover (OX)
crossoverOX <- function(sol1, sol2){
  
  n <- length(sol1)
  sol <- numeric(n)
  
  cut <- sample(1:(n-1), 1)
  
  sol[1:cut] <- sol1[1:cut]
  sol[(cut+1):n] <- sol2[-which(sol2  %in% sol1[1:cut])]
  
  return(sol)
}

#shift mutation operator (Murata)
shift <- function(sol){
  
  n <- length(sol)
  points <- sample(1:n,2)
  i <- points[1]
  j <- points[2]

  aux <- sol[i]
    
  if(i < j)    
    sol[i:(j-1)] <- sol[(i+1):j]
  else
   sol[(j+1):i] <- sol[j:(i-1)]

  sol[j] <- aux
  
  return(sol)
}

#genetic algorithm
GAFS <- function(problem, npop=10, iter=10, pmut=0.8, crOX=FALSE, elitist=TRUE, alpha=0, verbose=FALSE){
  
  m <- dim(problem)[1]
  n <- dim(problem)[2]
  
  bestfit <- Inf
  bestsol <- numeric(n)
  
  #initializing generation G and sons S
  G <- matrix(numeric(n*npop), npop, n)
  S <- matrix(numeric(n*npop), npop, n)
  for(i in 1:npop) G[i, ] <- sample(1:n, n)
  
  #compute makespan of G
  fit <- apply(G, 1, function(x) makespan(problem, x))
  
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
      
      if(pmut > runif(1)) S[i, ] <- shift(S[i, ])
      
    }
    
    if (elitist) S[sample(1:npop,1), ] <- bestsol
    
    G <- S
    
    #compute makespan of G
    fit <- apply(G, 1, function(x) makespan(problem, x))
    testfit <- min(fit)
    
    #reset counter if better solution found
    if(testfit < bestfit)
      T <- 0
    
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

#----shift moves for the FS----

shiftmove <- function(sol, i, j){
  
  n <- length(sol)

  if(i < j){
    
    aux <- sol[i]
    sol[i:(j-1)] <- sol[(i+1):j]
    sol[j] <- aux
    
  }else{
    
    aux <- sol[i]
    sol[(j+1):i] <- sol[j:(i-1)]
    sol[j] <- aux
    
  }
  
  return(sol)
}


#-----NEH heuristic-----

NEH <- function(Instance){
  
  n <- dim(Instance)[2]
  
  jobs <- order(apply(Instance, 2, sum), decreasing = TRUE)
  print(jobs)
  if( makespan(Instance, c(jobs[1], jobs[2])) < makespan(Instance, c(jobs[2], jobs[1])) )
    sol <- c(jobs[1], jobs[2])
  else
    sol <-c(jobs[2], jobs[1])
  
  print(sol)
  
  for(k in 3:n){
    
    fit <- Inf
    pos <- -1
    
    for(i in 1:k){
      if(i==1) test <- c(jobs[k], sol)
      if(i!= 1 & i!=k) test <- c(sol[1:(i-1)], jobs[k], sol[i:(k-1)])
      if(i==k) test <- c(sol, jobs[k])
      
      testfit <- makespan(Instance, test)
      
      if(testfit < fit){
        fit <- testfit
        pos <- i
      }
      
    }
    
    if(pos==1) sol <- c(jobs[k], sol)
    if(pos!= 1 & pos!=k) sol <- c(sol[1:(pos-1)], jobs[k], sol[pos:(k-1)])
    if(pos==k) sol <- c(sol, jobs[k])
    
    print(pos)
    print(sol)
  }
  
  fit <- makespan(Instance, sol)
  
  return(list(sol=sol, fit=fit))
  
  
}
#---- tabu search ----

#shift operator (Murata)
shift4move <- function(sol, i, j){
  
  aux <- sol[i]  
  
  if(i < j)
    sol[i:(j-1)] <- sol[(i+1):j]
  else
    sol[(j+1):i] <- sol[j:(i-1)]
    
  sol[j] <- aux
    
  return(sol)
}

#TS using the insertion operator
TSFS <- function(Instance, inisol, iter=100, tabu.size=5, eval=FALSE, asp=TRUE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(iter)
    evalbest <- numeric(iter)
  }
  
  #initialization
  sol <- inisol
  bestsol <- inisol
  bestfit <- makespan(Instance, sol)
  
  T <- 1
  Ttabu <- 1
  flag.tabu <- FALSE
  tabu.list <- matrix(numeric(2*tabu.size), tabu.size, 2)
  
  n <- length(sol)
  
  #generating list of moves
  
  moves <- matrix(numeric(2*(n-1)^2), (n-1)^2, 2)
  k <- 1
  for(i in 1:n)
    for(j in 1:n)
      if(i!=j & (i+1)!=j){
        moves[k, ] <- c(i,j)
        k <- k+1
      } 
  
    
  while (T<=iter){
    
    #find the best move
    fit <- Inf
    bestmove <- c(0,0)
    
    for(k in 1:dim(moves)[1]){
      
      testsol <- shift4move(sol, moves[k,1], moves[k,2])
      testfit <- makespan(Instance, testsol)
      
      if(testfit <= fit & sum(moves[k, ]%in% tabu.list)==0){
        fit <- testfit
        bestmove <- moves[k, ]
        flag.tabu <- TRUE
        }
      
      if(testfit < bestfit & sum(moves[k, ]%in% tabu.list)==2 & asp){
        fit <- testfit
        bestmove <- moves[k, ]
        flag.tabu <- FALSE
      }
    }
    
    print(bestmove)
    print(fit)
    
    #update sol
    sol <- shift4move(sol, bestmove[1], bestmove[2])
    print(sol)
    
    #update bestsol
    
    if(fit <= bestfit){
      bestfit <- fit
      bestsol <- sol
      
    }
    
    #update tabu list
    if(flag.tabu){
      tabu.list[Ttabu%%tabu.size+1, ] <- bestmove
      Ttabu <- Ttabu + 1
    }
    
    if(eval){
      evalfit[T] <- fit
      evalbest[T] <- bestfit
    }
    
    T <- T+1
  }
  
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest))
  else
    return(list(sol=bestsol, fit=bestfit))
}

swap4move <- function(sol, i, j){
  
  aux <- sol[i]
  sol[i] <- sol[j]
  sol[j] <- aux
  
  return(sol)
  
}

#TS using the swap operator
TSFS2 <- function(Instance, inisol, iter=100, tabu.size=5, eval=FALSE, asp=TRUE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(iter)
    evalbest <- numeric(iter)
  }
  
  #initialization
  sol <- inisol
  bestsol <- inisol
  bestfit <- makespan(Instance, sol)
  
  T <- 1
  Ttabu <- 1
  flag.tabu <- FALSE
  tabu.list <- matrix(numeric(2*tabu.size), tabu.size, 2)
  
  n <- length(sol)
  
  #generating list of moves
  
  moves <- matrix(numeric(n*(n-1)), n*(n-1)/2, 2)
  k <- 1
  for(i in 1:(n-1))
    for(j in (i+1):n){
        moves[k, ] <- c(i,j)
        k <- k+1
    }
  
  
  while (T<=iter){
    
    #find the best move
    fit <- Inf
    bestmove <- c(0,0)
    
    for(k in 1:dim(moves)[1]){
      
      testsol <- swap4move(sol, moves[k,1], moves[k,2])
      testfit <- makespan(Instance, testsol)
      
      if(testfit <= fit & sum(moves[k, ]%in% tabu.list)==0){
        fit <- testfit
        bestmove <- moves[k, ]
        flag.tabu <- TRUE
      }
      
      if(testfit < bestfit & sum(moves[k, ]%in% tabu.list)==2 & asp){
        fit <- testfit
        bestmove <- moves[k, ]
        flag.tabu <- FALSE
      }
    }
    
    print(bestmove)
    print(fit)
    
    #update sol
    sol <- shift4move(sol, bestmove[1], bestmove[2])
    print(sol)
    
    #update bestsol
    
    if(fit <= bestfit){
      bestfit <- fit
      bestsol <- sol
      
    }
    
    #update tabu list
    if(flag.tabu){
      tabu.list[Ttabu%%tabu.size+1, ] <- bestmove
      Ttabu <- Ttabu + 1
    }
    
    if(eval){
      evalfit[T] <- fit
      evalbest[T] <- bestfit
    }
    
    T <- T+1
  }
  
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest))
  else
    return(list(sol=bestsol, fit=bestfit))
}

#---- simulated annealing for the TSP based on 2opt----

SAFS <- function(Instance, inisol, Tmax=1000, mu=1, eval=FALSE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(Tmax)
    evalbest <- numeric(Tmax)
    evaltest <- numeric(Tmax)
  }
  
  #initialization
  sol <- inisol
  bestsol <- inisol
  fit <- makespan(Instance, sol)
  bestfit <- fit
  T <- Tmax
  n <- length(inisol)
  num.moves <- (n-1)^2
  
  #generating list of moves
  
  moves <- matrix(numeric(2*(n-1)^2), (n-1)^2, 2)
  k <- 1
  for(i in 1:n)
    for(j in 1:n)
      if(i!=j & (i+1)!=j){
        moves[k, ] <- c(i,j)
        k <- k+1
      } 
  
  while (T > 0) {
    
    test <- sample(1:num.moves, 1)
    i <- moves[test, 1]
    j <- moves[test, 2]
    
    testsol <- shift4move(sol, i, j)
    testfit <- makespan(Instance, testsol)
    
    if(testfit <= fit){
      
      sol <- testsol
      fit <- testfit
      
      if(testfit <= bestfit){
        bestsol <- sol
        bestfit <- testfit
      } 
      
    }else{
      if(exp(-mu*(testfit-fit)/T) > runif(1)){
        sol <- testsol
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

