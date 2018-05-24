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

#shift operator (Murata)
shift <- function(sol){
  
  n <- length(sol)
  points <- sample(1:n,2)
  i <- points[1]
  j <- points[2]
  
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

GAFS <- function(problem, npop=10, iter=100, pmut=0.8, elitist=TRUE, verbose=FALSE){
  
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
    probs <- (maxfit - fit)^2
    probs <- probs/sum(probs)
    
    if(verbose) print(probs)
    
    for(i in 1: npop){
      
      select <- sample(1:npop, 2, prob=probs)
      
      S[i, ] <- crossover2pointI(G[select[1], ], G[select[2], ])
      
      if(pmut > runif(1)) S[i, ] <- shift(S[i, ])
      
    }
    
    if (elitist) S[sample(1:npop,1), ] <- bestsol
    
    G <- S
    
    #compute makespan of G
    fit <- apply(G, 1, function(x) makespan(problem, x))
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

NEH <- function(Instance, verbose=FALSE){
  
  m <- dim(Instance)[1]
  n <- dim(Instance)[2]
  
  tasks <- order(apply(Instance, 2, sum), decreasing = TRUE)
  
  if(verbose) cat("ordered tasks:",tasks, "\n")
  
  if(makespan(Instance, c(tasks[1:2])) < makespan(Instance, c(tasks[2:1])))
    sol <- tasks[1:2]
  else
    sol <- tasks[2:1]
  
  if(verbose) cat("starting:", sol, "\n")
  
  for(i in 3:n){
    best <- Inf
    pos <- -1
    
    for(j in 1:i){
      if(j==1) test <- c(tasks[i], sol)
      if(j==i) test <- c(sol, tasks[i])
      if(j!=1 & j!=i) test <- c(sol[1:(j-1)], tasks[i], sol[j:(i-1)])
      
      if(makespan(Instance, test) < best){
        best <- makespan(Instance, test)
        pos <- j
      }
    }
    
    if(pos==1) sol <- c(tasks[i], sol)
    if(pos==i) sol <- c(sol, tasks[i])
    if(pos!=1 & pos!=i) sol <- c(sol[1:(pos-1)], tasks[i], sol[pos:(i-1)])
    
    if(verbose) cat("position:", pos, "\n")
    if(verbose) cat("solution:", sol, "\n")
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

