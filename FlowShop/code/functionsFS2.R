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

