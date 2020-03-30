#------ objective function: makespan -----

makespan <- function(M, sol){
  
  m <- dim(M)[1]
  n <- dim(M)[2]
  
  M <- M[ , sol]
  
  for(j in 2:n) M[1, j] <- M[1, (j-1)] + M[1, j]
  for(i in 2:m) M[i, 1] <- M[(i-1), 1] + M[i, 1]
  
  for(i in 2:m)
    for(j in 2:n)
      M[i,j] <- max(M[i-1,j], M[i, j-1]) + M[i, j]
  
  
  return(M[m, n])
}

#---- constructive heuristics: Palmer, trapezes and Johnson (m=2) ----

PalmerTrapezes <- function(M){
  
  m <- dim(M)[1]
  n <- dim(M)[2]
  S <- matrix(0, 2, n)
  
  for(j in 1:n)
    for(i in 1:m){
      S[1, j] <- S[1, j] + (m - i)*M[i, j]
      S[2, j] <- S[2, j] + (i - 1)*M[i, j]
    }
  
  palmer <- order(S[1, ] - S[2, ])
  trapezes <- Johnson(S)
  
  return(list(pal=palmer, tra=trapezes))
}

Johnson <- function(M){
  
  m <- dim(M)[1]
  n <- dim(M)[2]
  
  if(m!=2) return(success=FALSE)
  
  values <- apply(M, 2, min)
  
  mask <- M[1, ] == values
  
  sol <- numeric(n)
  sequence <- order(values)
  head <- 1
  tail <- n
  
  
  for(i in 1:n){
    if(mask[sequence[i]]==TRUE){
      sol[head] <- sequence[i]
      head <- head + 1
    }else{
      sol[tail] <- sequence[i]
      tail <- tail - 1
    }
  }
  
  return(sol)
  
}

#----- swap of two elements i and j -----

swap <- function(v, i , j){
  
  aux <- v[i]
  v[i] <- v[j]
  v[j] <- aux
  return(v)
}

#---- insertion of element i in position j -----

insertion <- function(v, i, j){
  
  n <- length(v)
  
  if(i!=j){
    aux <- v[i]
    if(i < j)    
      v[i:(j-1)] <- v[(i+1):j]
    else
      v[(j+1):i] <- v[j:(i-1)]
    v[j] <- aux
  }
  return(v)
}

#---- Hill climbing heuristic with swap and insertion -----

HCFS <- function(M, inisol, op="swap"){
  
  if(!op %in% c("swap", "insertion"))
    return(success==FALSE)
  
  n <- length(inisol)
  sol <- inisol
  obj <- makespan(M, sol)
  k <- TRUE
  
  while(k){
    #obtaining the solution of minimum value of o.f. from neighbourhood
    soltest <- numeric(n)
    objtest <- Inf
    
    if(op=="swap"){
      for(i in 1:(n-1))
        for(j in (i+1):n){
          test <- makespan(M, swap(sol, i, j))
          if(test < objtest){
            objtest <- test
            soltest <- swap(sol, i, j)
          }
        }
      
    }
    
    if(op=="insertion"){
      for(i in 1:n)
        for(j in 1:n)
          if(i!=j & i!=(j-1)){
            test <- makespan(M, insertion(sol, i, j))
            if(test < objtest){
              objtest <- test
              soltest <- insertion(sol, i, j)
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

#---- simulated annealing with swap and insertion -----

SAFS <- function(Instance, inisol, Tmax=1000, mu=1, op="swap", eval=FALSE){
  
  #checking operator
  if(!op %in% c("swap", "insertion"))
    return(success==FALSE)
  
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
  
  while (T > 0) {
    
    move <- sample(1:n, 2)
    
    if(op=="swap")
      testsol <- swap(sol, move[1], move[2])
    
    if(op=="insertion")
      testsol <- insertion(sol, move[1], move[2])
    
    testfit <- makespan(Instance, testsol)
    
    if(exp(-mu*(testfit-fit)/T) > runif(1)){
      sol <- testsol
      fit <- testfit
    }
    
    if(testfit <= bestfit){
      bestsol <- testsol
      bestfit <- testfit
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
    return(list(sol=bestsol, obj=bestfit))
}

#---- tabu search for the permutative flow shop ----

# https://stackoverflow.com/questions/32640682/check-whether-matrix-rows-equal-a-vector-in-r-vectorized

TSFS <- function(Instance, inisol, iter=100, tabu.size=5, op="swap", asp=TRUE, eval=FALSE, early=FALSE){
  
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
  
  check_tabu <- function(M, v) any(colSums(apply(M, 1, function(x) v == x))==2)
  
  n <- length(sol)
  
  
  while (T<=iter){
    
    found_best <- FALSE
    
    #find the best move
    fit <- Inf
    bestmove <- c(0,0)
    
    if(op=="swap"){
      
      #examining all swap moves
      for(i in 1:(n-1))
        for(j in (i+1):n){
          testsol <- swap(sol, i, j)
          testfit <- makespan(Instance, testsol)
          
          #improvement with non-tabu move
          if(testfit <= fit & check_tabu(tabu.list, c(i,j))==FALSE){
            fit <- testfit
            bestmove <- c(i, j)
            flag.tabu <- TRUE
          }
          
          #improvement with tabu move, and aspiration condition
          if(testfit < bestfit & check_tabu(tabu.list, c(i,j))==TRUE & asp==TRUE){
            fit <- testfit
            bestmove <- c(i, j)
            flag.tabu <- FALSE
          }
        }
      
      #obtain sol
      sol <- swap(sol, bestmove[1], bestmove[2])
      
    }
    
    
    if(op=="insertion"){
      
      for(i in 1:n)
        for(j in 1:n)
          if(i!=j & i!=(j-1)){
            
            testsol <- insertion(sol, i, j)
            testfit <- makespan(Instance, testsol)
            
            #improvement with non-tabu move
            if(testfit <= fit & check_tabu(tabu.list, c(i,j))==FALSE){
              fit <- testfit
              bestmove <- c(i, j)
              flag.tabu <- TRUE
            }
            
            #improvement with tabu move, and aspiration condition
            if(testfit < bestfit & check_tabu(tabu.list, c(i,j))==TRUE & asp==TRUE){
              fit <- testfit
              bestmove <- c(i, j)
              flag.tabu <- FALSE
            }
          }
      
      #obtain sol
      sol <- insertion(sol, bestmove[1], bestmove[2])
      
    }
    
    #update bestsol
    if(fit < bestfit){
      bestfit <- fit
      bestsol <- sol
      found_best <- TRUE
    }
    
    #update tabu list (works the same for both moves)
    if(flag.tabu){
      if(op=="swap")
        tabu.list[Ttabu%%tabu.size+1, ] <- bestmove
      
      if(op=="insertion")
        tabu.list[Ttabu%%tabu.size+1, ] <- rev(bestmove)
      
      Ttabu <- Ttabu + 1
    }
    
    if(eval){
      evalfit[T] <- fit
      evalbest[T] <- bestfit
    }
    
    
    if(early & found_best)
      T <- 0
    else
      T <- T + 1
  }
  
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest))
  else
    return(list(sol=bestsol, obj=bestfit))
}

plotSA <- function(SA, type="all"){
  library(tidyverse)
  if(type=="all")
    data <- data.frame(t=1:length(SA$evalfit), fit=SA$evalfit, test=SA$evaltest, best=SA$evalbest)
  if(type=="fit")
    data <- data.frame(t=1:length(SA$evalfit), fit=SA$evalfit, best=SA$evalbest)
  data %>% pivot_longer(-t, names_to = "sol", values_to = "fit") %>% ggplot(aes(t, fit, col=sol)) + geom_line()
}

plotTS <- function(TS){
  library(tidyverse)
  data.frame(t=1:length(TS$evalfit), fit=TS$evalfit, best=TS$evalbest) %>% pivot_longer(-t, names_to = "sol", values_to = "fit") %>% ggplot(aes(t, fit, col=sol)) + geom_line()
}

#----- Generator of starting solutions for a GRASP heuristic ----

# generator of starting solutions for a GRASP heuristic
PalmerGenerator <- function(M, rcl=4){
  
  m <- dim(M)[1]
  n <- dim(M)[2]
  S <- matrix(0, 2, n)
  
  for(j in 1:n)
    for(i in 1:m){
      S[1, j] <- S[1, j] + (m - i)*M[i, j]
      S[2, j] <- S[2, j] + (i - 1)*M[i, j]
    }
  
  
  diff <- S[1, ] - S[2, ]
  palmer <- order(diff)
  
  grasp <- numeric(n)
  
  chosen <- rep(FALSE, n)
  
  for(i in 1:n){
    k <- sample(1:min(rcl, n-i+1), 1)
    j <- which(chosen==FALSE)[k]
    grasp[i] <- palmer[j]
    chosen[j] <- TRUE
  }
  
  return(grasp)
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


#---- GRASP for the flowshop problem -----

GRASPFS <- function(M, rcl=4, iter=100, op="swap", opt="HC", ...){
  
  params <- list(...)
  
  if(!opt %in% c("HC", "TS", "SA"))
    return(success==FALSE)
  
  n <- dim(M)[2]
  bestsol <- numeric(n)
  bestfit <- Inf
  
  for(t in 1:iter){
    
    seed_sol <- PalmerGenerator(M, rcl=rcl)
    
    if(opt=="HC")
      test_sol <- HCFS(M, seed_sol, op=op)
    
    if(opt=="TS")
      test_sol <- TSFS(M, seed_sol, iter=25, op=op, early = TRUE)
    
    if(opt=="SA"){
      if(is.null(params$Tmax) & is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, op=op)
      if(!is.null(params$Tmax) & is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, Tmax=params$Tmax, op=op)
      if(is.null(params$Tmax) & !is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, mu=params$mu, op=op)
      if(!is.null(params$Tmax) & !is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, Tmax=params$Tmax, mu=params$mu, op=op)
    }
      
    
    if(test_sol$obj < bestfit){
      bestsol <- test_sol$sol
      bestfit <- test_sol$obj
    }
    
  }
  
  return(list(sol=bestsol, obj=bestfit))
}

#---- Perturbation of a solution for a ILS heuristic ----

# perturbation of a solution for a ILS heuristic
PerturbationInsertion <- function(v, ni=4){
  n <- length(v)
  for(i in 1:ni){
    ch <- sample(1:n, 2)
    v <- insertion(v, ch[1], ch[2])
  }
  return(v)
}

#----- Iterated local search for the flowshop ----

ILSFS <- function(M, ni=4, iter=100, op="swap", opt="HC", ...){
  
  params <- list(...)
  
  if(!opt %in% c("HC", "TS", "SA"))
    return(success==FALSE)
  
  n <- dim(M)[2]
  bestsol <- numeric(n)
  bestfit <- Inf
  
  sol <- PalmerTrapezes(M)$pal
  fit <- makespan(M, sol)
  
  for(i in 1:iter){
    
    seed_sol <- PerturbationInsertion(sol)
    
    if(opt=="HC")
      test_sol <- HCFS(M, seed_sol, op=op)
    
    if(opt=="TS")
      test_sol <- TSFS(M, seed_sol, iter=25, op=op, early = TRUE)
    
    if(opt=="SA"){
      if(is.null(params$Tmax) & is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, op=op)
      if(!is.null(params$Tmax) & is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, Tmax=params$Tmax, op=op)
      if(is.null(params$Tmax) & !is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, mu=params$mu, op=op)
      if(!is.null(params$Tmax) & !is.null(params$mu))
        test_sol <- SAFS(M, seed_sol, Tmax=params$Tmax, mu=params$mu, op=op)
    }
    
    if(test_sol$obj < fit){
      sol <- test_sol$sol
      fit <- test_sol$obj
    }
    
    if(test_sol$obj < bestfit){
      bestsol <- test_sol$sol
      bestfit <- test_sol$obj
    }
  }
  
  return(list(sol=bestsol, obj=bestfit))
  
}
