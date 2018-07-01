#----Generators of KP instances----

KPInstance01 <- function(R, size, h, seed=1212){
    set.seed(1212)
    w <- sample(1:R, size, replace=TRUE)
    u <- sapply(1:size, function(x) sample((w[x] + R/10 - R/20):(w[x] + R/10 + R/20), 1))
    W <- round(h/(101)*sum(w))
    return(list(w=w,u=u,W=W))
}

SubsetSumInstance01 <- function(R, size, h, seed=1212){
    set.seed(1212)
    w <- sample(1:R, size, replace=TRUE)
    u <- w
    W <- round(h/(101)*sum(w))
    return(list(w=w,u=w,W=W))
}

EasyInstance <- function(size, h){
  w <- sample(50:80, size, replace = TRUE)
  u <- sample(100:160, size, replace = TRUE)
  W <- round(h*sum(w))
  
  return(list(u=u, w=w, W=W))
}

#-----Relative utility heuristic----

HeuristicKP  <- function(Instance){
  
  solution <- rep(FALSE, length(Instance$u))
  RelUtil <- Instance$u/Instance$w
  OrderRelUtil <- order(RelUtil, decreasing = TRUE)
  
  Wleft <- Instance$W
  
  for(i in 1:length(Instance$u)){
    
    pos <- OrderRelUtil[i]
    if(Wleft >= Instance$w[pos]){
      solution[pos] <- TRUE
      Wleft <- Wleft - Instance$w[pos]
    }
  }
  
  fit <- FitnessKP(Instance, solution)
  
  return(list(sol=solution, fit=fit))
}

#---- simulated annealing for the KP ----

SAKP <- function(Inst, inisol, Tmax=1000, mu=1, eval=FALSE){
  
  sol <- inisol
  bestsol <- inisol
  
  if(eval){
    evalfit <- numeric(Tmax)
    evalbest <- numeric(Tmax)
    evaltest <- numeric(Tmax)
  }
  
  fit <- FitnessKP(Inst, sol)
  bestfit <- FitnessKP(Inst, bestsol)
  
  T <- Tmax
  
  while (T > 0) {
    
    testsol <- sol
    pos <- sample(1:length(Inst$u), 1)
    testsol[pos] <- !testsol[pos]
    testfit <- FitnessKP(Inst, testsol)
    
    if(testfit >= fit){
      
      sol <- testsol
      fit <- testfit
      
      if(testfit >= bestfit){
        bestsol <- testsol
        bestfit <- testfit
      } 
      
    }else{
      if(exp(-mu*(fit-testfit)/T) > runif(1)){
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

#----tabu search for KP----

TSKP <- function(Inst, inisol, iter=100, tabu.size=5, eval=FALSE, asp=TRUE){
  
  #tracking evaluation
  if(eval){
    evalfit <- numeric(iter)
    evalbest <- numeric(iter)
    move <- numeric(iter)
  }
  
  #initialization
  sol <- inisol
  bestsol <- inisol
  fit <- FitnessKP(Inst, sol)
  bestfit <- fit
  T <- 1
  Ttabu <- 1
  flag.tabu <- FALSE
  tabu.list <- numeric(tabu.size)
  
  n <- length(Inst$w)
  
  while (T<=iter){
    
    #find the best move
    fit <- -1
    ibest <- 0
    for(i in 1:n){
      
      testsol <- sol
      testsol[i] <- !testsol[i] #bit flip operator
      testfit <- FitnessKP(Inst, testsol)
      
      if(testfit >= bestfit & sum(tabu.list==i)==1 & asp){
        ibest <- i
        fit <- testfit
        flag.tabu <- TRUE
      }
      
      if(testfit >= fit & sum(tabu.list==i)==0){
        ibest <- i
        fit <- testfit
        flag.tabu <- FALSE
      }
    }
    
    #update sol and bestsol
    sol[ibest] <- !sol[ibest]
    if(fit >= bestfit){
      bestsol <- sol
      bestfit <- fit
    }
    
    #update tabu list
    
    if(!flag.tabu){
      Ttabu <- Ttabu + 1
      tabu.list[Ttabu%%tabu.size+1] <- ibest
    }
    
    if(eval){
      evalfit[T] <- fit
      evalbest[T] <- bestfit
      move[T] <- ibest
    }
    
    T <- T+1
  }
  
  if(eval)
    return(list(sol=bestsol, fit=bestfit, evalfit=evalfit, evalbest=evalbest, move=move))
  else
    return(list(sol=bestsol, fit=bestfit))
}



#-----A genetic algorithm for the KP----

FitnessKP <- function(dades, sol){
    if(sum(dades$w*sol) <=dades$W) fit <- sum(dades$u*sol)
    else fit <- 0
    return(fit)
}

CrossoverKP1point <- function(sol1, sol2){
    n <- length(sol1)
    point <- sample(1:(n-1), 1)
    son <- c(sol1[1:point], sol2[(point +1):n])
    
    return(son)
}

MutationKP3points <- function(sol){
    n <- length(sol)
    point <- sample(1:(n-2), 1)
    sol <- as.logical(sol)
    sol[point:(point+2)] <- !sol[point:(point+2)]
    sol <- as.numeric(sol)
    
    return(sol)
}

GeneticAlgorithmKP <- function(dades, npop, iterations, pmut, elitist =TRUE, verbose=FALSE){
    
    #problem size  
    n <- length(dades$w)
    
    #initial population
    rel.size <- min(dades$W/sum(dades$w),1)
    population.parent <- replicate(npop, sample(0:1, n, replace = TRUE, prob=c(1-rel.size, rel.size)), simplify=FALSE)
    
    #initializing next generation
    population.son <- replicate(npop, numeric(n), simplify=FALSE)
    
    best.obj <- 0
    best.sol <- numeric(n)
    
    count <- 1
    
    while(count <= iterations){
        
        #assess fitness of population parent
        
        fitness <- sapply(population.parent, function(x) FitnessKP(dades, x))
        
        #finding best solution
        iter.obj <- max(fitness)
        pos.max <- which(fitness == max(fitness))[1]
        iter.sol <- population.parent[[pos.max]]
        
        #update best solution
        if(iter.obj > best.obj){
            best.sol <- iter.sol
            best.obj <- iter.obj
        }
        
        #elitist strategy: replace worst solution by best solution so far
        if(elitist){
            pos.min <- which(fitness == min(fitness))[1]
            population.son[[pos.min]] <- best.sol 
            fitness[pos.min] <- best.obj
        }
        
        #verbose printing
        if(verbose & (count%%10==0)) {
            print(paste("population", count))
            sapply(1:npop, function(x) print(c(population.parent[[x]], fitness[x]), digits=0))
        }
        
        #computing probability of selection
        prob <- fitness^2
        prob <- prob/sum(prob)
        
        #creating son population
        for(i in 1:npop){
            #crossover
            parents <- sample(1:npop, 2, prob = prob)
            v1 <- population.parent[[parents[1]]]
            v2 <- population.parent[[parents[2]]]
            population.son[[i]] <- CrossoverKP1point(v1, v2)
            #mutation
            if(pmut > runif(1)) population.son[[i]] <- MutationKP3points(population.son[[i]])
        }
        
        
        #population son becomes parent
        population.parent <- population.son
        
        count <- count + 1
    }
    
    return(list(sol=best.sol, obj=best.obj))
    
}

#---solution with linear programming----

LPKP <- function(dades){
    library(Rglpk)
    n <- length(dades$u)
    obj <- dades$u
    mat <- matrix(dades$w, nrow=1)
    types <- rep("B", n)
    sol.lp <- Rglpk_solve_LP(obj=obj, mat=mat,  dir="<=", rhs=dades$W, types=types, max = TRUE)
    return(list(obj = sol.lp$optimum, sol = sol.lp$solution))
}


#---A branch and bound for the KP----

#n: number of items
#elements of a node of the branch and bound:
#cons: number of elements considered
#included: vector of logicals of size n: TRUE if element is in the subset of solutions, FALSE otherwise
#bound: upper bound of the solution
#cal_bound: function thatcalculates the bound for a node


#function that sorts an instance of KP

sort.KP <- function(instance){
    r <- instance$u/instance$w
    
    u <- instance$u[order(r, decreasing = TRUE)]
    w <- instance$w[order(r, decreasing = TRUE)]
    
    return(list(inst=list(u=u, w=w, W=instance$W), ord=order(r, decreasing = TRUE)))
    
}

#Creating a node computing its bound (works with a sorted instance)

nodegen <- function(instance, cons, included){
    
    weight <- sum(instance$w * included)
    
    if(weight > instance$W) return(list(success=FALSE)) #unfeasible node
    
    bound <- sum(instance$u * included)
    size <- length(instance$u)
    
    if(cons==size) #terminal node
        return(list(node=list(cons=cons, included=included, bound=bound), success=TRUE))
    
    loop <- TRUE
    k <- cons+1
    
    while(loop==TRUE & cons < size){
        if((weight + instance$w[k]) < instance$W){
            bound <- bound + instance$u[k]
            weight <- weight + instance$w[k]
            if(k==size)
                loop <- FALSE
            else
                k <- k+1
        }
        else #maximum capacity reached
            loop <- FALSE
    }
    
    if(k < size) #completing weight with fractional size
        bound <- bound + (instance$W - weight)*instance$u[k]/instance$w[k]
    
    return(list(node=list(cons=cons, included=included, bound=bound), success=TRUE))
}

#branch and bound function

BandBKP <- function(instance0, report.evol=FALSE){
    
    #sorting the instance
    instance <- sort.KP(instance0)$inst
    ord <- sort.KP(instance0)$ord
    
    size <- length(instance$w)
    
    #setting the lower bound
    best.obj <- 0
    best.sol <- logical(size)
    
    #initializing the list of nodes with initial node
    
    nodes <- list(nodegen(instance, 0, rep(FALSE, size))$node)
    
    #initializing evolution report
    
    if(report.evol){
      report.bounds <- numeric()
      report.nodes <- numeric()
      iter <- 1
      report.bounds[iter] <- best.obj
      report.nodes[iter] <- length(nodes)
    }
    
    #the Branch and Bound Loop
    
    while(length(nodes) > 0){
        
        #selecting the node to branch
        bounds <- sapply(nodes, function(x) x$bound)
        conss <- sapply(nodes, function(x) x$cons)
        select <- rep(0, length(nodes))
        select[which(bounds==max(bounds))] <- bounds[which(bounds==max(bounds))]
        select <- select*conss
        node.to.branch <- which(select==max(select))[1]
        
        #branching the node
        
        cons.to.branch <- nodes[[node.to.branch]]$cons + 1
        
        includedT <- nodes[[node.to.branch]]$included
        includedT[cons.to.branch] <- TRUE
        
        includedF <- nodes[[node.to.branch]]$included
        
        branch.T <- nodegen(instance, cons.to.branch, includedT)
        branch.F <- nodegen(instance, cons.to.branch, includedF)
        
        if(cons.to.branch == size){ #nodes are terminal
            if(branch.T$success == TRUE){
                node.T <- branch.T$node
                bound.T <- node.T$bound
                if(bound.T > best.obj){
                    best.obj <- bound.T
                    best.sol <- node.T$included
                }
            }
            
            node.F <- branch.F$node
            bound.F <- node.F$bound
            if(bound.F > best.obj){
                best.obj <- bound.F
                best.sol <- node.F$included
            }
            
        }else{ #nodes are not terminal: add to the list if feasible
            if(branch.T$success==TRUE)
                nodes[[length(nodes)+1]] <- branch.T$node
            
            nodes[[length(nodes)+1]] <- branch.F$node
            
            
        }
        
        #pruning the nodes that are worse than current solution
        
        bounds <- sapply(nodes, function(x) x$bound)
        
        pruned <- rep(FALSE, length(nodes))
        pruned[node.to.branch] <- TRUE
        pruned[which(bounds < best.obj)] <- TRUE
        
        
        new.nodes <- list()
        k <- 0
        
        for(i in 1:length(nodes)){
            if(!pruned[i]){
                k <- k+1
                new.nodes[[k]] <- nodes[[i]]
            }
        }
        
        nodes <- new.nodes
        
        if(report.evol){
          iter <- iter+1
          report.bounds[iter] <- best.obj
          report.nodes[iter] <- length(nodes)
        }
    }
    
    #undo the ordering to obtain the right solution
    
    sol.true <- numeric(size)
    
    for(i in 1:size) sol.true[ord[i]] <- best.sol[i]
    
    if(report.evol)
      return(list(solution=list(obj=best.obj, sol=sol.true), report=list(bounds=report.bounds, nodes=report.nodes)))
    else
      return(list(obj=best.obj, sol=sol.true))
    
}

