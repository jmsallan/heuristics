## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----KP instances function-----------------------------------------------
KPInstance01 <- function(R, size, h, seed=1212){
  set.seed(1212)
  w <- sample(1:R, size, replace=TRUE)
  u <- sapply(1:size, function(x) sample((w[x] + R/10 - R/500):(w[x] + R/10 + R/500), 1))
  W <- round(h/(101)*sum(w))
  return(list(w=w,u=u,W=W))
}

## ----KP fitness----------------------------------------------------------
FitnessKP <- function(dades, sol){
  if(sum(dades$w*sol) <=dades$W) fit <- sum(dades$u*sol)
    else fit <- 0
    return(fit)
}

## ----KP 1-point crossover------------------------------------------------
CrossoverKP1point <- function(sol1, sol2){
  n <- length(sol1)
  point <- sample(1:(n-1), 1)
  son <- c(sol1[1:point], sol2[(point +1):n])
  
  return(son)
}

## ----KP crossover test---------------------------------------------------
set.seed(1213)
a <- rep(1, 10)
b <- rep(0, 10)
for(i in 1:5) print(CrossoverKP1point(a,b))

## ----KP mutation---------------------------------------------------------
MutationKP3points <- function(sol){
  n <- length(sol)
  point <- sample(1:(n-2), 1)
  sol <- as.logical(sol)
  sol[point:(point+2)] <- !sol[point:(point+2)]
  sol <- as.numeric(sol)
  
  return(sol)
}

## ----KP mutation show----------------------------------------------------
ex <- rep(1, 10)
for(i in 1:5) print(MutationKP3points(ex))

## ----gaKP----------------------------------------------------------------
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
    if(verbose & (count%%10==0 | count==1)) {
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

## ----setting instancesJP-------------------------------------------------
InstanceKP10 <- KPInstance01(1000, 10, 40)
InstanceKP50 <- KPInstance01(1000, 50, 40, 2424)
InstanceKP100  <- KPInstance01(1000, 100, 40, 3333)

## ----LPKP----------------------------------------------------------------
LPKP <- function(dades){
  library(Rglpk)
  n <- length(dades$u)
  obj <- dades$u
  mat <- matrix(dades$w, nrow=1)
  types <- rep("B", n)
  sol.lp <- Rglpk_solve_LP(obj=obj, mat=mat,  dir="<=", rhs=dades$W, types=types, max = TRUE)
  return(list(sol = sol.lp$solution, obj = sol.lp$optimum))
}

## ----solve LPKP, message=FALSE-------------------------------------------
solPL.KP10 <- LPKP(InstanceKP10)
solPL.KP50 <- LPKP(InstanceKP50)

## ----runga01-------------------------------------------------------------
set.seed(4444)
solGA.KP10 <- GeneticAlgorithmKP(InstanceKP10, 10, 50, 0, elitist=FALSE, verbose=TRUE)

## ----runga02-------------------------------------------------------------
set.seed(4444)
solGA.KP10 <- GeneticAlgorithmKP(InstanceKP10, 10, 100, 0.2, elitist=TRUE, verbose=TRUE)

## ----runga04-------------------------------------------------------------
set.seed(4444)
solGA.KP10 <- GeneticAlgorithmKP(InstanceKP10, 10, 1000, 0.2, elitist=TRUE, verbose=FALSE)
solGA.KP10
solPL.KP10

## ----runga05-------------------------------------------------------------
set.seed(4444)
solGA.KP50 <- GeneticAlgorithmKP(InstanceKP50, 10, 1000, 0.2, elitist=TRUE, verbose=FALSE)
solGA.KP50
solPL.KP50

## ----runga06-------------------------------------------------------------
set.seed(4444)
solGA.KP100 <- GeneticAlgorithmKP(InstanceKP100, 10, 1000, 0.2, elitist=TRUE, verbose=FALSE)
solGA.KP100

