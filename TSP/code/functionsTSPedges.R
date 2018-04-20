#----from matrix distance to  ordered by distance edge list----

edgeList <- function(d){
  n <- dim(d)[1]
  
  #values of origin and destination nodes
  row <- rep(1:n, n)
  column <- as.vector(sapply(1:n, function(x) rep(x, n)))
  
  edges <- data.frame(as.vector(d), row, column)
  names(edges) <- c("distance", "ori", "des")
  
  edges <- edges[order(edges$distance), ]
  edges <- edges[-which(edges$ori == edges$des), ]
  rownames(edges) <- 1:nrow(edges)
  
  return(edges)
}

#----checks if an edge is compatible with the rest of edges----

#the edgeset list contains the sel node

compatibleEdge <- function(list, edgeset, sel){
  
  size <- max(list$ori) #number of vertices
  
  origin <- list[sel, ]$ori
  destin <- list[sel, ]$des
  #cat("origin ", origin, "\n")
  #cat("destin ", destin, "\n")
  
  #list of existing edges, excluding edge sel
  list.edges <- list[which(edgeset), ]
  list.edges <- list.edges[-which(as.numeric(rownames(list.edges))==sel), ]
  #print(list.edges)
  
  compatible <- TRUE
  check <- TRUE #TRUE if we need to check partial cycle
  
  #check fork
  fork1 <- nrow(list.edges[which(list.edges$ori==origin), ]) == 1
  fork2 <- nrow(list.edges[which(list.edges$des==destin), ]) == 1
  
  #cat("fork1: ", fork1, "fork2: ", fork2, "\n")
  
  if(fork1 | fork2){
    compatible <- FALSE
    check <- FALSE
  }
  
  #check partial cycle
  length.cicle <- 1
  
  while(check){
    destin <- list.edges[which(list.edges$ori==destin), ]$des
    #print(destin)
    
    if(length(destin)==0){
      compatible <- TRUE
      check <- FALSE
    } else if(destin==origin & length.cicle < (size-1)){
      compatible <- FALSE
      check <- FALSE
    } else{
      length.cicle <- length.cicle + 1
    }
  }
  
  return(compatible)
}

#----creates a node, checking if it is feasible and computes bound----

nodegen <- function(list, edgeset, cons){
  
  compatible <- compatibleEdge(list, edgeset, cons)
  
  if(!compatible)
    return(list(feasible=FALSE)) #unfeasible node
  
  size <- max(list$ori) #number of vertices
  
  edges.needed <- size - sum(edgeset)
  edges.available <- length(edgeset) - cons
  
  if(edges.needed > edges.available)
    return(list(feasible=FALSE)) #unfeasible node
  
  bound <- sum(edgeset*list$distance)
  
  if(edges.needed > 0)
    bound <- bound + sum(list$distance[(cons+1):(cons+edges.needed)])
  
  return(list(node=list(cons=cons, edgeset=edgeset, bound=bound), feasible=TRUE))
  
}

#----distances heuristic for the TSP------

HeuristicDistances <- function(d){
  
  #instance size (number of vertices)
  size <- dim(d)[1]
  
  #creating edge list from instance
  edges <- edgeList(d)
  
  edgeset <- c(TRUE, rep(FALSE, nrow(edges)-1))
  
  k <- 2
  
  while(sum(edgeset) < size){
    test.edgeset <- edgeset
    test.edgeset[k] <- TRUE
    if(compatibleEdge(edges, test.edgeset, k))
      edgeset <- test.edgeset  
    k <- k+1
  }
  
  obj <- sum(edgeset*edges$distance)
  
  sol <- edges[which(edgeset), ]
  
  return(list(obj=obj, sol=sol))
    
}

#----branch and bound for the TSP based on distances heuristic----

BandBTSPdist <- function(d, report.evol=FALSE){
  
  #instance size (number of vertices)
  size <- dim(d)[1]
  
  #creating edge list from instance
  edges <- edgeList(d)
  
  #setting the upper bound
  best.obj <- Inf
  best.sol <- logical(nrow(edges))
  
  #initializing node list with node 0
  
  node0 <- nodegen(edges, rep(FALSE, nrow(edges)), 0)
  nodes <- list(node0$node)
  
  if(report.evol){
    iter <- 0
    current.solution <- numeric()
    current.nodes <- numeric()
  }
  
  while(length(nodes) > 0){
    
    #selecting the node to branch
    bounds <- sapply(nodes, function(x) x$bound)
    conss <- sapply(nodes, function(x) x$cons)
    
    select <- rep(0, length(nodes))
    select[which(bounds==min(bounds))] <- bounds[which(bounds==min(bounds))]
    select <- select*conss
    
    node.to.branch <- which(select==max(select))[1]
    
    cons.to.branch <- nodes[[node.to.branch]]$cons + 1
    
    edgesetT <- nodes[[node.to.branch]]$edgeset
    edgesetT[cons.to.branch] <- TRUE
    
    edgesetF <- nodes[[node.to.branch]]$edgeset
    
    node.T <- nodegen(edges, edgesetT, cons.to.branch)
    node.F <- nodegen(edges, edgesetF, cons.to.branch)
    
    if(node.T$feasible & sum(edgesetT)==size){
      if(node.T$node$bound < best.obj){
        best.obj <- node.T$node$bound
        best.sol <- node.T$node$edgeset
      }
    }
    
    if(node.T$feasible & sum(edgesetT) < size)
      nodes[[length(nodes) + 1]] <- node.T$node
    
    if(node.F$feasible)
      nodes[[length(nodes) + 1]] <- node.F$node
    
    #pruning the nodes worse than current solution
    
    bounds <- sapply(nodes, function(x) x$bound)
    
    pruned <- rep(FALSE, length(nodes))
    pruned[node.to.branch] <- TRUE
    pruned[which(bounds >= best.obj)] <- TRUE
    
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
      current.solution[iter] <- best.obj
      current.nodes[iter] <- length(nodes)
    }
    
  }
  
  best.sol.edges <- edges[which(best.sol), ]
  
  if(report.evol)
    return(list(solution=list(sol=best.sol.edges, obj=best.obj), report=list(objs=current.solution, nodes=current.nodes)))
  else
    return(list(sol=best.sol.edges, obj=best.obj))
}  


#---- MST obtained with Prim algorithm ----

MSTPrim <- function(InputMatrix){
  
  #initialize size and get list of edges
  n <- dim(InputMatrix)[1]
  edges <- edgeList(InputMatrix)
  
  #if input matrix not symmetric Prim does not apply
  if(!isSymmetric(InputMatrix))
    return(success=FALSE)
  
  edges <- edges[which(edges$ori < edges$des), ]
  
  CoveredNodes <- rep(FALSE, n)
  CoveredNodes[1] <- TRUE
  
  sol <- numeric(0)
  
  while(sum(CoveredNodes) < n){
    k <- 1
    while(CoveredNodes[edges$ori[k]] + CoveredNodes[edges$des[k]] !=1)
      k <- k+1
    if(CoveredNodes[edges$ori[k]])
      CoveredNodes[edges$des[k]] <- TRUE
    else
      CoveredNodes[edges$ori[k]] <- TRUE
    
    sol <- c(sol, edges[k, ])
  }
  
  sol <- as.data.frame(matrix(unlist(sol), (n-1), 3, byrow = TRUE))
  colnames(sol) <- colnames(edges)
  return(list(success=TRUE, sol=sol))
}

#---- OneTree bound for the Symmetric TSP ---

OneTree <- function(InputMatrix, origin){
  
  #initialize size and edgelist
  n <- dim(InputMatrix)[1]
  edges <- edgeList(InputMatrix)
  
  #if input matrix not symmetric or origin too large 1-Tree does not apply
  if(!isSymmetric(InputMatrix) | origin > n)
    return(success=FALSE)
  
  #reduced matrix: supress row and column origin
  ReducedMatrix <- InputMatrix[-origin, -origin]
  
  #MST of the reduced matrix
  MST.RM <- MSTPrim(ReducedMatrix)$sol
  
  #MST has to be modified relabelling nodes larger than origin
  for(i in 1:(n-2)){
    if(MST.RM$ori[i] >= origin) MST.RM$ori[i] <- MST.RM$ori[i] + 1
    if(MST.RM$des[i] >= origin) MST.RM$des[i] <- MST.RM$des[i] + 1
  }
  
  #nodes from origin with edges of small distance
  NFO <- order(InputMatrix[origin, ])[2:3]
  
  #rows of data frame consisting of edges of smallest distance from origin
  edge1 <- sort(c(origin, NFO[1]))
  edge2 <- sort(c(origin, NFO[2]))
  edge1 <- c(InputMatrix[edge1[1], edge1[2]], edge1)
  edge2 <- c(InputMatrix[edge2[1], edge2[2]], edge2)
  
  #packing solution
  sol <- rbind(MST.RM, edge1, edge2)
  
  return(list(success=TRUE, sol=sol))
}










