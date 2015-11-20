#set the working directory (has to be changed for each computer)

setwd("~/Dropbox/00_TfDA1415q1/class04")

#loeading functions

source("tsp.R")

#testing NearestNeighbour

d1 <- matrix(c(0, 3, 4, 2, 7, 3, 0, 4, 6, 3, 4, 4, 0, 5, 8, 2, 6, 5, 0, 6, 7, 3, 8, 6, 0), 5, 5)

d2 <- matrix(c(0, 4, 6, 7, 4, 0, 8, 10, 6, 8, 0, 12, 7, 10, 12, 0), 4, 4)

NearestNeighbour(d1, 1)

for(i in 1:4)
    print(NearestNeighbour(d2, i))

#testing SwapNeighbourNodes

a <- 1:10

for(i in 1:10)
    print(SwapNeighbourNodes(a, i))


DistanceMatrixAtt <- function(data){
    n <- nrow(data)
    d <- matrix(numeric(), n, n)
    
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            xd <- data$xcoord[i] - data$xcoord[j]
            yd <- data$ycoord[i] - data$ycoord[j]
            
            rij <- sqrt((xd*xd + yd*yd) /10)
            
            tij <- round(rij) #caveat: rounding to nearest event number
            
            if(tij < rij){
                d[i, j] <- tij + 1}
            else{ 
                d[i, j] <- tij}
            
            d[j, i] <- d[i, j]
        }
    }
    
    for(i in 1:n)
        d[i, i] <- Inf
    
    return(d)
}

#reading ATT data48

data48 <- read.table(file="att48data.tsp", header=FALSE)
colnames(data48) <- c("node", "xcoord", "ycoord")

distance.att48 <- DistanceMatrixAtt(data48) 

#setting seed for (pseudo) random numbers

set.seed(1)

#solution by nearest neighbour

nn.att48 <- NearestNeighbour(distance.att48, 1)

#solution by simmulated annealing

sa.att48 <- SimulatedAnnealingTSP(distance.att48, 1000, initial.random=FALSE)
sa.att48b <- SimulatedAnnealingTSP(distance.att48, 1000, initial.random=TRUE)

#solution by genetic algorithm

ga.att48 <- GeneticAlgorithmTSP(distance.att48, 100, 200, 0.2, seed = FALSE)

time1 <- Sys.time()
ga.att48 <- GeneticAlgorithmTSP(distance.att48, 100, 100, 0.2, seed = FALSE)
time2 <- Sys.time()

time2 - time1
#solution by tabu search

time1 <- Sys.time()
ts.att48 <- TabuSearchTSP(distance.att48, 1000, initial.random = FALSE)
time2 <- Sys.time()

time2 - time1
#reading a280 data

data280 <- read.table(file="a280data.tsp", header=FALSE)
colnames(data280) <- c("node", "xcoord", "ycoord")

distance.att280 <- DistanceMatrixAtt(data280)

#optimal value of a280
opta280 <- read.table(file="a280data.opt.tour", header=FALSE)
obj.opt280 <- print(DistanceTSP(distance.att280, opta280[[1]]))

#different strategies for a280

nn.280 <- NearestNeighbour(distance.att280, 1)

sa.280 <- SimulatedAnnealingTSP(distance.att280, 1000, initial.random = FALSE)

ga.280 <- GeneticAlgorithmTSP(distance.att280, 20, 10, 0.2, seed = FALSE)

ga.280b <- GeneticAlgorithmTSP(distance.att280, 20, 10, 0.2, seed = TRUE)

ts.280 <- TabuSearchTSP(distance.att280, 100, initial.random = TRUE)

#assessing time of execution

time1 <- Sys.time()

ts.280b <- TabuSearchTSP(distance.att280, 100, initial.random = FALSE)

time2 <- Sys.time()

time2 - time1

time3 <- Sys.time()

nn.280 <- NearestNeighbour(distance.att280, 1)

time4 <- Sys.time()

time4 - time3