#setwd!

source("code/functionsTSP2.R")


#---- all 2opt moves and testgain function ----

n <- 50
D <- SampleTSP(n)
sol <- sample(1:n,n)
M <- matrix(numeric(n*(n-3)*4/2), n*(n-3)/2, 4)

T <- 1
for(i in 1:(n-2)){
  for(k in (i+1):min(n+i-3,n-1)){
    
    print(c(i,k))
    print(move2opt(sol, i, k))
    if(i==1)
      testgain <- D[sol[n], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[n] ,sol[i]] - D[sol[k], sol[k+1]]
    else
      testgain <- D[sol[i-1], sol[k]] +  D[sol[i], sol[k+1]] - D[sol[i-1] ,sol[i]] - D[sol[k], sol[k+1]]
    print(testgain)
    
    M[T, 1] <- TSP(D, sol)$d
    M[T, 2] <- testgain
    M[T, 3] <- TSP(D, move2opt(sol, i, k))$d
    M[T, 4] <- M[T, 1] + M[T, 2]

    T <- T+1
    
}
}
#---- a sample instance of size 10 ----

Instance10 <- SampleTSP(10)

sol110 <- TSP(Instance10, 1:10)
solnn <- NearestNeighbour(Instance10)


#---- testing Hill Climbing ----

TestSample <- SampleTSP(30, seed=1313)

sol <- 1:30
sol2 <- NearestNeighbour(TestSample)
TSP(TestSample, sol)

HC01 <- HillClimbing2opt(TestSample, sol)
HC02 <- HillClimbing2opt(TestSample, sol2$sol)

#----testing tabu search----

TS01 <- TSTSP2opt(TestSample, sol, asp=TRUE, eval=TRUE)

pdf("tabu.pdf", height=5, width = 10)
par(mfrow=c(1,2))
plot(1:100, TS01$evalfit, type = "l", col = "red", xlab="", ylab="")
lines(1:100, TS01$evalbest, col="blue")
plot(60:100, TS01$evalfit[60:100], type = "l", col = "red", xlab="", ylab="")
lines(60:100, TS01$evalbest[60:100], col="blue")
dev.off()

plot(1:100, TS01$evalgain, type = "l", col = "blue", lwd=2, xlab="", ylab="")

#----search without tabu list----

Silly01 <- SillyTSP2opt(TestSample, sol)

pdf("silly.pdf", height=5, width = 10)
par(mfrow=c(1,2))
plot(1:100, Silly01$evalfit, type = "l", col = "red", xlab="", ylab="")
lines(1:100, Silly01$evalbest, col="blue")
plot(60:100, Silly01$evalfit[60:100], type = "l", col = "red", xlab="", ylab="")
lines(60:100, Silly01$evalbest[60:100], col="blue")
dev.off()

plot(1:100, Silly01$evalgain, type = "l", col = "blue", lwd=2, xlab="", ylab="", main="gain")

#----testing simulated annealing----

set.seed(1313)
SA01 <- SATSP2opt(TestSample, sol, Tmax=3000, mu=10, eval = TRUE)

pdf("SAevol.pdf", height=5, width = 10)
par(mfrow=c(1,2))
plot(1:3000, SA01$evaltest, type = "l", col = "green", xlab="", ylab="")
lines(1:3000, SA01$evalbest, col="blue")
lines(1:3000, SA01$evalfit, col ="red")

plot(2500:3000, SA01$evaltest[2500:3000], type = "l", col = "green", xlab="", ylab="")
lines(2500:3000, SA01$evalbest[2500:3000], col="blue")
lines(2500:3000, SA01$evalfit[2500:3000], col ="red")
dev.off()

set.seed(1313)
SA02 <- SATSP2opt(TestSample, sol, Tmax=40500, mu=10, eval = FALSE)


#----testing GRASP----

for(i in 1:10) print(GRASP(TestSample, 2))

set.seed(1111)
GRASP01 <- GRASP2opt(TestSample, rcl.size = 2, tries = 10)

set.seed(1111)
GRASP02 <- GRASP2opt(TestSample, rcl.size = 2, ls="SA", tries = 10)

GRASP01$report
GRASP02$report
