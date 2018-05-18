#set up working directory

source("code/functionsKP.R")


#---- testing SAKP ----

#generating large instance
set.seed(111)
LargeInstance <- EasyInstance(1000, 0.4)

#starting solution
inisolLIEmpty <- sample(c(TRUE, FALSE), 1000, replace = TRUE, prob = c(0.1, 0.9))


#solving with heuristic
LIheuristic <- HeuristicKP(LargeInstance)

#optimal solution with LP
LILP <- LPKP(LargeInstance)


#simualated annealing behaviour

Tmax <- 1000

SA01 <- SAKP(Inst=LargeInstance, Tmax=Tmax, inisol = inisolLIEmpty, eval=TRUE)

plot(1:100, SA01$evalfit[1:100], type="l", col="red", xlab = "", ylab = "")
lines(1:100, SA01$evalbest[1:100], col="blue")
lines(1:100, SA01$evaltest[1:100], col="green")

#simulated annealing with same evals of objective function as TS

Tmax <- 500000
SA02 <- SAKP(Inst=LargeInstance, Tmax=500000, mu=1, inisol = inisolLIEmpty, eval=FALSE)


#tabu search behaviour

#aspiration condition
TS01 <- TSKP(Inst=LargeInstance, inisol = inisolLIEmpty, tabu.size = 100, iter=500, asp = TRUE, eval=TRUE)

plot(1:500, TS01$evalfit[1:500], type="l", col="red", xlab = "", ylab = "")
lines(1:500, TS01$evalbest[1:500], col="blue")

#no aspiration condition
TS02 <- TSKP(Inst=LargeInstance, inisol = inisolLIEmpty, tabu.size = 100, asp=FALSE, iter=500, eval=TRUE)

plot(1:500, TS02$evalfit, type="l", col="red", xlab = "", ylab = "")
lines(1:500, TS02$evalbest, col="blue")

#analyzing last iterations
plot(350:500, TS01$evalfit[350:500], type="l", col="red", xlab = "", ylab = "")
lines(350:500, TS01$evalbest[350:500], col="blue")
plot(350:500, TS02$evalfit[350:500], type="l", col="red", xlab = "", ylab = "")
lines(350:500, TS02$evalbest[350:500], col="blue")

#behaviour of last iterations with aspiration condition
TS01$move[450:500]
TS01$evalfit[450:500]
LargeInstance$u[c(781, 777, 380, 485, 721, 768, 679)]


TS03 <- TSKP(Inst=LargeInstance, inisol = LIheuristic$sol, tabu.size = 100, asp=FALSE, iter=500, eval=TRUE)
SA03 <- SAKP(Inst=LargeInstance, Tmax=500000, mu=1, inisol = LIheuristic$sol, eval=FALSE)
SA04 <- SAKP(Inst=LargeInstance, Tmax=500000, mu=1, inisol = inisolLIEmpty, eval=FALSE)

plot(1:500, TS03$evalfit, type="l", col="red", xlab = "", ylab = "")
lines(1:500, TS03$evalbest, col="blue")

plot(1:10, TS03$evalfit[1:10], type="l", col="red", xlab = "", ylab = "")
lines(1:10, TS03$evalbest[1:10], col="blue")
