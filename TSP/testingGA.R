load("testGA.RData")

source("code/functionsTSP.R")

#----reading att48-----

data48 <- read.table(file="instances/att48data.tsp", header=FALSE)
colnames(data48) <- c("node", "xcoord", "ycoord")

distance.att48 <- DistanceMatrixAtt(data48)

#----testing number of iterations----

iter <- 1000
ga.test <- GA(distance.att48, 50, iter, 0.2, report.evol=TRUE)

pdf("results/plotGA.pdf", width=7, height=5)
plot(1:iter, ga.test$evol, type="n", xlab="iter", ylab="dist", main="assessing iterations for the GA")
lines(1:iter, ga.test$evol, col="blue")
dev.off()

#parameters to play with:
#number of iterations: 2500, 1000, 500
#npop: 20, 50, 100
#pmut: 0, 0.2, 0.5, 0.8, 1
#number of iterations should be linked with population size, same eval of objective functions: 100.000
runs <- 10
muts <- c(0, 0.2, 0.5, 0.8, 1)
npops <- c(20, 50, 100)
iters <- c(2500, 1000, 500)
nmuts <- length(muts)
nnpops <- length(npops)

npop <- rep(rep(npops, each=nmuts), runs)
iter <- rep(rep(iters, each=nmuts), runs)
pmut <- rep(rep(muts, nnpops), runs)

set.seed(1212)
results <- mapply(function(x, y, z) GA(distance.att48, x, y, z, verbose=TRUE), npop, iter, pmut)

z <- unlist(results[(1:length(npop)*2)])

z.table <- data.frame(npop, iter, pmut, z)

#----linear models of objective function on parameters----

mod01 <- lm(z ~ factor(npop) + factor(pmut), data=z.table)
mod02 <- lm(z ~ factor(npop)*factor(pmut), data=z.table)

summary(mod01)
summary(mod02)
anova(mod01, mod02)

#----displaying results----

library(dplyr)

z.results <- z.table %>% group_by(npop, pmut) %>% summarise(mean=mean(z), max=max(z), min=min(z), sd=sd(z))

z.results

save.image("results/testGA.Rdata")


