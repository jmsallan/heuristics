load("experimentTSP.RData")
source("functionsTSP.R")

#----generating instances----

sizes <- rep(c(10, 20, 50, 100), each=5)
set.seed(3333)
seeds <- sample(1:1000, 20)
inst <- mapply(function(x, y) SampleTSP(x, y), sizes, seeds)

#----tabu search----

iter.ts <- 50
res.ts <- mapply(function(x, y, z) TS(1:x, y, z), sizes, inst, iter.ts)
z.ts <- unlist(res.ts[(1:20)*2])

#----simulated annealing----

iter.sa <- sizes*(sizes-1)*iter.ts/2
mu <- rep(c(100, 500, 1000), each=length(sizes))
runs.sa <- 5

set.seed(1313)
res.sa <- rep(mapply(function(x, y, z, t) SA(1:x, y, z, t), sizes, inst, iter.sa, mu), runs.sa)
z.sa <- unlist(res.sa[(1:(length(mu)*runs.sa))*2])

data.sa <- data.frame(z.sa, mu, rep(sizes, length(mu)/length(sizes)), rep(1:20, length(mu)/length(sizes)))
names(data.sa) <- c("z", "mu", "size", "inst")

#----data display----
#intro to dplyr package: https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html

library(dplyr)

mean.sa <- data.sa %>% group_by(inst, mu) %>% summarise(mean=mean(z), max=max(z), min=min(z), sd=sd(z))
names(mean.sa) <- c("inst", "mu", "mean", "max", "min", "sd")

#intro to reshape2: http://seananderson.ca/2013/10/19/reshape.html

library(reshape2)

mean.table <- dcast(mean.sa, inst ~ mu, value.var = c("mean"))
mean.table <- cbind(sizes, mean.table, z.ts)
names(mean.table)  <- c("size", "inst", "mu=100", "mu=500", "mu=1000", "tabu")

print(mean.table, digits=5)

save.image("experimentTSP.RData")
