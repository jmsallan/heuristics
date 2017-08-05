setwd("~/Dropbox (UPC)/quanti_docs/KP")

source("functionsKP.R")

easy <- list(w=c(30, 45, 15, 40, 35, 5), u=c(24, 108, 30, 100, 80.5, 3), W=100)
hard <- list(w=c(30, 45, 15, 70, 35, 5), u=c(24, 108, 30, 175, 80.5, 3), W=100)

InstanceKP10 <- KPInstance01(1000, 10, 40)
InstanceKP50 <- KPInstance01(1000, 50, 40, 2424)
InstanceKP100  <- KPInstance01(1000, 100, 40, 3333)

bandb.easy <- BandBKP(easy)
bandb.hard <- BandBKP(hard)

bandb.KP10 <- BandBKP(InstanceKP10)
bandb.KP50 <- BandBKP(InstanceKP50)
bandb.KP100 <- BandBKP(InstanceKP100)

bandb.KP50.report <- BandBKP(InstanceKP50, TRUE)
bandb.KP100.report <- BandBKP(InstanceKP100, TRUE)

LP.KP10 <- LPKP(InstanceKP10)
LP.KP50 <- LPKP(InstanceKP50)
LP.KP100 <- LPKP(InstanceKP100)


identical(bandb.KP100$sol, LP.KP100$sol)
identical(bandb.KP50$sol, LP.KP50$sol)
identical(bandb.KP10$sol, LP.KP10$sol)
