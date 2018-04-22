setwd("~/Dropbox (UPC)/quanti_docs/FlowShop")

#--- reading a set of Taillard flowshop instances ----

read.Taillard <- function(case){
    
    url <- paste0("http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/flowshop.dir/tai", case, ".txt")
    
    text <- readLines(url)
    
    text.split <- strsplit(text, " ")
    clean <- lapply(text.split, function(x) x[which(nchar(x)!=0)])
    
    lines <- length(clean)
    k <- 1
    num.instance <- 1
    instances <- list()
    
    while(k < lines){
      
      refs <- as.numeric(clean[[k+1]])
      n <- refs[1]
      m <- refs[2]
      seed <- refs[3]
      upper <- refs[4]
      lower <- refs[5]
      
      tij <- numeric(0)
      
      for(i in (k+3):(k+2+m)) tij <- c(tij, as.numeric(clean[[i]]))
      
      tij <- matrix(tij, m, n, byrow=TRUE)
      
      instances[[num.instance]] <- list(m=m, n=n, seed=seed, upper=upper, lower=lower, tij=tij)
      
      num.instance <- num.instance+1
      k <- k + m  + 3
    }
    
    return(instances)
}

#---- reading Taillard flowshop instances ----

tai20.5 <- read.Taillard("20_5")
tai20.10 <- read.Taillard("20_10")
tai20.20 <- read.Taillard("20_20")

tai50.5 <- read.Taillard("50_5")
tai50.10 <- read.Taillard("50_10")
tai50.20 <-  read.Taillard("50_20")

tai100.5 <- read.Taillard("100_5")
tai100.10 <- read.Taillard("100_10")
tai100.20 <-  read.Taillard("100_20")

tai200.10 <- read.Taillard("200_10")
tai200.20 <-  read.Taillard("200_20")

tai500.20 <-  read.Taillard("500_20")

save.image("instances/TaillardFS.RData")
