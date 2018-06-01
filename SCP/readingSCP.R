
url <- "http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/scp42.txt"



# The format of all of these 80 data files is:
# number of rows (m), number of columns (n)
# the cost of each column c(j),j=1,...,n
# for each row i (i=1,...,m): the number of columns which cover
# row i followed by a list of the columns which cover row i

readCSP <- function(url){
    a <- readLines(url)
    b <- lapply(a, function(x) unlist(strsplit(x, " ")))
    c <- lapply(b, function(x) as.numeric(x[-1]))
    d <- sapply(c, function(x) length(x))
    m <- c[[1]][1]
    n <- c[[1]][2]

    counter <- 2
    
    elements <- 0
    nextcounter <- counter
    while(elements < n){
      elements <- elements + d[nextcounter]
      nextcounter <- nextcounter + 1
    }
    
    costs <- numeric(0)
    for(i in counter:(nextcounter-1)) costs <- c(costs, c[[i]])
    
    counter <- nextcounter
    
    M <- matrix(rep(0, m*n), m, n)
    
    for(i in 1:m){
    
      elements <- 0
      nextcounter <- nextcounter+1
      while(elements < c[[counter]]){    
        elements <- elements + d[nextcounter]
        nextcounter <- nextcounter + 1
      }
      
      for(j in (counter+1):(nextcounter-1)){
        M[i , c[[j]]] <- 1
      }
      counter <- nextcounter
    }

return(list(M=M, costs=costs))
}

url <- "http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/scp42.txt"

paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", "scp42", ".txt")

urls <- lapply(c("scp41", "scp42", "scp43", "scp44", "scp45", "scp46", "scp47", "scp48", "scp49", "scp410"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

allelements <- lapply(urls, readCSP)

save.image("SCPinstances.RData")

matricesRDS <- saveRDS(allelements, "SCPinstances.RDS")
