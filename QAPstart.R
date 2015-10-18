#-----reading the Elshafei:77 file as an example---- 

file <- url("http://anjos.mgi.polymtl.ca/qaplib/data.d/els19.dat")

data <- scan(file)

close(file)

n <- data[1]

flow <- matrix(data[2:(n*n + 1)], n, n, byrow=TRUE)

distance <- matrix(data[(n*n + 2):(2*n*n + 1)], n , n)

optimum <- c(9,10,7,18,14,19,13,17,6,11,4,5,12,8,15,16,1,2,3)
value_optimum <- 17212548

#-------defining objective function------- 


perMatrix <- function(vec){
  n <- length(vec)
  m <- matrix(0, n , n)
  for(i in 1:n){
    m[i, vec[i]] <- 1
  }
  return(m)
}

objective.QAP <- function(flow, distance, permutation){
  pmatrix <- perMatrix(permutation)
  dist_perm <- pmatrix %*% distance %*% t(pmatrix)
  return(sum(flow * dist_perm))
}

#---------testing objective function------

obj.value <- objective.QAP(flow, distance, optimum)

obj.value == value_optimum
