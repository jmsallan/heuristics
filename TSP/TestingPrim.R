source("code/functionsTSPedges.R")
source("code/functionsTSP.R")

load("results/onetree.RData")

case1 <- matrix(c(0, 5, 12, 18, 12, 
                  5, 0, 5, 112, 6, 
                  12, 5, 0, 18, 12, 
                  18, 112, 18, 0, 11, 
                  12, 6, 12, 11, 0), 5, 5)

M <- sum(case1)*12

#all solutions
matrix0 <- case1
tree0 <- OneTree(matrix0, 2)$sol
bound0 <- sum(tree0$distance)

#branching 0 into 1, 2, 3

matrix1 <- matrix0
matrix1[1,5] <- M
matrix1[5,1] <- M
tree1 <- OneTree(matrix1, 2)$sol
bound1 <- sum(tree1$distance)

matrix2 <- matrix0
matrix2[3,5] <- M
matrix2[5,3] <- M
tree2 <- OneTree(matrix2, 2)$sol
bound2 <- sum(tree2$distance)

matrix3 <- matrix0
matrix3[4,5] <- M
matrix3[5,4] <- M
tree3 <- OneTree(matrix3, 2)$sol
bound3 <- sum(tree3$distance)

#branching 1 into 4, 5, 6

matrix4 <- matrix1
matrix4[1,3] <- M
matrix4[3,1] <- M
tree4 <- OneTree(matrix4, 2)$sol
bound4 <- sum(tree4$distance)

matrix5 <- matrix1
matrix5[2,3] <- M
matrix5[3,2] <- M
tree5 <- OneTree(matrix5, 2)$sol
bound5 <- sum(tree5$distance)

matrix6 <- matrix1
matrix6[5,3] <- M
matrix6[3,5] <- M
tree6 <- OneTree(matrix6, 2)$sol
bound6 <- sum(tree6$distance)

#branching 2 into 7, 8 , 9

matrix7 <- matrix2
matrix7[1,5] <- M
matrix7[5,1] <- M
tree7 <- OneTree(matrix7, 2)$sol
bound7 <- sum(tree7$distance)

matrix8 <- matrix2
matrix8[1,3] <- M
matrix8[3,1] <- M
tree8 <- OneTree(matrix8, 2)$sol
bound8 <- sum(tree8$distance)

matrix9 <- matrix2
matrix9[1,2] <- M
matrix9[2,1] <- M
tree9 <- OneTree(matrix9, 2)$sol
bound9 <- sum(tree9$distance)

#branching 5 into 10, 11, 12

matrix10 <- matrix5
matrix10[2,5] <- M
matrix10[5,2] <- M
tree10 <- OneTree(matrix10, 2)$sol
bound10 <- sum(tree10$distance)

matrix11 <- matrix5
matrix11[3,5] <- M
matrix11[5,3] <- M
tree11 <- OneTree(matrix11, 2)$sol
bound11 <- sum(tree11$distance)

matrix12 <- matrix5
matrix12[4,5] <- M
matrix12[5,4] <- M
tree12 <- OneTree(matrix12, 2)$sol
bound12 <- sum(tree12$distance)

#branching 9 into 13, 14, 15

matrix13 <- matrix9
matrix13[1,5] <- M
matrix13[5,1] <- M
tree13 <- OneTree(matrix13, 2)$sol
bound13 <- sum(tree13$distance)

matrix14 <- matrix9
matrix14[2,5] <- M
matrix14[5,2] <- M
tree14 <- OneTree(matrix14, 2)$sol
bound14 <- sum(tree14$distance)

matrix15 <- matrix9
matrix15[4,5] <- M
matrix15[5,4] <- M
tree15 <- OneTree(matrix15, 2)$sol
bound15 <- sum(tree15$distance)

#branching 3 into 16, 17, 18

matrix16 <- matrix3
matrix16[1,2] <- M
matrix16[2,1] <- M
tree16 <- OneTree(matrix16, 2)$sol
bound16 <- sum(tree16$distance)

matrix17 <- matrix3
matrix17[1,4] <- M
matrix17[4,1] <- M
tree17 <- OneTree(matrix17, 2)$sol
bound17 <- sum(tree17$distance)

matrix18 <- matrix3
matrix18[1,5] <- M
matrix18[5,1] <- M
tree18 <- OneTree(matrix18, 2)$sol
bound18 <- sum(tree18$distance)

#branching 17 into 19, 20, 21

matrix19 <- matrix17
matrix19[2,3] <- M
matrix19[3,2] <- M
tree19 <- OneTree(matrix19, 2)$sol
bound19 <- sum(tree19$distance)

matrix20 <- matrix17
matrix20[3,4] <- M
matrix20[4,3] <- M
tree20 <- OneTree(matrix20, 2)$sol
bound20 <- sum(tree20$distance)

matrix21 <- matrix17
matrix21[3,5] <- M
matrix21[5,3] <- M
tree21 <- OneTree(matrix21, 2)$sol
bound21 <- sum(tree21$distance)

#branching 19 into 22, 23, 24

matrix22 <- matrix19
matrix22[1,5] <- M
matrix22[5,1] <- M
tree22 <- OneTree(matrix22, 2)$sol
bound22 <- sum(tree22$distance)

matrix23 <- matrix19
matrix23[3,5] <- M
matrix23[5,3] <- M
tree23 <- OneTree(matrix23, 2)$sol
bound23 <- sum(tree23$distance)

matrix24 <- matrix19
matrix24[2,5] <- M
matrix24[5,2] <- M
tree24 <- OneTree(matrix24, 2)$sol
bound24 <- sum(tree24$distance)



#---- reproducing B and B for the small instance ----

SmallInstance <- matrix(c(0, 3, 4, 2, 12, 
                          3, 0, 4, 6, 3, 
                          4, 4, 0, 5, 8, 
                          2, 6, 5, 0, 6, 
                          12, 3, 8, 6, 0), 5, 5)

test <- MSTPrim(SmallInstance)


M <- sum(SmallInstance)*12

#all solutions
matrix0 <- SmallInstance
tree0 <- OneTree(matrix0, 1)$sol
bound0 <- sum(tree0$distance)

#branch node 0 excluding:
#(2,5) -> node 1 
#(2,3) -> node 2
#(1,2) -> node 3

matrix1 <- matrix0
matrix1[2,5] <- M
matrix1[5,2] <- M
tree1 <- OneTree(matrix1, 1)$sol
bound1 <- sum(tree1$distance)

matrix2 <- matrix0
matrix2[2,3] <- M
matrix2[3,2] <- M
tree2 <- OneTree(matrix2, 1)$sol
bound2 <- sum(tree2$distance)

matrix3 <- matrix0
matrix3[1,2] <- M
matrix3[2,1] <- M
tree3 <- OneTree(matrix3, 1)$sol
bound3 <- sum(tree3$distance)

#branch node 3 excluding:
#(2,3)
#(3,4)
#(1,3)

matrix4 <- matrix3
matrix4[2,3] <- M
matrix4[3,2] <- M
tree4 <- OneTree(matrix4, 1)$sol
bound4 <- sum(tree4$distance)

matrix5 <- matrix3
matrix5[3,4] <- M
matrix5[4,3] <- M
tree5 <- OneTree(matrix5, 1)$sol
bound5 <- sum(tree5$distance)

matrix6 <- matrix3
matrix6[1,3] <- M
matrix6[3,1] <- M
tree6 <- OneTree(matrix6, 1)$sol
bound6 <- sum(tree6$distance)

#---- sample exercices with five nodes ----

TSPSample <- function(n, range, euc=TRUE){
  
  x <- sample(1:range, n, replace = TRUE)
  y <- sample(1:range, n, replace = TRUE)
  
  coord <- as.data.frame(matrix(c(x, y), n, 2, byrow = TRUE))
  
  Distance <- Distancia(coord, euc)
  
  return(Distance)
}

set.seed(1313)

Caso1 <- TSPSample(5, 20, euc=FALSE)
Caso2 <- TSPSample(5, 20, euc=FALSE)

library(xtable)

xtable(Caso1, caption = "Assigment1, case 1", label = "tab:case1", digits = 0, align="cccccc")

xtable(Caso2, caption = "Assigment1, case 2", label = "tab:case2", digits = 0, align="cccccc")

#---- solving Caso1 ----

#all solutions
Caso1matrix0 <- Caso1
Caso1tree0 <- OneTree(Caso1matrix0, 1)$sol
Caso1bound0 <- sum(Caso1tree0$distance)

Caso1Matrix1 <- Caso1matrix0
Caso1Matrix1[1,5] <- M
Caso1Matrix1[5,1] <- M
Caso1tree1 <- OneTree(Caso1Matrix1, 1)$sol
Caso1bound1 <- sum(Caso1tree1$distance)

Caso1Matrix2 <- Caso1matrix0
Caso1Matrix2[2,5] <- M
Caso1Matrix2[5,2] <- M
Caso1tree2 <- OneTree(Caso1Matrix2, 1)$sol
Caso1bound2 <- sum(Caso1tree2$distance)

Caso1Matrix3 <- Caso1matrix0
Caso1Matrix3[4,5] <- M
Caso1Matrix3[5,4] <- M
Caso1tree3 <- OneTree(Caso1Matrix3, 1)$sol
Caso1bound3 <- sum(Caso1tree3$distance)

Caso1Matrix4 <- Caso1Matrix2
Caso1Matrix4[1,5] <- M
Caso1Matrix4[5,1] <- M
Caso1tree4 <- OneTree(Caso1Matrix4, 1)$sol
Caso1bound4 <- sum(Caso1tree4$distance)

Caso1Matrix5 <- Caso1Matrix2
Caso1Matrix5[3,5] <- M
Caso1Matrix5[5,3] <- M
Caso1tree5 <- OneTree(Caso1Matrix5, 1)$sol
Caso1bound5 <- sum(Caso1tree5$distance)

Caso1Matrix6 <- Caso1Matrix2
Caso1Matrix6[4,5] <- M
Caso1Matrix6[5,4] <- M
Caso1tree6 <- OneTree(Caso1Matrix6, 1)$sol
Caso1bound6 <- sum(Caso1tree6$distance)

Caso1Matrix7 <- Caso1Matrix1
Caso1Matrix7[1,2] <- M
Caso1Matrix7[2,1] <- M
Caso1tree7 <- OneTree(Caso1Matrix7, 1)$sol
Caso1bound7 <- sum(Caso1tree7$distance)

Caso1Matrix8 <- Caso1Matrix1
Caso1Matrix8[2,3] <- M
Caso1Matrix8[3,2] <- M
Caso1tree8 <- OneTree(Caso1Matrix8, 1)$sol
Caso1bound8 <- sum(Caso1tree8$distance)

Caso1Matrix9 <- Caso1Matrix1
Caso1Matrix9[5,2] <- M
Caso1Matrix9[2,5] <- M
Caso1tree9 <- OneTree(Caso1Matrix9, 1)$sol
Caso1bound9 <- sum(Caso1tree9$distance)

Caso1Matrix10 <- Caso1Matrix9
Caso1Matrix10[1,3] <- M
Caso1Matrix10[3,1] <- M
Caso1tree10 <- OneTree(Caso1Matrix10, 1)$sol
Caso1bound10 <- sum(Caso1tree10$distance)

Caso1Matrix11 <- Caso1Matrix9
Caso1Matrix11[2,3] <- M
Caso1Matrix11[3,2] <- M
Caso1tree11 <- OneTree(Caso1Matrix11, 1)$sol
Caso1bound11 <- sum(Caso1tree11$distance)

Caso1Matrix12 <- Caso1Matrix9
Caso1Matrix12[3,5] <- M
Caso1Matrix12[5,3] <- M
Caso1tree12 <- OneTree(Caso1Matrix12, 1)$sol
Caso1bound12 <- sum(Caso1tree12$distance)

Caso1Matrix15 <- Caso1Matrix8
Caso1Matrix15[4,5] <- M
Caso1Matrix15[5,4] <- M
Caso1tree15 <- OneTree(Caso1Matrix15, 1)$sol
Caso1bound15 <- sum(Caso1tree15$distance)

library(gtools)
perm5 <- permutations(4, 4, 2:5)
sols <- apply(perm5, 1, function(x) TSP(Caso1, c(1,x))$d)
#---- solving Caso2 ----

#all solutions
Caso2matrix0 <- Caso2
Caso2tree0 <- OneTree(Caso2matrix0, 1)$sol
Caso2bound0 <- sum(Caso2tree0$distance)

#solution 12345, with total distance 58

save.image("results/onetree.RData")
