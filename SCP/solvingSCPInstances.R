#set wd

source("code/functionsSCP.R")

#---- reading and solving SCP instances ----

url <- "http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/scp42.txt"

urls41_410 <- lapply(c("scp41", "scp42", "scp43", "scp44", "scp45", "scp46", "scp47", "scp48", "scp49", "scp410"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instances41_410 <- lapply(urls41_410, readCSP)

SolutionInstances41_410 <- lapply(Instances41_410, solveLP_SCP)

urls51_510 <- lapply(c("scp51", "scp52", "scp53", "scp54", "scp55", "scp56", "scp57", "scp58", "scp59", "scp510"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instances51_510 <- lapply(urls51_510, readCSP)

SolutionInstances51_510 <- lapply(Instances51_510, solveLP_SCP)

urls61_65 <- lapply(c("scp61", "scp62", "scp63", "scp64", "scp65"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instances61_65 <- lapply(urls61_65, readCSP)

SolutionInstances61_65 <- lapply(Instances61_65, solveLP_SCP)

urlsa1_a5 <- lapply(c("scpa1", "scpa2", "scpa3", "scpa4", "scpa5"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instancesa1_a5 <- lapply(urlsa1_a5, readCSP)

SolutionInstancesa1_a5 <- lapply(Instancesa1_a5, solveLP_SCP)

urlsb1_b5 <- lapply(c("scpb1", "scpb2", "scpb3", "scpb4", "scpb5"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instancesb1_b5 <- lapply(urlsb1_b5, readCSP)

SolutionInstancesb1_b5 <- lapply(Instancesb1_b5, solveLP_SCP)

urlsc1_c5 <- lapply(c("scpc1", "scpc2", "scpc3", "scpc4", "scpc5"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instancesc1_c5 <- lapply(urlsc1_c5, readCSP)

SolutionInstancesc1_c5 <- lapply(Instancesc1_c5, solveLP_SCP)

urlsd1_d5 <- lapply(c("scpd1", "scpd2", "scpd3", "scpd4", "scpd5"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instancesd1_d5 <- lapply(urlsd1_d5, readCSP)

SolutionInstancesd1_d5 <- lapply(Instancesd1_d5, solveLP_SCP)

urlse1_e5 <- lapply(c("scpe1", "scpe2", "scpe3", "scpe4", "scpe5"), function(x)  paste0("http://people.brunel.ac.uk/~mastjjb/jeb/orlib/files/", x, ".txt"))

Instancese1_e5 <- lapply(urlse1_e5, readCSP)

SolutionInstancese1_e5 <- lapply(Instancese1_e5, solveLP_SCP)


save.image("results/SCPresults.Rdata")
