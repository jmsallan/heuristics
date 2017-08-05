# heuristics

This is the gitHub repository of heuristics of Jose M Sallan. Here you will find quite simple heuristics to attack combinatorial problems, as the code here is set for teaching purposes.

At this moment this repository contains the following material

##Travelling Salesperson Problem (TSP)

* **functionsTSP.R** contains several functions and algorithms to tackle the TSP using R, including an implementation of the nearest neighbor heuristic, and implementations of hill climbing, tabu search, simulated annealing and genetic algorithms. The file also contains a function to generate random instances of a given size.
* **demoTSP.R** contains examples of use of the functions of the previous file, one on a randomly generated instance of size 20, and another on the **att48** instance.
* **experimentTSP.R** contains a experiment to compare the performance of SA and TS heuristics (SA wins!). Several interesting utilities of R are used here, e.g. the *apply* family functions and functions from *dplyr* and *reshape2* packages. For convenience, results of the code are stored in **experimentTSP.RData**.
* **testingGA.R** tests the application of the genetic algorithm to the instance **att48**. The results of running the code are stored in **testGA.RDAata**.


## Quadratic assignment problem (QAP)

* In **QAPfunctions.R** are listed several functions and algorithms to tackle the QAP using R, based on the ones defined for the TSP.
* **QAPresults** solves several instances of QAP using the heuristics defined in **QAPfunctions.R**. Illustrates how to get instances from the Internet, and how to use markdown to write reproducible reports.
* **QAPtypewriter.R** evaluates several large instances of QAP as in the previous file. Here I have used the stargazer package to build text and LaTeX outputs.

##Knapsack Problem (KP)

I have developed a genetic algorithm to solve the Knapsack problem. It can be easily tackled through linear programming, but this GA is illustrative about how to deal with a maximum problem with unfeasible solutions. The heuristic and the instances are packed in a single script.

In the **functionsKP.R** can be found several instance generators of the KP and of the subset sum (SS) problem. The SS is a variant of the KP where all items have an utility equal to its weight. There are two heuristics implemented: a genetic algorithm (which works somewhat well for small instances) and a branch and bound, which returns the optimal solution.


