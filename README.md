# heuristics

This is the gitHub repository of heuristics of [Jose M Sallan](http://josemsallan.blogspot.com/). Here you will find quite simple heuristics to attack combinatorial problems, as the code here is set for teaching purposes.

At this moment this repository contains the following material:

## Travelling Salesperson Problem (TSP)

###code

The folder **code** contains two sets of functions:

* **functionsTSP.R** contains several functions and algorithms to tackle the TSP using R, including an implementation of the nearest neighbor heuristic, and implementations of hill climbing, tabu search, simulated annealing and genetic algorithms. All these implementations code instances of the TSP as a permutation starting from node 1. The file also contains a function to generate random instances of a given size.
* **functionsTSPedges.R** contains functions to solve the TSP based on a edge list representation. Contains a function implementing an heuristic based on picking compatible edges on a decreasing number of distances, and a branch and bound heuristic. The later works in a reasonable time with instances of around 10 edges.

###demos

In the root of the folder there are some scripts using the functions of **code**. The results (data files and graphs) of those implementations are stored in the folder **results**.

* **demoTSP.R** contains examples of use of the functions of the previous file, one on a randomly generated instance of size 20, and another on the **att48** instance.
* **demoTSPedges.R** contains examples of the implementation of the branch and bound heuristic implemented in **functionsTSPedges.R**. The results are stored in **resultsdemoTSPedges.RData**.
* **experimentTSP.R** contains a experiment to compare the performance of SA and TS heuristics (SA wins!). Several interesting utilities of R are used here, e.g. the *apply* family functions and functions from *dplyr* and *tidyr* packages. For convenience, results of the code are stored in **experimentTSP.RData**.
* **testingGA.R** tests the application of the genetic algorithm to the instance **att48**. The results of running the code are stored in **testGA.RDAata**.
* **TestingPrim.R** is a preliminary implementation of a lower bound of the TSP based on 1-tree.

## Quadratic assignment problem (QAP)

* In **QAPfunctions.R** are listed several functions and algorithms to tackle the QAP using R, based on the ones defined for the TSP.
* **QAPresults** solves several instances of QAP using the heuristics defined in **QAPfunctions.R**. Illustrates how to get instances from the Internet, and how to use markdown to write reproducible reports. Results are stored in **QAPresults.RData**.
* **QAPtypewriter.R** evaluates several large instances of QAP as in the previous file. Here I have used the stargazer package to build text and LaTeX outputs. Results can be found in **QAPresultsTypewriters.RData**.

## Knapsack Problem (KP)

* I have developed a genetic algorithm and a branch and bound method to solve the Knapsack problem (KP). KP can be easily tackled through linear programming, but I have used GA to illustrate about how to deal with a maximum problem with unfeasible solutions. The heuristic and the instances are packed in a single script. The script **functionsKP.R** contains the functions that execute both algorithms, together with instance generators for the KP and the subset sum(SS) problem. The SS is a variant of the KP where all items have an utility equal to its weight. The script **demoKP.R** demonstrates how do these functions work.

* The **2017-03-22gaKP** files are a RMarkdown document explaining how the GA works.

