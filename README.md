# Matlab-Scripts

This repository contains scripts used in the project The Evolution of Sex Ratios and Sex Allocation: a
Theoretical Study, submitted on April 2021 by candidate 1048690 of the Oxford Masters in Mathematical Sciences.

The main scripts are:
 - RunOptimal.m : finds optimal sex and investment allocation solving Optimizaiton Program 1 (Basic) and generates figure 2.1. (Section 2.2 of the maint text).
 - RunOptimalUnc.m: finds optimal investment allocation solving Optimizaiton Program 5 (with Uncertainty) and generates figures 2.2 and
   2.3 (Section 2.3 of the maint text).
 - RunOptimalInv.m: finds optimal investment allocation solving Optimizaiton Program 6 (Investment Only) and generates figure 2.4 (Section 2.4 of the main text).
 - MultipleGens.m: simulates the evolution of a population in which mothers determine sex allocation according to Optimizaiton Program 5 (with Uncertainty). (Section 3.2.3 of the maint text).
 - MultipleGensInv.m: simulates the evolution of a population in which mothers determine investment allocation according to Optimizaiton Program 6 (Investment Only). (Section 3.2.3 of the maint text).

The folders supp_functions.m and supp_data.m contain the supplementary scripts and data used in the main scripts.

To change the parameters and/or functional forms of the Fixed-Brood model, modify the scripts co.m (offspring condition),
invs.m (investemnt capability) or f_male.m and f_female.m (male and female fitness). These are all in the supp_functions folder.

To change parameters such as population size, number of iterations, etc., see the main scripts.

All scripts are commented to be user friendly but we refer to the main text of the project for definitions of parameters, functions, etc. 

All scripts were developed by the author except for limitdist.m, used to find the limiting distribution from the 
probability transition matrix of a markov chain. This was taken from http://www.math.wustl.edu/~feres/Math450Lect04.pdf .

