# 6/21/19  
First attempts at creating an R package as part of a project workflow. Also trying this is as a better way to keep track of progress. This document is a lab notebook for this project, idea being that it can contain stream-of-consiousness thoughts/musings/ideas as this project develops outside of the more formal analysis folder and the other project infrastructure. 

Main goal for today was to transfer key functions from the [elimination feasibility project](https://github.com/cmhoove14/EliminationFeasibility) into this project in an attempt to better organize functions from that project that will be reused here. These functions include the dynamic model itself, the function to estimate $R_eff$, and helper functions used to estimate density dependencies

# 6/26/19  
Continued transferring over functions and simulations from previous work and started working on a document (`Analysis/Model_Sims.Rmd`) to demonstrate functionality and keep track of the functions implemented in the package and what they do. Now able to produce generic $R_eff$ curve and simulate the model through time. Next want to use `gganimate` to show how $R_eff$ estimates change over the course of different intervention campaigns, e.g. annual MDA, MDA+snail control, MDA+prawn intervention, etc.

# 6/27/19  
Transferred over code for the stochastic model implemented in `adaptivetau`

# 7/1/2019  
Added deterministic and stochastic versions of an age stratified model with non-linear man-to-snail force of infection as described in [Gurarie et al 2018](https://doi.org/10.1371/journal.pntd.0006514) and functions to simulate it. Added functionality to simulate models with or without events. Modified base model to include clumping parameter responsive to mean worm burden. Incorporated clumping parameter function derived from Senegal data relationship between log of mean worm burden and clumping parameter

# To-dos  
### Incorporate dynamic clumping parameter into $R_{eff}$ formulation
### Derive $R_{eff}$ for age stratified model  
### Document deterministic models (equations, base parameters) in vignettes
#### Document derivation of $R_{eff}$  
### Fit deterministic versions of models in .Rmds  
