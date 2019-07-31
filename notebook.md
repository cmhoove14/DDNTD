# 6/21/19  
First attempts at creating an R package as part of a project workflow. Also trying this is as a better way to keep track of progress. This document is a lab notebook for this project, idea being that it can contain stream-of-consiousness thoughts/musings/ideas as this project develops outside of the more formal analysis folder and the other project infrastructure. 

Main goal for today was to transfer key functions from the [elimination feasibility project](https://github.com/cmhoove14/EliminationFeasibility) into this project in an attempt to better organize functions from that project that will be reused here. These functions include the dynamic model itself, the function to estimate $R_eff$, and helper functions used to estimate density dependencies

# 6/26/19  
Continued transferring over functions and simulations from previous work and started working on a document (`Analysis/Model_Sims.Rmd`) to demonstrate functionality and keep track of the functions implemented in the package and what they do. Now able to produce generic $R_eff$ curve and simulate the model through time. Next want to use `gganimate` to show how $R_eff$ estimates change over the course of different intervention campaigns, e.g. annual MDA, MDA+snail control, MDA+prawn intervention, etc.

# 6/27/19  
Transferred over code for the stochastic model implemented in `adaptivetau`

# 7/1/2019  
Added deterministic and stochastic versions of an age stratified model with non-linear man-to-snail force of infection as described in [Gurarie et al 2018](https://doi.org/10.1371/journal.pntd.0006514) and functions to simulate it. Added functionality to simulate models with or without events. Modified base model to include clumping parameter responsive to mean worm burden. Incorporated clumping parameter function derived from Senegal data relationship between log of mean worm burden and clumping parameter

# 7/2/2019  
Added explorations of mating probability over different values of the clumping parameter to assess the impact of variable $\kappa$ on the breakpoint. Initial simulations of age structured model. Some restructuring/renaming of variables  

# 7/3/2019  
Added `Model_animations.Rmd` which uses the model functions to produce simulations and visualizations as animations of worm burdens and $R_{eff}$ estimates through time over the course of different types of intervention efforts. Added a few functions/tweaks to existing functions to enable these sims/viz as well.

# 7/5/2019  
Worked on revamping and updating simulations with the simple model presented in [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128). New document `Garch_Mod.md` reflects these efforts. 

# 7/8/2019  
Worked on presentation for NIMBioS July 2019 talk. Realized that discrete time model isn't necessary to build transition and utility matrices for MDPtoolbox and extracting relevant values of $I_t$ and $W_t$ in the Garchitorena et al model at relevant time steps to construct these matrices actually gives way more flexibility in terms of implementing interventions

# 7/9/2019  
Ran simulations using differential equation model to identify optimal treatment allocations over values of $\mathcalP$, $R_0$, $T$, and $M$. Added function to simulate transmission through time with these inputs and return best treatment allocation $A_t$. Added density-dependent version of Garchitorena et al model with human to environment transmission a parabolic function of prevalence in the human population such that it decreases towards 0 at 0 prevalence and at 1 prevalence and peaks at 0.5 prevalence, simulating effects of both positive and negative density dependence

# 7/14/2019  
Finished running some simulations of $P(e)$ with the schisto model. Right now can only run the stochastic `adaptiveTau` model with permanent parameter reductions and single variable pulses like MDA. Need to spend some time making model more generalizable for easy implementation of different types/combinations of variables.  

Similar story with Garchitorena et al model. Since parameter interventions are implemented as forcing functions, the function passed to ode and the function simulating intervention decisions are hard coded for a particular intervention variable. Need to find a way to pass these options as function options

# 7/17-18/19  
Worked on deriving $R_{eff}$ from the basic schistosomiasis model, documented in the `Reff_derivation.Rmd` document. Spent lots of time just trying to simplify the expression down into something more analytically tractable and interpretable, but that's starting to feel like a bit of a losing battle.  

# 7/19/19  
Created outline for schisto $R_{eff}$ paper, did a little bit of coding for the age-stratified model  

# 7/22/19  
Did more work on the analystic expression of $R_{eff}$ and on the manuscript outline  

# 7/29/19  
Started writing the introduction and methods for the manuscript. Did some more scribbling to figure out different ways of expressing $R_{eff}$

# 7/30/19  
Worked on paper introduction and additional structure more. Finally settled on a final version of the model and corresponding $R_{eff}$ and $W_{bp}$ expressions, started working on coding these up.

# 7/31/19  
Worked on ESA presentation corresponding to this work  

# To-dos  
### Generalizable model structure for interventions in schisto stochastic model
### Generiaizable parameter intervention structure for Garchitorena et al model
### Derive $R_{eff}$ for age stratified model  
  This should be fairly straightforward since the force of infection in the snail population is a function of the total infectious input from the human population, so rather than the $M(W)$ term we would get the sum of infectious input over all populations considered
  
### Stochastic age stratified model  
### Fit base model and age-stratified models to data  
### Estimate BBR from Senegal human data  
### Document deterministic models (equations, base parameters) in vignettes  
#### Document derivation of $R_{eff}$  
