Optimal control of neglected tropical diseases using a simple dynamic model
================
Chris Hoover
May 18, 2018

------------------------------------------------------------------------

Base case: Garchitorena et al simple model
==========================================

Using a simple, generalizable model of NTD transmission presented in [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128), two interventions will be considered: 1) drug administration, implemented as a pulse reduction in the state variable, *I*, that reduces the prevalence of the disease in the population and 2) environmental remediation (e.g. improvement in sanitation, vector control), implemented as a permanent alteration of a model parameter, that reduces the transmission of the disease.

The model consists of two state variables, *I* and *W*, that correspond to the prevalence of infection in the human population (i.e. proportion infected at time=*t*) and the degree of contamination of the environment with the disease-causing agent, respectively. We make the simplifying assumption that individuals can only be susceptible, *S*, or infected, *I*, meaning *S* + *I* = 1 and eliminating the need for a recovered *R* compartment as is typical of SIR models but would complicate things here. The model equations then are:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdI%7D%7Bdt%7D%3D(%5Cbeta_EW+%5Cbeta_DI)(1-I)-%5Cgamma%20I)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5COmega+V%5Csigma%5Clambda%20I-%5Crho%20W)

For most environmentally mediated infectious diseases, transmission to humans is exclusively from the environment and there is no exogenous production of infectious agents in the environment (i.e. *Ω* = 0 and *β*<sub>*D*</sub> = 0. We can derive *R*<sub>0</sub> for this system quite simply since transmission is only determined by the environmental component:

![](http://latex.codecogs.com/gif.latex?R_0%3D%5Cfrac%7BV%5Csigma%5Clambda%5Cbeta_E%7D%7B%5Cgamma%5Cdelta%7D)

with parameter definitions and values in table 1, we have *R*<sub>0</sub>≈ 4

|                   | Value   | Description                                                       |
|:------------------|:--------|:------------------------------------------------------------------|
| *β*<sub>*E*</sub> | 4e-05   | Transmission rate from environment to human population            |
| *β*<sub>*D*</sub> | 0       | Human to human transmission rate                                  |
| *γ*               | 0.00091 | Rate of recovery from infected back to susceptible                |
| *Ω*               | 0       | Recruitment rate of infectious agents in the environment          |
| *V*               | 1       | Abundance of vectors/intermediate hosts/suitable environment      |
| *λ*               | 1       | Recruitment rate of infectious agents by infectious individuals   |
| *σ*               | 1       | Fraction of infectious agents produced that reach the environment |
| *ρ*               | 0.01111 | Mortality rate of infectious agents in the environment            |

Effects of parameter changes on *R*<sub>0</sub>
-----------------------------------------------

Even within this simple model, we have a alot of options for simulating different interventions. Drug mased treatments in the human population can be implemented as instantaneous reductions in *I*, improvements in sanitation or education that reduces environmental contamination can be modeled as reductions in variable *σ*, education may also reduce transmission from environment to people (e.g. people learn to avoid exposure) modeled as reductions in *β*<sub>*E*</sub>, vector control or intermediate host control can be modeled as a temporary reduction in parameter *V* or as an instantaneous reduction in state variable *W*. Let's see how some of these parameter changes affect *R*<sub>0</sub>

``` r
#vector of relative parameter reductions
rel_par_reds <- seq(0,1,0.05)

#Function to estimate r0 with change in parameter value
est_r0_change_par <- function(pars, par_change, change_amount){
  pars_use <- pars
  pars_use[par_change] <- pars_use[par_change]*(1-change_amount)
  
  garch_r0(pars_use)
}

#Estimate r0 with relative changes in parameter reductions
garch_r0_pars <- data.frame(rel_par = rel_par_reds,
                            r0_beta_e = map_dbl(rel_par_reds, est_r0_change_par, 
                                                pars = garch_pars,
                                                par_change = "beta_e"),
                            r0_sigma = map_dbl(rel_par_reds, est_r0_change_par, 
                                               pars = garch_pars,
                                               par_change = "sigma"),
                            r0_V = map_dbl(rel_par_reds, est_r0_change_par, 
                                           pars = garch_pars,
                                           par_change = "beta_e"))

garch_r0_pars %>% 
  gather("Parameter Change", "R0", r0_beta_e:r0_V) %>% 
  ggplot(aes(x = rel_par, y = R0, col = `Parameter Change`)) +
    geom_line() +
    theme_classic()
```

![](Garch_Mod_files/figure-markdown_github/int_r0-1.png)

Should've seen that one coming actually, since all of these parameters are in the same place in the *R*<sub>0</sub> expression, relative changes in each will have the same affect on *R*<sub>0</sub>. So the true question here is what is the cost of reducing each?

MDA Intervention
----------------

Now want to investigate dynamics under a routine MDA campaign with and without an accompanied environmental intervention such as a reduction in the vector/intermediate host population

``` r
#Get equilibrium estimates for state variables
garch_eq <- runsteady(y = c(I = 0.5, W = 20), func = DDNTD::garch_mod,
                      parms = DDNTD::garch_pars)[["y"]]

#Events representing interventions to implement in the model
n_years <- 20                  #Simulate 10 years of annual intervention
drug_efficacy <- 0.75            #75% reduction in prevalence at each drug treatment
run_time <- c(1:(365*n_years)) #Run model for 10 years

#data frame with event times for drug administration
drugs <- data.frame(var=rep('I', times = n_years/2),
                    time = c(1:(n_years/2))*365,
                    value = rep((1-drug_efficacy), times = n_years/2),
                    method = rep("mult", times = (n_years/2)))
    
#Run model with annual drug administration intervention
drug_only <- sim_garch_mod(garch_eq, run_time, garch_mod, garch_pars,
                           events_df = drugs)

#Reduce vector/intermediate host population
pars_env <- garch_pars
pars_env["V"] <- garch_pars["V"] * 0.8 

#Run model with annual drug administration and environmental intervention 
drug_env <- sim_garch_mod(garch_eq, run_time, garch_mod, pars_env,
                           events_df = drugs)

#Plot results
rbind(drug_env, drug_only) %>% 
  mutate(Intervention = c(rep("Drug + Env", length(run_time)),
                          rep("Drug", length(run_time)))) %>% 
  ggplot(aes(x = time/365, y = I, lty = Intervention)) + 
    theme_bw() + theme(legend.position = c(0.5,0.85)) +
    geom_line(size = 0.75) + xlab("time (years)") + ylab("Prevalence") + ylim(c(0,0.75))
```

![](Garch_Mod_files/figure-markdown_github/mod_prep-1.png)

Clear that the combined intervention performs better, but in a resource constrained environment, how do we determine where to allocate resources?

Optimal control framework
=========================

Translate the model into a Markov Decision Process (MDP) solvable with stochastic dynamic programming (SDP)
-----------------------------------------------------------------------------------------------------------

Next we want to translate our system of continuous time differential equations into a Markov Decision Process (MDP) consisting of **1)** a Markov chain in which the state of the system at time *t* + 1 is dependent only on the current state of the system (i.e. the state at *t*) and **2)** a decision or control action that is being made at each state transition (i.e. from *t* to *t* + 1). We therefore want a single, discrete time equation that captures about the same dynamics of our simple disease system.

We'll start with the assumption that the dynamics of the infectious agents in the environment are faster than the dynamics of the prevalence in the human population, therefore they reach a steady state equilibrium:

![](http://latex.codecogs.com/gif.latex?W%5E*%3D%5Cfrac%7BV%5Csigma%5Clambda%20I%7D%7B%5Crho%7D)

Again, substituting we get:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdI%7D%7Bdt%7D%3D%5CBig(%5Cfrac%7B%5Cbeta_EV%5Csigma%5Clambda%20I%7D%7B%5Crho%7D%5CBig)%5CBig(1-I%5CBig)-%5Cgamma%20I)

which can be translated to a discrete time model as:

![](http://latex.codecogs.com/gif.latex?I_%7Bt+1%7D%3DI_t+%5CBig(%5Cfrac%7B%5Cbeta_EV%5Csigma%5Clambda%20I%7D%7B%5Crho%7D%5CBig)%5CBig(1-I_t%5CBig)-%5Cgamma%20I_t)

let's check to make sure that the dynamics of the continuous time model hold in this simplified, discrete time version of the model.

``` r
discrete_sim <- sim_discrete_mod(I_0 = garch_eq["I"],
                                 time = max(run_time),
                                 parameters = garch_pars,
                                 events_df = drugs) %>% 
  mutate(Intervention = "Drug")

discrete_env_sim <- sim_discrete_mod(I_0 = garch_eq["I"],
                                     time = max(run_time),
                                     parameters = pars_env,
                                     events_df = drugs) %>% 
  mutate(Intervention = "Drug + Env")

bind_rows(discrete_sim, discrete_env_sim) %>% 
  ggplot(aes(x = t, y = I, lty = Intervention)) +
    theme_classic() + 
    theme(legend.position = c(0.5,0.85)) +
    geom_line(size = 1.2) + 
    labs(x = "time", y = "Prevalence") +
    ylim(c(0,1))
```

![](Garch_Mod_files/figure-markdown_github/discrete_time-1.png)

Now we'll work through the "Six steps of stochastic dynamic programming" as described in [Marescot et al](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12082) to get our MDP and set ourselves up for optimization

1.  **Define the optimization objective of the problem.**
    Our objective is to minimize the costs, *C*, associated with 1) infection of individuals and 2) implementing interventions to reduce their infections

2.  **Define the set of states that represent the configuration of the system at time *t***
    We define the state variable, *X*<sub>*t*</sub>, as the prevalence of infection in the human population, *I*<sub>*t*</sub> ∈ \[0, 1\].

3.  **Define the decision variable, *A*<sub>*t*</sub> that is the component of the system to be controlled to meet the objective**
    We define the decision variable *A*<sub>*t*</sub> as the proportion of capital, *M*, committed to drug administration with the rest allocated towards environmental interventions.

4.  **Build a transition model describing the system's dynamics as a function of the decision variable and the system state in the prior time step**
    Beginning with the discrete time model above, we have:

![](http://latex.codecogs.com/gif.latex?I_%7Bt+1%7D%3DI_t+%5CBig(%5Cfrac%7B%5Cbeta_EV%5Csigma%5Clambda%20I%7D%7B%5Crho%7D%5CBig)%5CBig(1-I_t%5CBig)-%5Cgamma%20I_t)

Incorporating interventions based on our decision variable, *A*<sub>*t*</sub> gives us:

![](http://latex.codecogs.com/gif.latex?I_%7Bt+1%7D%3DI_t+%5CBig(%5Cfrac%7B%5Cbeta_EV%5Csigma%5Clambda%20IM(1-A_t)%5Cmu%7D%7B%5Crho%7D%5CBig)%5CBig(1-I_t%5CBig)-I_t(%5Cgamma+MA_t%5Ctheta))

Where *θ* is a constant that converts capital, *M*, to units of prevalence and *μ* is a constant that converts *M* to the same units as the parameter being intervened on with possibilities being *β*<sub>*E*</sub>, *V*, *σ*.

1.  **Define the utility function *U*<sub>*t*</sub>(*X*<sub>*t*</sub>, *A*<sub>*t*</sub>) representing the desirability of acting in a given state**

![](http://latex.codecogs.com/gif.latex?%5Cmax_%7BA%7DC%3D%5Cmax_%7BA%7D%5Csum_0%5ET%5Cfrac%7B%7B-%5CPi(I_t)%7D%7D%7B%5Cdelta_t%7D)

where *Π*(*I*<sub>*t*</sub>)=*d**I*<sub>*t*</sub> + *M*, *d* is the cost associated with having prevalence of *I*<sub>*t*</sub> and *M* is the capital spent on intervention. We assume capital spent each year is the same and therefore only consider the cost of infection, *d**I*<sub>*t*</sub>. We're also maximizing the negative costs, aka minimizing the costs here.

1.  **Determine the optimal solution of the optimization problem**
    Step 6 is Part 3 below

Find the optimal solution using SDP in the `MDPtoolbox` R package
-----------------------------------------------------------------

Let's translate the problem outlined above into code and also define some new parameters that control how the interventions act in the model. For the environmental intervention, we'll parameterize the effect based on the conversion parameter, *μ*, depending on the parameter we choose to target (*β*<sub>*E*</sub>, *V*, *σ*), therefore we actually have three different parameters, *μ*<sub>*β*<sub>*E*</sub></sub>, *μ*<sub>*V*</sub>, *μ*<sub>*σ*</sub>. These parameters simply convert capital allocated towards the intervention into the same units as the parameter in the model on which the intervention acts (i.e. converts capital to units of snail mortality or snail carrying capacity). We'll also need a paramter that does this for the MDA intervention, we'll assign *θ* for this purpose.

``` r
#Vector of all possible states, i.e. prevalence ranging from 0 to 1  
states <- seq(0, 1, 0.01)

#Vector of possible actions, A_t, ranging from total investment in drugs (M=1) 
#  to total investment in environmental intervention (M=0)  
decisions <- seq(0, 1, 0.01)

#Intervention costs, available capital, and conversion parameters
H <- 1000      # number of people
M <- 1.5*H     # available capital, here modeled as $1.50 per person

theta <- H*0.01*1.8  #cost of 1% reduction in prevalence (if M goes towards reducing prevalence, 80% reduction)
mu <- theta        #cost of 1% reduction in W*
  
  
#Parameters for the model and utility function: 
mdp_pars <- c(garch_pars, 
              "H" = H,        #Number of people in community
              "d" = 100,      # Arbitrary cost of having prevalence of I_t
              "M" = M,        # Capital available to spend on MDA
              "theta" = theta, # Scaling of capital spent on MDA to reduction in prevalence, 
              "mu" = mu,    # Scaling of environmental intervention to cost (impact per cost)
              "delta" = 0.95)  #Arbitrary discounting rate  

#Initialize transition matrix across I states and A actions
transition <- array(0, dim = c(length(states), length(states), length(decisions)))

# Initialize utility matrix
utility <- array(0, dim = c(length(states), length(decisions)))
```

``` r
test_choices <- as.data.frame(sapply(c(0,0.5,1), sim_int_choice,
                                     I_0 = garch_eq["I"],
                                     int_par = "V",
                                     time = 365,
                                     parameters = mdp_pars)) %>% 
  mutate("time" = c(1:367))
colnames(test_choices) <- c("A0", "A0.5", "A1", "time")

test_choices %>% 
  gather("A", "Val", A0:A1) %>% 
  ggplot(aes(x = time, y = Val, col = A)) +
    theme_classic() +
    labs(x = "time (days)",
         y = "prevalence",
         title = "Prevalence trajectories after different decisions (A)") +
    geom_line(size = 1.2)
```

![](Garch_Mod_files/figure-markdown_github/opt_test_viz-1.png)

``` r
# Fill in the transition and utility matrices
# Loop on all states 
for (i in 1:length(states)) {

    # Loop on all actions
    for (a in 1:length(decisions)) {

# Calculate the transition state at the next step, given the current state i and the intervention, a
        nextpop <- round(sim_int_choice(A_t = decisions[a], 
                                        I_0 = states[i],
                                        int_par = "V",
                                        time = 365,
                                        parameters = mdp_pars)[365+2], 2)
        nextind <- pmatch(nextpop, states) #replace which[] with pmatch to ensure proper indexing
        
#Since this model is deterministic, assign probability of 1 to value of I_t+1, everything else =0       
        transition[i, nextind, a] <- 1

# Compute utility as exponentially increasing function of prevalence
    # Implying higher prevalence values are increasingly worse
        utility[i,a] <- -1*(mdp_pars["d"]*nextpop)^1.5

    } # end of action loop
} # end of state loop

#Check that transition and utility matrices are valid
mdp_check(transition, utility)
```

    ## [1] ""

Now we have the transition and utility matrices. The transition matrix indicates for a given starting prevalence, *I*<sub>*t*</sub>, and action *A*<sub>*t*</sub>, the probability of the state at *I*<sub>*t* + 1</sub>. Since this model is deterministic, the transition matrix contains a whole bunch of 0s and some 1s that indicate the value *I*<sub>*t* + 1</sub>|*I*<sub>*t*</sub>, *A*<sub>*t*</sub>. Future models will incorporate stochasticity in transmission to spread the probability around. The utility matrix simply represents the utility of each action, *A*<sub>*t*</sub> applied to each prevalence value, *I*<sub>*t*</sub>.

We then use the `MDPtoolbox::mdp_finite_horizon` function to identify the optimal policy for each starting condition. This function returns a list with three elements:
\* The valuation matrix, **V**, in which rows represent starting states, columns represent the time steps at which actions occur, and cell values indicate the value of the objective function
\* The policy matrix, **P**, in which rows again represent starting states, columns again represent the time steps at which actions occur, and cell values indicate the *index* of the decision vector (i.e. if the cell value is 50, the optimal decision variable at time step \[column\] for starting state \[row\] is the 50th entry of the decision vector)

``` r
#optimize over 10 years
t_opt <- 10

#Get optimal intervention strategy given t_opt
opt1 <- mdp_finite_horizon(transition, utility, 
                           mdp_pars["delta"], t_opt)

#Implement example optimal intervention
opt_p50 <- sim_opt_choice(A_vec = decisions,
                          opt_list = opt1,
                          states = states,
                          p_start = 75,
                          t_per_step = 365,
                          int_par = "V",
                          parameters = mdp_pars)

opt_p50 %>% 
  ggplot(aes(x = time, y = I)) +
    geom_line() + 
    theme_classic() +
    ylim(c(0,1))
```

![](Garch_Mod_files/figure-markdown_github/get_opt_finite-1.png)

``` r
#Get optimal intervention strategy using value iteration
opt2 <- mdp_value_iteration(transition, utility, 
                            discount = mdp_pars["delta"], epsilon = 0.001)
```

    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"

``` r
opt2["policy"]
```

    ## $policy
    ##   [1]  99  78  99  90  98  92  98  94  98 101  98 100  97 100  97  99 101
    ##  [18]  99 100  98 100 101  99 100 101 100 101  99 100 101 100 100 101 100
    ##  [35] 101 101 100 101 101 100 101 101 100 101 101 100 101 101 100 101 101
    ##  [52] 100 101 101 100 101 101 100 101 101 101 100 101 101 100 101 101 101
    ##  [69] 101 101 101 100 101 101 101 100 101 101 101 101 101 101 100 101 101
    ##  [86] 101 100 101 101 101 100 101 101 101 101 101 101 101 101 101

Explore optimal interventions across values of *R*<sub>0</sub> and 𝒫
--------------------------------------------------------------------

``` r
opt_pol_fx <- function(beta_e, script_P, pars, decisions, states){
  pars["beta_e"] <- beta_e
  pars["mu"] <- pars["theta"]*script_P
  
  r0 <- garch_r0(pars)
  
  eq <- runsteady(y = c(I = 0.5, W = 20), func = DDNTD::garch_mod,
                  parms = pars)[["y"]]
  
  #Initialize transition matrix across I states and A actions
  trans_mat <- array(0, dim = c(length(states), length(states), length(decisions)))
  
  # Initialize utility matrix
  util_mat <- array(0, dim = c(length(states), length(decisions)))
  
  
  for (i in 1:length(states)) {
    
    # Loop on all actions
    for (a in 1:length(decisions)) {
      
      # Calculate the transition state at the next step, given the current state i and the intervention, a
      nextpop <- round(sim_int_choice(A_t = decisions[a], 
                                      I_0 = states[i],
                                      int_par = "V",
                                      time = 365,
                                      parameters = pars)[365+2], 2)
      nextind <- pmatch(nextpop, states) #replace which[] with pmatch to ensure proper indexing
      
      #Since this model is deterministic, assign probability of 1 to value of I_t+1, everything else =0     
      trans_mat[i, nextind, a] <- 1
      
      # Compute utility as exponentially increasing function of prevalence
      # Implying higher prevalence values are increasingly worse
      util_mat[i,a] <- -1*(pars["d"]*nextpop)^1.5
      
    } # end of action loop
  } # end of state loop
  
  #Check that transition and utility matrices are valid
  if(mdp_check(trans_mat, util_mat) != ""){
    
    stop("Transition or utility matrix error")
    
  } else {
    #Get optimal intervention strategy using value iteration
    opt_pols <- mdp_value_iteration(transition, utility, 
                                    discount = pars["delta"], epsilon = 0.001)
    
    opt_pol <- opt_pols[["policy"]][round(eq[["I"]],2)*100]
    
  }
  
  return(c("R0" = r0, 
           "I_eq" = eq[["I"]], 
           "W_eq" = eq[["W"]], 
           "A_t" = decisions[opt_pol]))
}

opt_pol_df <- as.data.frame(expand.grid("beta_e" = seq(mdp_pars["gamma"] * mdp_pars["rho"]+1e-6, 
                                                       (mdp_pars["gamma"] * mdp_pars["rho"])*20,
                                                       length.out = 10),
                                        "script_P" = seq(0.01, 100, length.out = 10)))

#test_opt <- opt_pol_fx(opt_pol_df$beta_e[10], opt_pol_df$script_P[10], mdp_pars, decisions, states)

opt_pol_slns <- mapply(opt_pol_fx, opt_pol_df$beta_e, opt_pol_df$script_P, 
                       MoreArgs = list(pars = mdp_pars,
                                       decisions = decisions,
                                       states = states))
```

    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
    ## [1] "MDP Toolbox: iterations stopped, epsilon-optimal policy found"
