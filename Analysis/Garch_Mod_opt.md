Garchitorena et al model optimize across {*R*<sub>0</sub>, P} parameter space
================

-   [Garchitorena et al simplified model](#garchitorena-et-al-simplified-model)
    -   [Example Dynamics](#example-dynamics)
        -   [Without dd](#without-dd)
        -   [With dd](#with-dd)
-   [Optimal control framework](#optimal-control-framework)
    -   [Identify optimal control strategy with no density dependence](#identify-optimal-control-strategy-with-no-density-dependence)
    -   [Identify optimal control strategy with density dependence](#identify-optimal-control-strategy-with-density-dependence)
        -   [Plot optimal choices across variables](#plot-optimal-choices-across-variables)

Garchitorena et al simplified model
===================================

Using a simple, generalizable model of NTD transmission presented in [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128):

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

Example Dynamics
----------------

### Without dd

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
    theme_classic() + theme(legend.position = c(0.5,0.85)) +
    geom_line(size = 0.75) + 
  labs(x = "time (years)", 
       y = "Prevalence",
       title = "Model dynamics with annual MDA",
       subtitle = "with and without environmental intervention")
```

![](Garch_Mod_opt_files/figure-markdown_github/mod_example-1.png)

### With dd

``` r
#Get equilibrium estimates for state variables
garch_eq_dd <- runsteady(y = c(I = 0.5, W = 20), func = DDNTD::garch_mod_dd,
                      parms = DDNTD::garch_pars)[["y"]]

#PLot Reff curve
plot(seq(0,1,0.01), sapply(seq(0,1,0.01), 
                           function(I) garch_r0(garch_pars) * (4*I*(1-I))),
     xlab = "Prevalence (I)", ylab = expression(R[eff]), type = "l")
  lines(seq(0,1,0.01), sapply(seq(0,1,0.01), 
                              function(I) garch_r0(pars_env) * (4*I*(1-I))),
        col = 2)
  abline(h = 1, lty = 2)
```

![](Garch_Mod_opt_files/figure-markdown_github/dd_mod_reff-1.png)

``` r
#Run model with annual drug administration intervention
drug_only_dd <- sim_garch_mod(garch_eq_dd, run_time, garch_mod_dd, garch_pars,
                              events_df = drugs[1:3,])

#Run model with annual drug administration and environmental intervention 
drug_env_dd <- sim_garch_mod(garch_eq_dd, run_time, garch_mod_dd, pars_env,
                             events_df = drugs[1:3,])

#Plot results
rbind(drug_env_dd, drug_only_dd) %>% 
  mutate(Intervention = c(rep("Drug + Env", length(run_time)),
                          rep("Drug", length(run_time)))) %>% 
  ggplot(aes(x = time/365, y = I, lty = Intervention)) + 
    theme_classic() + theme(legend.position = c(0.5,0.85)) +
    geom_line(size = 0.75) + 
  labs(x = "time (years)", 
       y = "Prevalence",
       title = "Model dynamics with annual MDA and density dependence",
       subtitle = "with and without environmental intervention")
```

![](Garch_Mod_opt_files/figure-markdown_github/example_dd-1.png)

Optimal control framework
=========================

``` r
#Intervention costs, available capital, and conversion parameters
H <- 1000      # number of people
M <- H     # available capital, here modeled as $1.00 per person
  
#Parameters for the model and utility function: 
mdp_pars <- c(garch_pars, 
              "d" = 100,      # Arbitrary cost of having prevalence of I_t
              "M" = M,        # Capital available to spend on MDA
              "theta" = 0.8/M, # Scaling of capital spent on MDA to reduction in prevalence, 
              "mu" = 0.1*(0.8/M),    # Scaling of environmental intervention to cost (impact per cost)
              "delta" = 0.95)  #Arbitrary discounting rate  

test_grid <- expand.grid(R0 = seq(1.5, 10, 0.5),
                         script_P = exp(seq(log(0.001), 
                                            log(1), 
                                            length.out = 10)),
                         T_frame = c(5, 10, 20)*365,
                         freq = 365,
                         M = M,
                         A = seq(0,1,0.05))
```

``` r
make_force_fx <- function(t_max, burn_in, freq, par_set, par, par_effect){
  force_fx <- approxfun(as.data.frame(list(
    times = c(1:t_max),
    par_val = c(rep(par_set[par], times = burn_in),
                sapply(c(1:(t_max/freq - 1)),
                       function(t) rep(par_set[par]*(par_effect^t), 
                                       freq)))
  )), rule = 2)
  
  return(force_fx)
  
}
garch_mod_forceV = function(t, n, p){
  beta_e <- p["beta_e"]
  beta_d <- p["beta_d"]
  gamma <- p["gamma"]
  omega <- p["omega"]
  V <- V_force_fx(t)
  sigma <- p["sigma"]
  lambda <- p["lambda"]
  rho <- p["rho"]
  
    I = n[1]
    W = n[2]

    dIdt = (beta_e*W + beta_d*I)*(1-I) - gamma*I
    dWdt = omega + V*sigma*lambda*I - rho*W

    return(list(c(dIdt, dWdt)))
}

garch_mod_dd_forceV = function(t, n, p){
  beta_e <- p["beta_e"]
  beta_d <- p["beta_d"]
  gamma <- p["gamma"]
  omega <- p["omega"]
  V <- V_force_fx(t)
  sigma <- p["sigma"]
  lambda <- p["lambda"]
  rho <- p["rho"]
  
    I = n[1]
    W = n[2]

    dIdt = (beta_e*W + beta_d*I)*(1-I) - gamma*I
    dWdt = omega + V*sigma*lambda*I*(4*I*(1-I)) - rho*W

    return(list(c(dIdt, dWdt)))
}

sim_w_decision <- function(R0, script_P, 
                           T_frame, freq, 
                           M, A, 
                           pars, par_target, 
                           base_mod, int_mod){
  
  pars["beta_e"] <- (pars["gamma"] * pars["rho"]) * R0
  pars["mu"] <- pars["theta"]*script_P

  #Get reduction in prevalence based on capital and allocation towards MDA
  M_A_I <- 1 - M*A*pars["theta"]      #Total percent reduction in I
  M_A_par <- 1 - M*(1-A)*pars["mu"]       #per year reduction in target par

  #Create events dataframes based on capital allocation decisions
  I_events <- data.frame(var = rep("I", times = round(T_frame/freq)),
                         time = c(1:round(T_frame/freq))*freq,
                         value = rep(M_A_I, time = round(T_frame/freq)),
                         method = rep("mult", times = round(T_frame/freq)))
  
  #generate forcing function for target parameter intervened on
  V_force_fx <<- make_force_fx(t_max = T_frame,
                              burn_in = freq,
                              freq = freq,
                              par_set = pars,
                              par = par_target,
                              par_effect = M_A_par)
  
  eq_vals <- runsteady(y = c(I = 0.5, W = 20), func = base_mod,
                       parms = pars)[["y"]]

  sim <- as.data.frame(ode(y = eq_vals, 
                           times = c(1:T_frame), 
                           func = int_mod, 
                           parms = pars,
                           events = list(data = I_events)))

  return(sim)

}

sim_w_a_get_u <- function(R0, script_P, 
                          T_frame, freq, 
                          M, A, 
                          pars, par_target, 
                          base_mod, int_mod){
  sim <- sim_w_decision(R0, script_P, 
                        T_frame, freq, 
                        M, A, 
                        pars, par_target, 
                        base_mod, int_mod)
  U <- -sum((sim %>% 
              filter(time > freq) %>% 
              pull(I) *pars["d"])^1.5)
  
  return(Utility = U)
}
```

``` r
test_sim_no_dd <- sim_w_decision(R0 = 4, script_P = 0.1,
                            T_frame = 3650,
                            freq = 365,
                            M = M,
                            A = 1,
                            pars = mdp_pars,
                            par_target = "V",
                            base_mod = garch_mod,
                            int_mod = garch_mod_forceV) %>% 
  mutate("Dens_Dep" = "No")

test_sim_dd <- sim_w_decision(R0 = 4, script_P = 0.1,
                               T_frame = 3650,
                               freq = 365,
                               M = M,
                               A = 1,
                               pars = mdp_pars,
                               par_target = "V",
                               base_mod = garch_mod_dd,
                               int_mod = garch_mod_dd_forceV) %>% 
  mutate("Dens_Dep" = "Yes")

  rbind(test_sim_no_dd, test_sim_dd) %>% 
    ggplot(aes(x = time, y = I, col = Dens_Dep)) +
      geom_line(size = 1.2) +
      theme_classic() +
      scale_x_continuous(breaks = c(0:(3650/365))*365,
                         labels = c(0:(3650/365))) +
      theme(legend.position = c(0.8,0.8)) +
      labs(x = "Time (years)",
           y = "Prevalence (I)",
           title = "Influence of PDD on effect of intervention",
           col = "PDD")
```

![](Garch_Mod_opt_files/figure-markdown_github/test_decision_dd-1.png)

Identify optimal control strategy with no density dependence
------------------------------------------------------------

``` r
plan(multiprocess)
#tic()
test_Us <- future_pmap_dbl(as.list(test_grid), sim_w_a_get_u, 
                           pars = mdp_pars, par_target = "V",
                           base_mod = garch_mod, int_mod = garch_mod_forceV)
#toc()
```

Identify optimal control strategy with density dependence
---------------------------------------------------------

``` r
plan(multiprocess)
#tic()
test_Us_dd <- future_pmap_dbl(as.list(test_grid), sim_w_a_get_u, 
                              pars = mdp_pars, par_target = "V",
                              base_mod = garch_mod_dd, int_mod = garch_mod_dd_forceV)
#toc()
```

### Plot optimal choices across variables

``` r
cbind(test_grid, Utility = test_Us) %>% 
  group_by(R0, script_P, T_frame, M) %>% 
  summarise(max_u = max(Utility),
            A_opt = A[which(Utility == max_u)]) %>% 
  mutate(C_labs = paste0("C = ", M),
         T_labs = paste0("T = ", T_frame/365),
         T_labs = factor(T_labs, levels = c("T = 5", "T = 10", "T = 20"))) %>% 
  ggplot(aes(x = R0, y = script_P, z = A_opt)) +
    facet_grid(. ~ T_labs) +
    geom_tile(aes(fill = A_opt)) +
    theme_classic() +
    scale_fill_viridis(option = "magma", direction = 1) +
    scale_x_continuous(breaks = c(1.5, 3, 6, 9)) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.001,0.01,0.1,1),
                       labels = c(0.001,0.01,0.1,1)) +
    labs(x = expression(Transmission~Intensity~(italic(R[0]))),
         y = expression(Relative~Intervention~Cost~(italic(Rho))),
         fill = expression(Intervention~Allocation~(italic(A))))
```

![](Garch_Mod_opt_files/figure-markdown_github/opt_no_dd_plot-1.png)

``` r
cbind(test_grid, Utility = test_Us_dd) %>% 
  group_by(R0, script_P, T_frame, M) %>% 
  summarise(max_u = max(Utility),
            A_opt = A[which(Utility == max_u)]) %>% 
  mutate(C_labs = paste0("C = ", M),
         T_labs = paste0("T = ", T_frame/365),
         T_labs = factor(T_labs, levels = c("T = 5", "T = 10", "T = 20"))) %>% 
  ggplot(aes(x = R0, y = script_P, z = A_opt)) +
    facet_grid(. ~ T_labs) +
    geom_tile(aes(fill = A_opt)) +
    theme_classic() +
    scale_fill_viridis(option = "magma", direction = 1) +
    scale_x_continuous(breaks = c(1.5, 3, 6, 9)) +
    scale_y_continuous(trans = "log",
                       breaks = c(0.001,0.01,0.1,1),
                       labels = c(0.001,0.01,0.1,1)) +
    labs(x = expression(Transmission~Intensity~(italic(R[0]))),
         y = expression(Relative~Intervention~Cost~(italic(Rho))),
         fill = expression(Intervention~Allocation~(italic(A))))
```

![](Garch_Mod_opt_files/figure-markdown_github/opt_dd_plot-1.png)
