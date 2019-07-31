Age structured schistosomiasis model
================

Age structured model
====================

``` r
H = 300
area = 200

age_strat_pars <- c(
  ##standard snail parameters
    f_N=0.10, # recruitment rate (from sokolow et al)
    C=50*area, # carrying capacity corresponding to 50 snails per square meter
    mu_N=1/60, #Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=1/40, #Transition rate from exposed to infected (assuming pre-patency period of 40 days)
    mu_I=1/10 - 1/60, # Increased mortality rate of infected snails

  #Adult Worm, Miracidia and Circariae Parameters
    mu_W = 1/(3.3*365), # death rate of adult worms
    m = 5.2,            # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = 0.08,           # mean egg viability of eggs shed into environment

  #Density dependence parameters
    gamma = 5e-3,       # parameter of fecundity reduction function
    xi = 2.8e-3,        # parameter for acquired immunity funtion

  #Human parameters
    H = H,          # Total number of people
    prop_SAC = 0.5,       #Percent of people that are school age children (SAC)
    prop_adult = 0.5, #percent of people that are not school age children (assumed here to be adults)
    cvrg_SAC = 0.9,  #MDA coverage in SAC population
    cvrg_adult = 0,  #MDA coverage in adult population
    u_SAC = 50, # mL urine per SAC/day/10mL assumed to be half of adult (approximate, ranges from 20 - 100)
    u_adult = 100, # mL urine per adult/day/10mL (approximate, ranges from 80 - 200)
    rho_SAC = 0.3, # relative number of eggs shed by SAC that make it to snail habitat (~sanitation)
    rho_adult = 0.015, # relative number of eggs shed by adults that make it to snail habitat (~sanitation); 5% of SAC (truscott et al)
    omega_SAC = 1,     # relative infection risk of SAC (related to clean water access/education/water contact)
    omega_adult = 0.1,   # relative infection risk of adults (related to clean water access/education/water contact)

  #Transmission parameters
    lambda=1.2e-3, # snail-to-man transmission
    beta=1.6e-3  # man-to-snail transmission
)

schisto_age_strat_mod <- function(t, n, parameters) {
  with(as.list(parameters),{

    S=n[1]
    E=n[2]
    I=n[3]
    Wt_SAC=n[4]
    Wu_SAC=n[5]
    Wt_adult=n[6]
    Wu_adult=n[7]

    #Total snail population
      N=S+E+I

    #weighting treated and untreated populations among SAC
      W_SAC = (cvrg_SAC*Wt_SAC) + ((1-cvrg_SAC)*Wu_SAC)

    #weighting treated and untreated populations among adults
      W_adult = (cvrg_adult*Wt_adult) + ((1-cvrg_adult)*Wu_adult)

    #Weighting SAC and adult populations
      W_tot = W_SAC*prop_SAC + W_adult*prop_adult

    #Update clumping parameter, k from estimate of eggs burden per 10mL estimate
      k_t_SAC = k_from_log_W(Wt_SAC)
      k_u_SAC = k_from_log_W(Wu_SAC)
      k_t_adult = k_from_log_W(Wt_adult)
      k_u_adult = k_from_log_W(Wu_adult)

    #Estimate mating probability within each strata
      phi_Wt_SAC = phi_Wk(W = Wt_SAC, k = k_t_SAC)  #Mating probability in treated SAC population
      phi_Wu_SAC = phi_Wk(W = Wu_SAC, k = k_u_SAC)  #Mating probability in untreated SAC population
      phi_Wt_adult = phi_Wk(W = Wt_adult, k = k_t_adult)  #Mating probability in treated adult population
      phi_Wu_adult = phi_Wk(W = Wu_adult, k = k_u_adult)  #Mating probability in untreated adult population

    #Estimate mean eggs produced per person in each strata as product of
    # mated female worms,
    # eggs produced per female worm per 10mL urine,
    # and reduction in fecundity due to crowding
      eggs_Wt_SAC = 0.5*(Wt_SAC*phi_Wt_SAC) * m * f_Wgk(Wt_SAC, gamma, k_t_SAC)

      eggs_Wu_SAC = 0.5*(Wu_SAC*phi_Wu_SAC) * m * f_Wgk(Wu_SAC, gamma, k_u_SAC)

      eggs_Wt_adult = 0.5*(Wt_adult*phi_Wt_adult) * m * f_Wgk(Wt_adult, gamma, k_t_adult)

      eggs_Wu_adult = 0.5*(Wu_adult*phi_Wu_adult) * m * f_Wgk(Wu_adult, gamma, k_u_adult)

    #Estimate miracidia produced by each strata as product of
    # mean eggs per 10 mL urine for individuals in each strata,
    # number of people in each strata,
    # egg viability
    # mean mL urine produced by an average individual in each group/10,
    # contamination coefficient for SAC/adults,

      M_Wt_SAC = eggs_Wt_SAC * ((H*prop_SAC)*cvrg_SAC) * v * u_SAC * rho_SAC

      M_Wu_SAC = eggs_Wu_SAC * ((H*prop_SAC)*(1-cvrg_SAC)) * v * u_SAC * rho_SAC

      M_Wt_adult = eggs_Wt_adult * ((H*prop_adult)*cvrg_adult) * v * u_adult * rho_adult

      M_Wu_adult = eggs_Wu_adult * ((H*prop_adult)*(1-cvrg_adult)) * v * u_adult * rho_adult

    # Estimate total miracidia entering snail habitat
      M_tot = M_Wt_SAC + M_Wu_SAC + M_Wt_adult + M_Wu_adult
      
      #if(t %% 100 == 0) print(1-exp(-M_tot/N))

    # Snail infection dynamics
      dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*(1-exp(-M_tot/N))*S #Susceptible snails

      dEdt= beta*(1-exp(-M_tot/N))*S - (mu_N+sigma)*E #Exposed snails

      dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human populations
      dWt_SACdt= (omega_SAC*lambda*I*R_Wv(Wt_SAC, xi)) - (mu_W*Wt_SAC)
      dWu_SACdt= (omega_SAC*lambda*I*R_Wv(Wu_SAC, xi)) - (mu_W*Wu_SAC)
      dWt_adultdt= (omega_adult*lambda*I*R_Wv(Wt_adult, xi)) - (mu_W*Wt_adult)
      dWu_adultdt= (omega_adult*lambda*I*R_Wv(Wu_adult, xi)) - (mu_W*Wu_adult)

    return(list(c(dSdt,dEdt,dIdt,
                  dWt_SACdt,dWu_SACdt,
                  dWt_adultdt, dWu_adultdt)))
  })
}
```

``` r
years <- 20
age_time <- c(1:(365*years))

age_start <- c(S=5000, E=0, I=0, Wt_SAC=10, Wu_SAC=10, Wt_adult=10, Wu_adult=10)

#Run to equibrium with base parameter set
age_eqbm <- runsteady(y = age_start, func = schisto_age_strat_mod,
                      parms = age_strat_pars)[["y"]]

schisto_age_sim <- sim_schisto_mod(nstart = age_eqbm, 
                                   time = age_time, 
                                   model = schisto_age_strat_mod,
                                   parameters = age_strat_pars,
                                   events_df = NA)
```

Simulate MDA in 75% of SAC populaition
--------------------------------------

``` r
eff <- 0.93
age_strat_pars["cvrg_SAC"] <- 0.75

sac_mda <- data.frame(var = rep('Wt_SAC', years/2),
                      time = c(1:(years/2))*365,
                      value = rep((1 - eff), years/2),
                      method = rep('mult', years/2))

schisto_mda_sim <- sim_schisto_mod(nstart = age_eqbm, 
                                   time = age_time, 
                                   model = schisto_age_strat_mod,
                                   parameters = age_strat_pars,
                                   events_df = sac_mda)
```

    ## Warning in if (is.na(events_df)) {: the condition has length > 1 and only
    ## the first element will be used

``` r
schisto_mda_sim %>% 
  gather("Treatment", "Worm Burden",  Wt_SAC:Wu_adult) %>% 
  ggplot(aes(x = time, y = `Worm Burden`, lty = Treatment)) +
    geom_line() +
    theme_classic() +
    theme(legend.position = c(0.6,0.8)) +
    ylim(c(0,100)) +
    scale_x_continuous(breaks = c(0:years)*365,
                       labels = c(-1:(years-1)))
```

![](age_structured_schisto_model_files/figure-markdown_github/sim_mda_sac-1.png)