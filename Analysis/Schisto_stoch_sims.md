Schisto stochastic mode Probability of Elimination
================
Chris Hoover
July 11, 2019

``` r
years <- 10
mda.eff <- 0.93
snail.eff <- 0.85

#Base time and starting values for state variables
base_start <- c(S=5000, E=0, I=0, Wt=10, Wu=10)
base_time <- c(0:(365*(years*2)))


#Intervention parameter sets 
  # 80 % MDA coverage
    base_pars -> pars_cvrg80 ; pars_cvrg80["cvrg"] = .80
  
  # 80% MDA coverage and 20% reduction in snail habitat
    pars_cvrg80 -> pars_cvrg80_C20 ; pars_cvrg80_C20["C"] = pars_cvrg80["C"]*0.8

  # 60 % MDA coverage
    base_pars -> pars_cvrg60 ; pars_cvrg60["cvrg"] = .60
  
  # 60% MDA coverage and 20% reduction in snail habitat
    pars_cvrg60 -> pars_cvrg60_C20 ; pars_cvrg60_C20["C"] = pars_cvrg60["C"]*0.8

#Events data frames    
mda.events <- data.frame(var = rep('Wt', years),
                         time = c(1:years)*365,
                         value = rep((1 - mda.eff), years),
                         method = rep('mult', years))

#Snail control dataframes
snail.annual.events <- data.frame(var = rep(c('S', 'E', 'I'), years),
                                  time = rep(c(1:years)*365+1, each = 3),
                                  value = rep((1 - snail.eff), years*3),
                                  method = rep('mult', years*3))

#MDA and snail control events together
mda.snail.events <- rbind(mda.events, snail.annual.events) %>% 
  arrange(time)

#Run to equibrium with base parameter set
base_eqbm <- runsteady(y = base_start, func = schisto_base_mod,
                       parms = pars_cvrg80)[["y"]]

#Run to equibrium with C20_red parameter set
C20red_eqbm <- runsteady(y = base_start, func = schisto_base_mod,
                       parms = pars_cvrg80_C20)[["y"]]
```

``` r
#simulate annual MDA 
schisto_sims <- bind_rows(sim_schisto_mod(nstart = base_eqbm, 
                                          time = base_time, 
                                          model = schisto_base_mod,
                                          parameters = pars_cvrg80,
                                          events_df = NA) %>% 
                            mutate(W = Wt*pars_cvrg80["cvrg"] + Wu*(1-pars_cvrg80["cvrg"]),
                                   Sim = "Nothing"),
                          sim_schisto_mod(nstart = base_eqbm, 
                                          time = base_time, 
                                          model = schisto_base_mod,
                                          parameters = pars_cvrg80,
                                          events_df = mda.events) %>% 
                            mutate(W = Wt*pars_cvrg80["cvrg"] + Wu*(1-pars_cvrg80["cvrg"]),
                                   Sim = "Annual MDA 80% coverage"),
                          sim_schisto_mod(nstart = base_eqbm, 
                                          time = base_time, 
                                          model = schisto_base_mod,
                                          parameters = pars_cvrg80,
                                          events_df = snail.annual.events) %>% 
                            mutate(W = Wt*pars_cvrg80["cvrg"] + Wu*(1-pars_cvrg80["cvrg"]),
                                   Sim = "Annual Snail Control"),
                          sim_schisto_mod(nstart = base_eqbm, 
                                          time = base_time, 
                                          model = schisto_base_mod,
                                          parameters = pars_cvrg60,
                                          events_df = mda.snail.events) %>% 
                            mutate(W = Wt*pars_cvrg60["cvrg"] + Wu*(1-pars_cvrg60["cvrg"]),
                                   Sim = "Annual MDA 60% + Snail Control"),
                          sim_schisto_mod(nstart = base_eqbm, 
                                          time = base_time, 
                                          model = schisto_base_mod,
                                          parameters = pars_cvrg60_C20,
                                          events_df = mda.events) %>% 
                            mutate(W = Wt*pars_cvrg60_C20["cvrg"] + Wu*(1-pars_cvrg60_C20["cvrg"]),
                                   Sim = "Annual MDA 60% + Habitat reduction"),
                          sim_schisto_mod(nstart = base_eqbm, 
                                          time = base_time, 
                                          model = schisto_base_mod,
                                          parameters = pars_cvrg60_C20,
                                          events_df = mda.snail.events) %>% 
                            mutate(W = Wt*pars_cvrg60_C20["cvrg"] + Wu*(1-pars_cvrg60_C20["cvrg"]),
                                   Sim = "Annual MDA 60% + Habitat reduction + Snail Control"))

schisto_sims %>% 
    ggplot(aes(x = time, y = W, col = Sim)) +
    annotate("rect", xmin = 365, xmax = max(mda.events$time), ymin = -Inf, ymax = Inf,
             alpha = .2) +
      geom_line(size = 1.2) +
      scale_x_continuous(breaks = c(0:(years*2))*365,
                         labels = c(-1:(years*2-1))) +
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(title = "Human infection dynamics under different interventions",
           x = "Time (Years)",
           y = expression(Mean~Worm~Burden~italic((W))))
```

![](Schisto_stoch_sims_files/figure-markdown_github/deterministic_model_sims-1.png)

``` r
set.seed(10)

schisto_stoch_sims <- bind_rows(lapply(c(1:10), sim_schisto_stoch, 
                                       pars = pars_cvrg80, events = mda.events))

schisto_stoch_sims %>% 
  group_by(sim) %>% 
  summarise(I_fin = I[which(time == max(time))],
            Wt_fin = Wt[which(time == max(time))],
            Wu_fin = Wu[which(time == max(time))],
            elim = if_else(I_fin + Wt_fin + Wu_fin == 0, 1, 0))
```

    ## # A tibble: 10 x 5
    ##      sim I_fin Wt_fin Wu_fin  elim
    ##    <int> <dbl>  <dbl>  <dbl> <dbl>
    ##  1     1   122     22     20     0
    ##  2     2   135     16     21     0
    ##  3     3   140     21     22     0
    ##  4     4    75     10      8     0
    ##  5     5   138     21     16     0
    ##  6     6    93     14     11     0
    ##  7     7   127     21     24     0
    ##  8     8    26      4      2     0
    ##  9     9   170     37     40     0
    ## 10    10   141     23     25     0

``` r
schisto_stoch_sims_plot <- schisto_stoch_sims %>% 
  ggplot(aes(x = time, y = W, col = as.factor(sim))) +
    annotate("rect", xmin = 365, xmax = max(mda.events$time), ymin = -Inf, ymax = Inf,
             alpha = .2) +
    geom_line(size = 1.1) +
    scale_x_continuous(breaks = c(0:(years*2))*365,
                       labels = c(-1:((years*2)-1))) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "time (years)",
         y = expression(mean~worm~burden~(italic(W))),
         title = "Human infection dynamics from stochastic model", 
         subtitle = paste0("Anual MDA for ", years/2, " years at ", pars_cvrg80["cvrg"]*100, "% coverage"))

schisto_stoch_sims_plot
```

![](Schisto_stoch_sims_files/figure-markdown_github/stoch_sims_plot_example-1.png)

``` r
n_sims <- 100
n_reps <- 20

set.seed(430)
# 80% MDA coverage only
pe_mda80 <- replicate(n_reps, sum(replicate(n_sims, sim_schisto_stoch_elim_1_0(pars_cvrg80, 
                                                                               base_eqbm, 
                                                                               mda.events))))
mean(pe_mda80)
```

    ## [1] 7.8

``` r
sd(pe_mda80)
```

    ## [1] 2.117595

``` r
#60% MDA coverage only
pe_mda60 <- replicate(n_reps, sum(replicate(n_sims, sim_schisto_stoch_elim_1_0(pars_cvrg60, 
                                                                               base_eqbm, 
                                                                               mda.events))))
mean(pe_mda60)
```

    ## [1] 0

``` r
sd(pe_mda60)
```

    ## [1] 0

``` r
#60% MDA coverage and annual snail control
#pe_mda60_snails <- replicate(n_reps, sum(replicate(n_sims, sim_schisto_stoch_elim_1_0(pars_cvrg60, 
 #                                                                   base_eqbm, 
 #                                                                   mda.snail.events))))
#mean(pe_mda60_snails)
#sd(pe_mda60_snails)

#60% MDA coverage and 20% snail habitat reduction
pe_mda60_Cred <- replicate(n_reps, sum(replicate(n_sims, sim_schisto_stoch_elim_1_0(pars_cvrg60_C20, 
                                                                                    C20red_eqbm, 
                                                                                    mda.events))))
mean(pe_mda60_Cred)
```

    ## [1] 4.35

``` r
sd(pe_mda60)
```

    ## [1] 0

``` r
#60% MDA coverage and 20% snail habitat reduction and annual snail control
#pe_mda60_snails_Cred <- replicate(n_reps, sum(replicate(n_sims, sim_schisto_stoch_elim_1_0(pars_cvrg60_C20, 
#                                                                         base_eqbm, 
#                                                                         mda.snail.events))))
#mean(pe_mda60_snails_Cred)
#sd(pe_mda60_snails_Cred)
```

``` r
data.frame("Intervention" = c("Annual MDA 80% Coverage",
                              "Annual MDA 60% Coverage",
                              "Annual MDA 60% Coverage + 20% Habitat reduction"),
           "Mean P(e)" = round(c(mean(pe_mda80), mean(pe_mda60), mean(pe_mda60_Cred)), 2),
           "St Dev P(e)" = round(c(sd(pe_mda80), sd(pe_mda60), sd(pe_mda60_Cred)), 2)) %>% 
knitr::kable(format = "markdown")
```

| Intervention                                    |  Mean.P.e.|  St.Dev.P.e.|
|:------------------------------------------------|----------:|------------:|
| Annual MDA 80% Coverage                         |       7.80|         2.12|
| Annual MDA 60% Coverage                         |       0.00|         0.00|
| Annual MDA 60% Coverage + 20% Habitat reduction |       4.35|         2.18|
