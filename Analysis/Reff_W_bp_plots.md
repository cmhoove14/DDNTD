*R*<sub>*e**f**f*</sub>*a**n**d**W*<sub>*b**p*</sub>*e**s**t**i**m**a**t**i**o**n*
================

``` r
W = 60
Lambda = 0.025
fit_pars <- base_pars
fit_pars["cvrg"] = 0

  M = 0.5*W*base_pars["H"]*phi_Wk(W, k_from_log_W(W))*rho_Wk(W, base_pars["zeta"], k_from_log_W(W))*base_pars["omega"]*base_pars["U"]*base_pars["m"]*base_pars["v"]  
  
  N_eq = Lambda_get_N_eq(Lambda, base_pars["K"], base_pars["mu_N"], base_pars["r"], base_pars["sigma"])
  I_eq = (base_pars["sigma"]*N_eq)/(base_pars["mu_I"]*(base_pars["mu_N"]+base_pars["sigma"])/Lambda+base_pars["mu_I"]+base_pars["sigma"])
  
  fit_pars["beta"] = Lambda*M/N_eq
  
  fit_pars["alpha"] = (W*(base_pars["mu_H"] + base_pars["mu_W"]))/(base_pars["omega"]*I_eq*base_pars["theta"])

W_seq <- exp_seq(1e-4, 200, 200)

Reff_R0 <- data.frame(W = W_seq) %>% 
  mutate(dWdt = sapply(W, dWdt_W, pars = fit_pars),
         Reff = sapply(W, R_eff_W_I_P, pars = fit_pars, I_P = 0.12,
                       PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1)[1,],
         R0 = sapply(W, R_eff_W_I_P, pars = fit_pars, I_P = 0.12,
                       PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1)[2,],
         Reff_W = sapply(W, Reff_W, pars = fit_pars)[1,],
         I_P = sapply(W, Reff_W, pars = fit_pars)[2,],
         Lambda = sapply(W, Reff_W, pars = fit_pars)[3,])

Reff_R0 %>% 
  ggplot(aes(x = W, y = Reff_W)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))  +
    scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       limits = c(1e-4, 200)) +
    labs(x = "Mean Worm Burden (W)",
         y = expression(R[eff]),
         title = expression(The~Effective~Reproductive~Rate))
```

![](Reff_W_bp_plots_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
Reff_R0 %>% 
  gather("R", "Value", Reff:Reff_W) %>% 
  ggplot(aes(x = W, y = Value, col = R)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() +
        scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       limits = c(1e-4, 200)) +
    labs(x = "Mean Worm Burden (W)",
         y = "Reproduction Number",
         title = expression(R[eff]~R[0]~Comparison))
```

![](Reff_W_bp_plots_files/figure-markdown_github/unnamed-chunk-1-2.png)

Makes sense that these don't match up exactly since the infected snail prevalence is constant in the *I*<sub>*P*</sub> derived *R*<sub>*e**f**f*</sub>

``` r
Reff_R0 %>% 
  ggplot(aes(x = I_P, y = Reff)) + geom_line() + theme_classic() +
    scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1),
                       labels = c("0.0001", "0.001", "0.01", "0.1"),
                       limits = c(1e-4, 0.25)) +
  labs(title = "Reff across infected snail prevalence")
```

    ## Warning: Removed 92 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/Reff_IP_comp-1.png)

``` r
Reff_R0 %>% 
  ggplot(aes(x = Lambda, y = Reff)) + geom_line() + theme_classic()# +
```

![](Reff_W_bp_plots_files/figure-markdown_github/Reff_Lambda_comp-1.png)

``` r
    scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1),
                       labels = c("0.0001", "0.001", "0.01", "0.1"),
                       limits = c(1e-4, 0.2))
```

    ## <ScaleContinuousPosition>
    ##  Range:  
    ##  Limits: -9.21 -- -1.61

``` r
W_get_N <- function(W, pars){
    ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

  # Get miracidial density as function of worm burden
    #Update clumping parameter
      kap = k_from_log_W(W)

    #Estimate mating probability
      phi = phi_Wk(W = W, k = kap)

    #Estimate density dependent fecundity
      rho = rho_Wk(W, zeta, k = kap)

    # Estimate total miracidia entering snail habitat
      M_tot = 0.5*H*omega*v*m*W*phi*rho*U

  # Get man-to-snail FOI and snail population sizeas solution given M_tot and other parameters
      #Estimate N
        N <- uniroot.all(function(N) K*(1-(mu_N+(beta*M_tot/N))/(r*(1+(beta*M_tot/(N*(mu_N+sigma))))))-N,
                         interval = c(0,K))
        
    I_P <- (N*sigma)/((N*mu_I*(mu_N+sigma)/beta*M_tot) + mu_I + sigma)    
        
    return(c(N, I_P, M_tot) )   
}

N_from_W <- data.frame(W = W_seq) %>% 
  mutate(N_W = sapply(W, W_get_N, pars = fit_pars)[1,],
         I_P = sapply(W, W_get_N, pars = fit_pars)[2,],
         M = sapply(W, W_get_N, pars = fit_pars)[3,])

N_from_W %>% 
  ggplot(aes(x = W, y = N_W)) + geom_line() + theme_classic()
```

![](Reff_W_bp_plots_files/figure-markdown_github/N_W_plot-1.png)

``` r
N_from_W %>% 
  ggplot(aes(x = W, y = I_P)) + geom_line() + theme_classic() +
          scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       limits = c(1e-4, 200)) 
```

![](Reff_W_bp_plots_files/figure-markdown_github/N_W_plot-2.png)

``` r
N_from_W %>% 
  ggplot(aes(x = W, y = M)) + geom_line() + theme_classic() +
    ylim(c(0,10)) +
          scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       limits = c(1e-4, 200)) 
```

    ## Warning: Removed 92 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/N_W_plot-3.png)

``` r
main_plot <- Reff_R0 %>% 
  ggplot(aes(x = W, y = dWdt*365)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 0, lty = 2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       limits = c(1e-4, 200)) +
    labs(x = "Mean Worm Burden (W)",
         y = expression(italic(frac(dW, dt))))

main_plot
```

![](Reff_W_bp_plots_files/figure-markdown_github/dwdt_plot-1.png)

``` r
inset <- Reff_R0 %>% 
  ggplot(aes(x = W, y = dWdt*365)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 0, lty = 2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_blank(),
          axis.title.y = element_blank()) +
    ylim(c(-1e-3,1e-3)) +
    scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                       limits = c(1e-4, 1))

inset
```

    ## Warning: Removed 108 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/dwdt_plot-2.png)

``` r
dd_product <- data.frame(W = W_seq) %>% 
  mutate(phi_rho_gam = sapply(W, 
                              function(W) phi_Wk(W, k_from_log_W(W)) * rho_Wk(W, age_strat_pars["zeta"], k_from_log_W(W)) * gam_Wxi(W, age_strat_pars["xi"])),
         phi_rho_z10_gam = sapply(W, 
                                  function(W) phi_Wk(W, k_from_log_W(W)) * rho_Wk(W, age_strat_pars["zeta"]*10, k_from_log_W(W)) * gam_Wxi(W, age_strat_pars["xi"])),
         M = sapply(W, 
                    function(W) phi_Wk(W, k_from_log_W(W)) * rho_Wk(W, age_strat_pars["zeta"], k_from_log_W(W)) * W*fit_pars["H"]*fit_pars["m"]*fit_pars["v"]*fit_pars["omega"]*fit_pars["U"]*0.5))

  dd_product %>% 
    ggplot(aes(x = W, y = phi_rho_gam)) + 
      geom_line() +
      geom_line(aes(x = W, y = phi_rho_z10_gam), col = 2) +
      theme_classic() +
      labs(y = expression(Phi~rho~gamma))
```

![](Reff_W_bp_plots_files/figure-markdown_github/dd_plots-1.png)

``` r
  dd_product %>% 
    ggplot(aes(x = W, y = M)) + 
      geom_line() +
      theme_classic() +
      labs(y = "M") +
    ylim(c(0,10)) +
      scale_x_continuous(trans = "log",
                         breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1),
                         labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                         limits = c(1e-4, 1))
```

    ## Warning: Removed 92 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/dd_plots-2.png)

``` r
Reff_IP <- data.frame(W = W_seq) %>% 
  mutate(Reff_I01 = sapply(W, R_eff_W_I_P, pars = fit_pars, I_P = 0.1,
                       PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1)[1,],
         Reff_I005 = sapply(W, R_eff_W_I_P, pars = fit_pars, I_P = 0.05,
                       PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1)[1,],
         Reff_I001 = sapply(W, R_eff_W_I_P, pars = fit_pars, I_P = 0.01,
                       PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1)[1,])

Reff_IP %>% 
  gather("I Prev", "Reff", Reff_I01:Reff_I001) %>% 
  ggplot(aes(x = W, y = Reff, col = `I Prev`)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 1, lty = 2) +
    theme_classic() +
        scale_x_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100"),
                       limits = c(1e-4, 200)) +
    labs(x = "Mean Worm Burden (W)",
         y = expression(R[eff]),
         title = expression(R[eff]~I[P]~Relationship))
```

![](Reff_W_bp_plots_files/figure-markdown_github/Reff_I_P-1.png)

``` r
W_bp_IP <- data.frame(I_P = seq(0.001,0.1,0.001)) %>% 
  mutate(Wbp = sapply(I_P, W_bp_I_P, pars = fit_pars,
                      PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1))

W_bp_IP %>% 
  mutate(Wbp = if_else(Wbp > 100, NA_real_, Wbp)) %>% 
  ggplot(aes(x = I_P, y = Wbp)) +
    geom_line(size = 1.2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))  +
    scale_y_continuous(trans = "log",
                       breaks = c(1e-4, 1e-3, 1e-2, 0.1, 1, 10),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10"),
                       limits = c(1e-4, 10)) +
    scale_x_continuous(breaks = c(0,0.025,0.05,0.075,0.1),
                       labels = c("0", "2.5", "5.0", "7.5", "10")) +    
    labs(y = expression(Worm~Burden~Breakpoint~(W[bp])),
         x = expression(paste("Infected Snail Prevalence (", I[P], ", %)")),
         title = "Breakpoint as function of snail prevalence")
```

![](Reff_W_bp_plots_files/figure-markdown_github/W_bp_I_P-1.png)

``` r
W_bp_R0 <- data.frame(R0 = seq(0.95,5,0.05)) %>% 
  mutate(Wbp = sapply(R0, W_bp_R_0, 
                      PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1,
                      zeta = age_strat_pars["zeta"], xi = age_strat_pars["xi"]))

W_bp_R0 %>% 
  ggplot(aes(x = R0, y = Wbp)) +
    geom_line(size = 1.2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))  +
    scale_x_continuous(breaks = c(1:5),
                       limits = c(1,5)) +    
    scale_y_continuous(breaks = seq(0.2, 1.6, by = 0.2),
                       limits = c(0.2,1.6)) +    
    labs(y = expression(italic(W[bp])),
         x = expression(italic(R[0])))
```

    ## Warning: Removed 1 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/Wbp_R0-1.png)

``` r
W_bp_R0 <- W_bp_R0 %>% 
  mutate(bp_cvrg = sapply(Wbp, cvrg_from_W_bp, W_t = W, epsilon = 0.94))

W_bp_R0 %>% 
  ggplot(aes(x = R0, y = bp_cvrg)) +
    geom_line(size = 1.2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))  +
    scale_x_continuous(breaks = c(1:4),
                       limits = c(1,4)) +    
    labs(y = expression(MDA~Coverage~to~reach~breakpoint),
         x = expression(R[0]),
         title = "Coverage to reach breakpoint as function of basic reproductive rate")
```

    ## Warning: Removed 21 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/cvrg-1.png)

``` r
test_omegas <- seq(0.5,0.99,0.01)
result_Wbps <- numeric(length = length(test_omegas))

for(i in test_omegas){
  fit_pars -> use_pars ; use_pars["omega"] = fit_pars["omega"]*i
  
  result_Wbps[which(test_omegas == i)] <- W_bp(use_pars, PDD = phi_Wk, DDF=rho_Wk, DDI = nil_1)
}
```

    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"
    ## [1] "R_0 less than 1, no breakpoint exists"

``` r
W_bp_omega <- data_frame(omega_red = 100 - test_omegas*100,
                         W_bp = result_Wbps)
```

    ## Warning: `data_frame()` is deprecated, use `tibble()`.
    ## This warning is displayed once per session.

``` r
W_bp_omega %>% 
  ggplot(aes(x = omega_red, y = W_bp)) +
    geom_line(size = 1.2) +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))  +
    scale_x_continuous(breaks = c(0, 10, 20,30),
                       limits = c(0, 31)) +    
    labs(y = expression(italic(W[bp])),
         x = expression(paste("% reduction in exposure/contamination ", (italic(omega)))))
```

    ## Warning: Removed 21 rows containing missing values (geom_path).

![](Reff_W_bp_plots_files/figure-markdown_github/W_bp_omega-1.png)

``` r
W_bp_omega <- W_bp_omega %>% 
  mutate(bp_cvrg = sapply(W_bp, cvrg_from_W_bp, W_t = W, epsilon = 0.94))

W_bp_omega %>% 
  ggplot(aes(x = omega_red, y = bp_cvrg*100)) +
    geom_line(size = 1.2) +
    theme_classic() +
    geom_hline(yintercept = 100, lty = 2) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))  +
    scale_x_continuous(breaks = c(0, 10, 20,30),
                       limits = c(0, 31)) +    
    labs(y = expression(MDA~Coverage~to~reach~W[bp]),
         x = expression(paste("% reduction in exposure/contamination ", (italic(omega))))) +
    annotate("rect", xmin = 0, xmax = 42, ymin = 100.01, ymax = Inf, 
             alpha = 0.2, col = "red", fill = "red")
```

    ## Warning: Removed 21 rows containing missing values (geom_path).

    ## Warning: Removed 1 rows containing missing values (geom_rect).

![](Reff_W_bp_plots_files/figure-markdown_github/cvrg_omega-1.png)

``` r
W_bp_omega_Wt <- expand.grid(W_bp = W_bp_omega$W_bp,
                             W_t = exp_seq(1,100,30)) %>% 
  mutate(omega_red = rep(W_bp_omega$omega_red, times = 30),
         bp_cvrg = mapply(cvrg_from_W_bp, W_t, W_bp, MoreArgs = list(epsilon = 0.94))*100,
         bp_cvrg = case_when(bp_cvrg > 0 & bp_cvrg <100 ~ bp_cvrg,
                             bp_cvrg < 0 ~ NA_real_,
                             bp_cvrg > 100 ~ NA_real_))


W_bp_omega_Wt %>% 
  ggplot(aes(y = W_t, x = omega_red, z = bp_cvrg)) + 
    geom_tile(aes(fill = bp_cvrg)) + 
    scale_x_continuous(breaks = c(0, 10, 20,30),
                       limits = c(0, 31)) +
    scale_y_continuous(trans = "log",
                       breaks = c(1,10,100)) +
    scale_fill_viridis(option = "magma", limits = c(0,100), direction = -1) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) +
    labs(y = expression(paste("Pre-MDA mean worm burden, ", W[t])),
         x = expression(paste("% reduction in exposure/contamination ", (italic(omega)))),
         fill = expression(italic(cvrg))) #+ annotate("text", x = 15, y = 18, label = expression(paste(W[bp], " not achievable \nwith MDA")), col = "red", size = 6)
```

    ## Warning: Removed 570 rows containing missing values (geom_tile).

![](Reff_W_bp_plots_files/figure-markdown_github/unnamed-chunk-2-1.png)
