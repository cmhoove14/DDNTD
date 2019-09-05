Fitting W and kappa to egg burden and prevalence data, implications for estimation of mating probability
================

Clumping parameter as function of infection intensity
-----------------------------------------------------

S. haematobium infection data consisting of egg burden and prevalence estimates was collected from a literature review. We restrict estimation to those studies which report an arithmetic mean egg burden and prevalence estimate

The mean worm burden, *W*, and clumping parameter, *κ*, are estimated as functions of the mean egg burden, ℰ, and prevalence of **signs of infection**, which can be interpreted as the prevalence of mated worm pairs, *Ω*:

![](http://latex.codecogs.com/gif.latex?%5Cmathcal%7BE%7D%3D0.5W%5Cphi(W,%5Ckappa)%5Crho(W,%5Ckappa,%5Czeta)m)

![](http://latex.codecogs.com/gif.latex?%5COmega%3D1-2%5Cbig(1+%5Cfrac%7BW%7D%7B2%5Ckappa%7D%5Cbig)%5E%7B-k%7D+%5Cbig(1+%5Cfrac%7BW%7D%7B%5Ckappa%7D%5Cbig)%5E%7B-k%7D)

### Naive estimation

Estimate aggregation parameter from egg intensity and prevalence assuming mean worm burden is equivalent to egg burden times 2 (for assumed 1:1 sex ratio) divided by estimate of egg output per mated female worm (assuming no DD fecundity and no mating function)
![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/kap_from_intensity-1.png)

### Estimation incorporating mating probability

Uses rootsolving to solve for mean worm burden and aggregation from equations for egg output per mated female worm and prevalence of mated pairs

    ## diagonal element is zero 
    ## [1] 2

![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/fit_kap_w-1.png)

### Estimation incorporating mating probability and fecundity reduction due to crowding

    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 1
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2

    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 1
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2

    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 1
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2

    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2

    ## Error in integrate(mate_integral, 0, 2 * pi): the integral is probably divergent

    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2
    ## diagonal element is zero 
    ## [1] 2

So this isn't working for some reason... Let's look at estimates of W and

### Comp results

![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/comp_naive_fit-1.png)

In general, it appears that the clumping parameter is most affected by including density dependence in the estimation. Kappa tends to increase with this method, as does W to some extent. Also appears W increases somewhat proportionally with its own magnitude (i.e. bigger Ws means bigger differences in W between naive and rootsolve estimations).

### kappa as a non-linear (saturating) function of W

![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/k_w_sat_fit-1.png)

### kappa as a linear function of W

![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/k_w_lin_fit-1.png)

### kappa as an exponential function of W

![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/k_w_exp_fit-1.png)

Implications for mating probability
-----------------------------------

![](Fitting_W_and_kappa_to_data_implications_for_mating_prob_files/figure-markdown_github/mate_prob_kappa-1.png)

So if the clumping parameter approaches 0 as the worm burden approaches 0 (as with the saturating function, black line above, or for a linear function with fixed intercept at *κ*=0) the mating probability stays quite high even as the worm burden decreases. This could have serious implications for the breakpoint which could only exist at extremely small worm burdens as a result.
