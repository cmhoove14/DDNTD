Critical Transitions and Schistosomiasis control
================
October 11, 2019

Critical transitions general theory
-----------------------------------

Critical transitions thought of in terms of a general process
*d**x* = *f*(*x*, *θ*)*g*(*x*, *θ*)*d**W*
 where x is the state of the system, *f*(*x*, *θ*) describes the deterministic part of the system, and *g*(*x*, *θ*) describes how stochasticity interacts with the system. Changes in *θ* can move the system closer to a threshold where a transition is likely to occur. In a disease setting, particularly for parasitic diseases where transmission is a function of infection intensity, *f*(*x*, *θ*)≈*R*<sub>*e**f**f*</sub>, *θ* ≈ *R*<sub>0</sub>, *x* ≈ *W* (system state = infection intensity), and *g*(*x*, *θ*) is related to stochastic components of transmission such as heterogeneities in egg shedding, susceptibility, and unkown inputs from external sources.

For schistosomiasis and other helminths, *R*<sub>*e**f**f*</sub> can be measured as the product of *R*<sub>0</sub> and the magnitude of density dependencies at a particular system state, measured in terms of worm burden, *W*. Density dependencies common to helminths include positive density dependent mate limitation,*Φ*, and negative density dependent fecundity due to crowding, *ρ*. Both of these quantities can be estimated from data on the distribution of parasites among definitive hosts, often assumed to be negative binomially distributed and measured in terms of the mean worm burden, *W*, and its dispersion among the population, *κ*, thus we have:
*R*<sub>*e**f**f*</sub>(*W*)=*R*<sub>0</sub>*Φ*(*W*, *κ*)*ρ*(*W*, *κ*, *ζ*)
 Before any sort of intervention such as PZQ administration (e.g. perturbation to the system), we can assume the system is at its endemic equilibrium where *R*<sub>*e**f**f*</sub> = 1, thus we can estimate *R*<sub>0</sub> = (*Φ*(*W*, *κ*)*ρ*(*W*, *κ*, *ζ*))<sup>−1</sup>.

#### Model and connection to epidemiological data

We can connect this to a simplified helminth model as presented in \[1\] with two state variables: the prevalence of infection in the environmental reservoir, *y*, (e.g. snail infection prevalence for schistosomiasis) and the mean parasite burden in the definitive human host population, *W*:
$$\\frac{dW}{dt}=\\alpha y-\\mu\_WW$$

$$\\frac{dy}{dt}=\\beta 0.5W\\Phi(W,\\kappa)\\rho(W,\\kappa)(1-y)-\\mu\_yy$$
We can then use equilibrated mean worm burden, $W^\*=\\frac{\\alpha yN}{\\mu\_W}$, and $R\_0=\\frac{\\alpha\\beta}{2\\mu\_W\\mu\_y}$ to reduce the snail equation to:

$$\\frac{dy}{dt}=y\[R\_0\\Phi(W^\*,\\kappa)\\rho(W^\*,\\kappa)(1-y)-1\]$$

We can also express snail infection dynamics in terms of egg output from infected individuals, estimated as ℰ = 0.5*W**Φ*(*W*, *κ*)*ρ*(*W*, *κ*) to give:
$$\\frac{dy}{dt}=\\beta\\mathcal{E}(W,\\kappa)(1-y)-\\mu\_yy$$

Methods for detecting "early warning signals" of critical transitions abound and typically draw from the theory that the rate of return to an equilibirum is decreased as the system approaches a transition (i.e. perturb an unstable system far from the equilibrium and it takes a long time to return). This can be measured e.g. via increased autocorrelation in dense time series.

This also implies that systems that return rapidly to the pre-purturbation equilibrium are very stable and are far from a critical transition. Persistent hotspot phenomena in which community infection intensities return to pre-treatment levels even before the next treatment can be thought of as a highly stable system, far from a critical transition. This is why perturbations to the system (MDA) have little long term impact. Communities with lower bounce back rates are less stable and therefore take longer to return to their endemic (pre-treatment) equilibirum. This implies lower transmission potential, here embodied by *R*<sub>0</sub>.

Rebound can be measured by the bounce back rate, *B**B**R*, an empirical estimator of *R*<sub>*e**f**f*</sub>. It is expected to be high in areas with conditions suggestive of a high *R*<sub>0</sub>, indicative of a stable system far from a tipping point. This could be tested in an empirical statistical model with sufficient data pertaining to the determinants of *R*<sub>0</sub> through time and longitudinal measurement of *B**B**R*. For instance, Spear et al \[2\] express *R*<sub>0</sub> in terms of invariant biological components of schistosomiasis transmission, *P*<sub>*b*</sub>, and site-specific parameters related to the behavioral and environmental determinants of transmission, *P*<sub>*s*</sub>. Parameters included in *P*<sub>*s*</sub> relate to human water contact, uninfected snail density, contamination behaviors related to sanitation, the amount of snail habitat, the amount of surface water, and seasonally varying factors such as temperature, rainfall, vegetation indices, and seasonal water contact patterns.

References
----------

1. Gurarie D, King CH. Population biology of schistosoma mating, aggregation, and transmission breakpoints: More reliable model analysis for the end-game in communities at risk. Munderloh UG, editor. PLoS ONE. 2014;9: e115875. doi:[10.1371/journal.pone.0115875](https://doi.org/10.1371/journal.pone.0115875)

2. Spear R, Zhong B, Liang S. Low transmission to elimination: Rural development as a key determinant of the end-game dynamics of schistosoma japonicum in china. Tropical Medicine and Infectious Disease. 2017;2: 35. doi:[10.3390/tropicalmed2030035](https://doi.org/10.3390/tropicalmed2030035)
