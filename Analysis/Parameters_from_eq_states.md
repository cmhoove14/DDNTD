Fitting procedure
================

Man-to-snail FOI, *Λ*, from equilibrium solutions to snail infection dynamics and input infected snail prevalence
-----------------------------------------------------------------------------------------------------------------

Beginning with a basic schistosomiasis model with *S* − *E* − *I* infection dynamics and logistic population growth among the intermediate host snail population we have three ODEs and one simple relation between each infection class and the total snail population, N:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdS%7D%7Bdt%7D%3Dr%5CBig(1-%5Cfrac%7BN%7D%7BK%7D%5CBig)%5CBig(S+E%5CBig)-(%5Cmu_N+%5CLambda)%20S)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdE%7D%7Bdt%7D%3D%5CLambda%20S-(%5Cmu_N+%5Csigma)E)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdI%7D%7Bdt%7D%3D%5Csigma%20E%20-%20%5Cmu_I%20I)

![](http://latex.codecogs.com/gif.latex?N%3DS+E+I)

Where *Λ* is the man-to-snail force of infection (FOI). At equilibrium, ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdS%7D%7Bdt%7D%3D%5Cfrac%7BdE%7D%7Bdt%7D%3D%5Cfrac%7BdI%7D%7Bdt%7D%3D0) we have:

![](http://latex.codecogs.com/gif.latex?E%5E*%3D%5Cfrac%7B%5CLambda%20S%5E*%7D%7B%5Cmu_N+%5Csigma%7D)

![](http://latex.codecogs.com/gif.latex?I%5E*%3D%5Cfrac%7B%5Csigma%20E%5E*%7D%7B%5Cmu_I%7D%3D%5Cfrac%7B%5Csigma%5CLambda%20S%5E*%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D)

and therefore:

![](http://latex.codecogs.com/gif.latex?N%5E*%3DS%5E*%5CBig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5CLambda%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5CBig))

In addition, we can solve equation 1 representing susceptible snail dynamics for *N*<sup>\*</sup> in terms of *Λ* as:

![](http://latex.codecogs.com/gif.latex?N%5E*(%5CLambda)%3DK%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig))

Which gives:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig)%7D%7B%5CBig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5CLambda%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5CBig)%7D)

and

![](http://latex.codecogs.com/gif.latex?I%5E*%3D%5Cfrac%7BK%5Csigma%5CLambda%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig)%7D%7B%5CBig(%5Cmu_I(%5Cmu_N+%5Csigma)%5CBig)%5CBig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5CLambda%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5CBig)%7D)

Which simplifies to:

![](http://latex.codecogs.com/gif.latex?I%5E*%3D%5Cfrac%7BK%5Csigma%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig)%7D%7B%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda%7D+%5Cmu_I+%5Csigma%5CBig)%7D%3D%5Cfrac%7B%5Csigma%20N%5E*%7D%7B%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda%7D+%5Cmu_I+%5Csigma%5CBig)%7D)

Therefore with infected snail prevalence, *I*<sub>*P*</sub> = *I*<sup>\*</sup>/*N*<sup>\*</sup>, and other snail population and infection parameters as inputs, we can estimate *Λ* as:

![](http://latex.codecogs.com/gif.latex?I_P%3D%5Cfrac%7B%5Csigma%7D%7B%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda%7D+%5Cmu_I+%5Csigma%7D)

![](http://latex.codecogs.com/gif.latex?%5CLambda%3D%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D)

Contamination fraction, *Ω*, and snail FOI parameters, (*β* or *Λ*<sub>0</sub>), as function of equilibrium worm burden in children and adults
----------------------------------------------------------------------------------------------------------------------------------------------

We first assume that differences in child and adult infection rates arise predominately from variation in behavior that affects both exposure and contribution to infection in the same manner. Parameter *Ω* therefore represents the relative exposure/contamination of children to adults, *Ω* = *ω*<sub>*C*</sub>/*ω*<sub>*A*</sub>, which can be estimated from the equilibrium infection ratio *Ω* = *W*<sub>*C*</sub><sup>\*</sup>/*W*<sub>*A*</sub><sup>\*</sup>.

With this parameter, we can estimate the equilibirum miracidial density *M* as:

![](http://latex.codecogs.com/gif.latex?M%3D0.5%5Cmathbf%7BH%7Dm%5Comega_A(W_C%5E*%20h_C%5CPhi(W_C%5E*)%5Crho(W_C%5E*)U_C%5COmega+W_A%5E*%20h_A%5CPhi(W_A%5E*)%5Crho(W_A%5E*)U_A))

### Saturating man-to-snail FOI

Now, since *Λ* has been estimated above, *N*<sup>\*</sup> can be estimated as a function of *Λ*, and we have an equilibrium estimate of miracidial density, *M*, we can estimate the baseline miracidial invasion rate, *Λ*<sub>0</sub>, of the saturating form of the man-to-snail FOI:

![](http://latex.codecogs.com/gif.latex?%5CLambda_0%3D%5Cfrac%7B%5CLambda%7D%7B(1-e%5E%7B-M/N%5E*(%5CLambda)%7D)%7D)

Except in rare circumstances, *Λ* ≈ *Λ*<sub>0</sub> because of the high *M*/*N* ratio in transmission settings that have not been intervened upon (i.e. equilibrium settings).

### Linear man-to-snail FOI

As an alternative, we can use the more common (but perhaps less accurate) linear form of the man-to-snail FOI:

![](http://latex.codecogs.com/gif.latex?%5CLambda%3D%5Cfrac%7B%5Cbeta%20M%7D%7BN%5E*(%5CLambda)%7D)

to estimate the probability of snail infection given miracidial density, *β* as:

![](http://latex.codecogs.com/gif.latex?%5Cbeta%3D%5Cfrac%7B%5CLambda%20N%5E*(%5CLambda)%7D%7BM%7D)

Worm establishment rate, *α*, as function of equilibrium worm burden
--------------------------------------------------------------------

Assuming the observed mean worm burden in each population fraction prior to intervention is approximately at equilibrium, the worm acquisition rate, *λ*<sub>*i**j*</sub>, for each group can be estimated as a function of the exposure fractions, *ω*<sub>*a*</sub> and *ω*<sub>*c*</sub>, the infected snail density, *I*<sup>\*</sup> = *I*<sub>*P*</sub>*N*<sup>\*</sup>, the snail cercarial shedding rate, *θ*, and the probability of cercarial establishment per exposure, *α*, with *α* being the only unknown and estimated as:

![](http://latex.codecogs.com/gif.latex?%5Calpha%3D%5Cfrac%7BW_i%5E*(%5Cmu_W+%5Cmu_H_i)%7D%7B%5Comega_i%20I%5E*%5Ctheta%7D)

*R*<sub>*e**f**f*</sub> derivation with linear snail FOI
--------------------------------------------------------

The effective reproduction number, *R*<sub>*e**f**f*</sub>, for schistosomiasis is defined as the number of mated adult female worms produced by a single adult female worm over the course of her lifetime. Unlike teh basic reproduction number, *R*<sub>0</sub>, *R*<sub>*e**f**f*</sub> changes as a function of the current level of infection, here measured in terms of the mean worm burden, *W*. From Anderson and May **Infectious Diseases of Humans**, we can estimate the effective reproduction number from the rate of worm burden change as:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5Clambda-(%5Cmu_W+%5Cmu_H)W)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D(%5Cmu_W+%5Cmu_H)W%5CBig(%5Cfrac%7B%5Clambda%7D%7B(%5Cmu_W+%5Cmu_H)W%7D-1%5CBig))

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D(%5Cmu_W+%5Cmu_H)W%5CBig(R_%7Beff%7D(W)-1%5CBig))

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W)%3D%5Cfrac%7B%5Clambda%7D%7B(%5Cmu_W+%5Cmu_H)W%7D)

Where *λ* is the worm acquisition rate or snail-to-man FOI and is a function of the exposure coefficient, *ω*, probability of worm acquisition given exposure, *α*, and cercarial concentration estimated as a function of cercarial shedding and the infected snail population assuming fast snail infection dynamics: ![](http://latex.codecogs.com/gif.latex?C%3D%5Ctheta%20I%5E*). If we incorporate linear snail FOI and the infected snail population in terms of the total snail population, *N*, we have:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W)%3D%5Cfrac%7B%5Calpha%5Comega%5Ctheta%7D%7BW%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%7D%5Ctimes%5Cfrac%7B%5Csigma%20N%7D%7B%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda%7D+%5Cmu_I+%5Csigma%5CBig)%7D)

Given that N can be estimated as a function of *Λ* which is a function of *W*, we can estimate *N* and ![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D) from input W by solving for each from:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W)%3D%5Cfrac%7B%5Calpha%5Comega%5Ctheta%5Csigma%20N(W)%7D%7BW%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda(W)%7D+%5Cmu_I+%5Csigma%5CBig)%7D)

![](http://latex.codecogs.com/gif.latex?N(W)%3DK%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda(W)%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda(W)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig))

Or:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W)%3D%5Cfrac%7BK%5Calpha%5Comega%5Ctheta%5Csigma%5Cbig(1+%5Cfrac%7B%5CLambda(W)%7D%7B%5Cmu_N+%5Csigma%7D-%5Cfrac%7B%5Cmu_I%7D%7Br%7D-%5Cfrac%7B%5CLambda(W)%7D%7Br%7D%5Cbig)%7D%7BW%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda(W)%7D+%5Cmu_I+%5Csigma%5CBig)%5CBig(1+%5Cfrac%7B%5CLambda(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%7D)

Alternatively, we can make one more substitution based on the fact that infected snail prevalence can be estimated as a function of snail FOI:

![](http://latex.codecogs.com/gif.latex?I_P%3D%5Cfrac%7B%5Csigma%7D%7B%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda%7D+%5Cmu_I+%5Csigma%7D)

therefore we can substitute in terms of infected snail prevalence, *I*<sub>*P*</sub>, to get:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W,%20I_P)%3D%5Cfrac%7B%5Calpha%5Comega%5Ctheta%20N%20I_P%7D%7BW%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%7D)

And then with linear FOI, ![](http://latex.codecogs.com/gif.latex?%5CLambda%3D%5Cbeta%20M/N)substitute for N in terms of *I*<sub>*P*</sub>:

![](http://latex.codecogs.com/gif.latex?N%3D%5Cfrac%7B%5Cbeta%20M%5Cbig(%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%5Cbig)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%3D%5Cfrac%7B%5Cbeta%5Comega%20WHmv%5CPhi(W)%5Crho(W)U%5Cbig(%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%5Cbig)%7D%7B2%5Cmu_I(%5Cmu_N+%5Csigma)%7D)

and get:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W,%20I_P)%3D%5Cfrac%7B%5Calpha%5Comega%5E2%5Ctheta%5Cbeta%5Cphi(W)%5Crho(W)HmvU(%5Csigma-I_P%5Cmu_I-I_P%5Csigma)%7D%7B2%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5Cbig(%5Cmu_I(%5Cmu_N+%5Csigma)%5Cbig)%7D)

### *R*<sub>*e**f**f*</sub> for stratified worm burdens

While the mean worm burden is often modeled as a single state variable, it is more accurate to further stratify the compartment to more accurately capture the skewed distribution of worms amongst the human population or to more accurately model interventions such as MDA that are often only administered to SAC or that systematically miss particular portions of the population. We therefore wish to express ![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D) as a function of the population mean worm burden, ![](http://latex.codecogs.com/gif.latex?%5Cbar%7BW%7D), derived from multiple treatment groups, *i*, and age groups, *j*, where:

![](http://latex.codecogs.com/gif.latex?%5Cbar%7BW%7D%3D%5Csum_i%5Csum_j%20h_%7Bij%7DW_%7Bij%7D)

The net ![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D) can be thought of similarly as some weighted average of contributions from all treatment groups:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(%5Cbar%7BW%7D)%3D%5Csum_i%5Csum_jh_%7Bij%7D%5Cfrac%7B%5Calpha%5Comega_i%5Ctheta%5Csigma%20N(%5Cbar%7BW%7D)%7D%7BW_%7Bij%7D%5Cbig(%5Cmu_W+%5Cmu_%7BH_i%7D%5Cbig)%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda(%5Cbar%7BW%7D)%7D+%5Cmu_I+%5Csigma%5CBig)%7D)

![](http://latex.codecogs.com/gif.latex?N(%5Cbar%7BW%7D)%3DK%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda(%5Cbar%7BW%7D)%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda(%5Cbar%7BW%7D)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig))

Which we can express as:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(%5Cbar%7BW%7D)%3D%5Cfrac%7B%5Calpha%5Ctheta%5Csigma%20N(%5Cbar%7BW%7D)%7D%7B%5Cbar%7BW%7D%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda(%5Cbar%7BW%7D)%7D+%5Cmu_I+%5Csigma%5CBig)%7D%5Csum_i%20h_%7Bi%7D%5Cfrac%7B%5Comega_i%7D%7B%5Cbig(%5Cmu_W+%5Cmu_%7BH_i%7D%5Cbig)%7D)

![](http://latex.codecogs.com/gif.latex?N(%5Cbar%7BW%7D)%3DK%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda(%5Cbar%7BW%7D)%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda(%5Cbar%7BW%7D)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig))

Worm burden breakpoint estimation
---------------------------------

The breakpoint worm burden population size, ![](http://latex.codecogs.com/gif.latex?W_%7Bbp%7D) is an unstable equilibrium at which each female worm replaces herself, i.e. ![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D%3D1) therefore:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D(W_%7Bbp%7D,%20I_P)%3D1%3D%5Cfrac%7B%5Calpha%5Comega%5E2%5Ctheta%5Cbeta%5Cphi(W_%7Bbp%7D)%5Crho(W_%7Bbp%7D)HmvU%20I_P%7D%7B2%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5Cbig(%5Cmu_I(%5Cmu_N+%5Csigma)%5Cbig)%7D)

From this, it's clear that interplay between the density dependence functions is key in determining ![](http://latex.codecogs.com/gif.latex?W_%7Bbp%7D), particularly:

![](http://latex.codecogs.com/gif.latex?%5Cphi(W_%7Bbp%7D)%5Crho(W_%7Bbp%7D)%3D%5Cfrac%7B2%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5Cbig(%5Cmu_I(%5Cmu_N+%5Csigma)%5Cbig)%7D%7B%5Calpha%5Comega%5E2%5Ctheta%5Cbeta%20HmvU%20I_P%7D)

### MDA coverage necessary to reach the breakpoint

Another useful outcome is the MDA coverage necessary to reach the breakpoint. We can model MDA as a reduction in the mean worm burden at the following time step in terms of the MDA coverage and drug efficacy, *ϵ*:

![](http://latex.codecogs.com/gif.latex?W_%7Bt+1%7D%3D(1-%5Cepsilon)cvrgW_t+(1-cvrg)W_t)

If our goal is to reduce the worm burden to or below the breakpoint (i.e. ![](http://latex.codecogs.com/gif.latex?W_%7Bt+1%7D%5Cleq%20W_%7Bbp%7D)) with MDA, we can estimate the coverage required as:

![](http://latex.codecogs.com/gif.latex?cvrg%3D%5Cfrac%7B%5Cfrac%7BW_%7Bbp%7D%7D%7BW_t%7D-1%7D%7B-%5Cepsilon%7D)

NEEDS UPDATING
==============

![](http://latex.codecogs.com/gif.latex?W_%7Bbp%7D%3D%5Cfrac%7B%5Calpha%5Comega%5Ctheta%20K%5Csigma%5CLambda%5CBig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D-%5Cfrac%7B%5Cmu_N%7D%7Br%7D-%5Cfrac%7B%5CLambda%7D%7Br%7D%5CBig)%7D%7B%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5Cbig(%5Cmu_I(%5Cmu_N+%5Csigma+%5CLambda)+%5Csigma%5CLambda%5Cbig)%5Cbig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D)

Obviously the man-to-snail FOI, *Λ* is a major determinant of the breakpoint. As noted previously, we can estimate it given infected snail prevalence,

![](http://latex.codecogs.com/gif.latex?%5CLambda%3D%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D)

And express the breakpoint in terms of infected snail prevalence as:

![](http://latex.codecogs.com/gif.latex?W_%7Bbp%7D(I_P)%3D%5Cfrac%7B%5Calpha%5Comega%5Ctheta%20K%5Csigma%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D%5CBig(1+%5Cfrac%7B%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D%7D%7B%5Cmu_N+%5Csigma%7D-%5Cfrac%7B%5Cmu_N%7D%7Br%7D-%5Cfrac%7B%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D%7D%7Br%7D%5CBig)%7D%7B%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5Cbig(%5Cmu_I(%5Cmu_N+%5Csigma+%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D)+%5Csigma%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D%5Cbig)%5Cbig(1+%5Cfrac%7B%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D)

Which simplifies (somewhat) to:

![](http://latex.codecogs.com/gif.latex?W_%7Bbp%7D(I_P)%3D%5Cfrac%7B%5Calpha%5Comega%5Ctheta%20KI_P%5Cmu_I(%5Cmu_N+%5Csigma)%5CBig(1+%5Cfrac%7BI_P%5Cmu_I%7D%7B%5Csigma-I_P(%5Cmu_I+%5Csigma)%7D-%5Cfrac%7B%5Cmu_N%7D%7Br%7D-%5Cfrac%7BI_P%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7Br(%5Csigma-I_P(%5Cmu_I+%5Csigma))%7D%5CBig)%7D%7B%5Cbig(1-%5Cfrac%7BI_P%5Cmu_I%7D%7B%5Csigma%7D-I_P%5Cbig)%5Cbig(%5Cmu_W+%5Cmu_H%5Cbig)%5Cbig(%5Cmu_I(%5Cmu_N+%5Csigma+%5Cfrac%7BI_P%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Csigma-I_P(%5Cmu_I+%5Csigma)%7D)+%5Cfrac%7BI_P%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B1-%5Cfrac%7BI_P%5Cmu_I%7D%7B%5Csigma%7D-I_P%7D%5Cbig)%5Cbig(1+%5Cfrac%7BI_P%5Cmu_I%7D%7B%5Csigma-I_P(%5Cmu_I+%5Csigma)%7D%5Cbig)%7D)
