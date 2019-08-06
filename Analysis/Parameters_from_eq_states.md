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

![](http://latex.codecogs.com/gif.latex?I%5E*%3D%5Cfrac%7BK%5Csigma%5CBig(1-%5Cfrac%7B%5Cmu_N+%5CLambda%7D%7Br%5Cbig(1+%5Cfrac%7B%5CLambda%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%5CBig)%7D%7B%5CBig(%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5CLambda%7D+%5Cmu_I+%5Csigma%5CBig)%7D)

Therefore with infected snail prevalence, *I*<sub>*P*</sub> = *I*<sup>\*</sup>/*N*<sup>\*</sup>, and other snail population and infection parameters as inputs, we can estimate *Λ* as:

![](http://latex.codecogs.com/gif.latex?%5CLambda%3D%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Cfrac%7B%5Csigma%7D%7BI_P%7D-%5Cmu_I-%5Csigma%7D)

Contamination fraction, *Ω*, and miracidia invasion rate, *Λ*<sub>0</sub>, as function of equilibrium worm burden in children and adults
----------------------------------------------------------------------------------------------------------------------------------------

We first assume that differences in child and adult infection rates arise predominately from variation in behavior that affects both exposure and contribution to infection in the same manner. Parameter *Ω* therefore represents the relative exposure/contamination of children to adults, *Ω* = *ω*<sub>*C*</sub>/*ω*<sub>*A*</sub>, which can be estimated from the equilibrium infection ratio *Ω* = *W*<sub>*C*</sub><sup>\*</sup>/*W*<sub>*A*</sub><sup>\*</sup>.

With this parameter, we can estimate the equilibirum miracidial density *M* as:

![](http://latex.codecogs.com/gif.latex?M%3D0.5%5Cmathbf%7BH%7Dm%5Comega_A(W_C%5E*%20h_C%5CPhi(W_C%5E*)%5Crho(W_C%5E*)U_C%5COmega+W_A%5E*%20h_A%5CPhi(W_A%5E*)%5Crho(W_A%5E*)U_A))

Now, since *Λ* has been estimated above, *N*<sup>\*</sup> can be estimated as a function of *Λ*, and we have an equilibrium estimate of miracidial density, *M*, we can estimate the baseline miracidial invasion rate, *Λ*<sub>0</sub>:

![](http://latex.codecogs.com/gif.latex?%5CLambda_0%3D%5Cfrac%7B%5CLambda%7D%7B(1-e%5E%7B-M/N%5E*(%5CLambda)%7D)%7D)

Except in rare circumstances, *Λ* ≈ *Λ*<sub>0</sub> because of the high *M*/*N* ratio in transmission settings that have not been intervened upon (i.e. equilibrium settings)

Worm establishment rate, *α*, as function of equilibrium worm burden
--------------------------------------------------------------------

Assuming the observed mean worm burden in each population fraction prior to intervention is approximately at equilibrium, the worm acquisition rate, *λ*<sub>*i*</sub>*j*, for each group can be estimated as a function of the exposure fractions, *ω*<sub>*a*</sub> and *ω*<sub>*c*</sub>, the infected snail density, *I*<sup>\*</sup> = *I*<sub>*P*</sub>*N*<sup>\*</sup>, the snail cercarial shedding rate, *θ*, and the probability of cercarial establishment per exposure, *α*, with *α* being the only unkown and estimated as:

![](http://latex.codecogs.com/gif.latex?%5Calpha%3D%5Cfrac%7BW_i%5E*(%5Cmu_W+%5Cmu_H_i)%7D%7B%5Comega_i%20I%5E*%5Ctheta%7D)
