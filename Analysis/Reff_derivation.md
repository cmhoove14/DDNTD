*R*<sub>*e**f**f*</sub> Derivation
================

Basic Schistosomiasis Model
---------------------------

Beginning with a basic schistosomiasis model with *S* − *E* − *I* infection dynamics and logistic population growth among the intermediate host snail population and human infection modeled via the state variable *W* representing the mean worm burden across the human population, assumed to be negative binomially distributed with clumping parameter *κ*. We have four ODEs:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdS%7D%7Bdt%7D%3Df_N%5Cbig(1-%5Cfrac%7BN%7D%7BK%7D%5Cbig)%5Cbig(S+E%5Cbig)-%5Cmu_N%20S-%5Cbeta%20M(W)S)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdE%7D%7Bdt%7D%3D%5Cbeta%20M(W)S-(%5Cmu_N+%5Csigma)E)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdI%7D%7Bdt%7D%3D%5Csigma%20E%20-%20%5Cmu_I%20I)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5Clambda%5Crho(W,%20%5Ckappa,%20%5Cxi)%20I-(%5Cmu_W+%5Cmu_H)W)

Where *N* = *S* + *E* + *I*, *ρ*(*W*, *κ*, *ξ*) is an estimate of the density dependent probability of establishment of adult worms due to acquired immunity, and *M*(*W*) is an estimate of miracidia released into the population and is estimated as:

![](http://latex.codecogs.com/gif.latex?M(W)%3D0.5W%5Cphi(W,%20%5Ckappa)%5Cgamma(W,%20%5Ckappa,%20%5Calpha)H)

with *ϕ*(*W*, *κ*) an estimate of the mating probability and *γ*(*W*, *κ*, *α*) an estimate of density dependent fecundity of adult worms.

Parameter symbology, values and definitions are shown in table 1.

<table>
<colgroup>
<col width="10%" />
<col width="7%" />
<col width="16%" />
<col width="66%" />
</colgroup>
<thead>
<tr class="header">
<th align="left"></th>
<th align="left">Value</th>
<th align="left">Units</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><span class="math inline"><em>f</em><sub><em>N</em></sub></span></td>
<td align="left">10</td>
<td align="left"><span class="math inline"><em>S</em><em>d</em><sup>−1</sup></span></td>
<td align="left">Snail fecundity rate</td>
</tr>
<tr class="even">
<td align="left"><span class="math inline"><em>K</em></span></td>
<td align="left">50</td>
<td align="left"><span class="math inline"><em>N</em><em>m</em><sup>−2</sup></span></td>
<td align="left">Snail environmental carrying capacity</td>
</tr>
<tr class="odd">
<td align="left"><span class="math inline"><em>μ</em><sub><em>N</em></sub></span></td>
<td align="left">60</td>
<td align="left"><span class="math inline"><em>N</em><em>d</em><sup>−1</sup></span></td>
<td align="left">Snail mortality rate</td>
</tr>
<tr class="even">
<td align="left"><span class="math inline"><em>σ</em></span></td>
<td align="left">40</td>
<td align="left"><span class="math inline"><em>E</em><em>d</em><sup>−1</sup></span></td>
<td align="left">Pre-patent period</td>
</tr>
<tr class="odd">
<td align="left"><span class="math inline"><em>μ</em><sub><em>I</em></sub></span></td>
<td align="left">10</td>
<td align="left"><span class="math inline"><em>I</em><em>d</em><sup>−1</sup></span></td>
<td align="left">Excess mortality of infected snails</td>
</tr>
<tr class="even">
<td align="left"><span class="math inline"><em>μ</em><sub><em>W</em></sub></span></td>
<td align="left">3.3</td>
<td align="left"><span class="math inline"><em>W</em><em>y</em><sup>−1</sup></span></td>
<td align="left">Adult parasite mortality rate</td>
</tr>
<tr class="odd">
<td align="left"><span class="math inline"><em>H</em></span></td>
<td align="left">1.5</td>
<td align="left"><span class="math inline"><em>H</em><em>m</em><sup>−2</sup></span></td>
<td align="left">Human host population</td>
</tr>
<tr class="even">
<td align="left"><span class="math inline"><em>μ</em><sub><em>H</em></sub></span></td>
<td align="left">60</td>
<td align="left"><span class="math inline"><em>H</em><em>y</em><sup>−1</sup></span></td>
<td align="left">Human host mortality rate</td>
</tr>
<tr class="odd">
<td align="left"><span class="math inline"><em>λ</em></span></td>
<td align="left">4820</td>
<td align="left"><span class="math inline"><em>W</em><em>I</em><sup>−1</sup><em>d</em><sup>−1</sup></span></td>
<td align="left">Man-to-snail transmission parameter</td>
</tr>
<tr class="even">
<td align="left"><span class="math inline"><em>β</em></span></td>
<td align="left">625000</td>
<td align="left"><span class="math inline"><em>E</em><em>S</em><sup>−1</sup><em>d</em><sup>−1</sup></span></td>
<td align="left">Snail-to-man transmission parameter</td>
</tr>
<tr class="odd">
<td align="left"><span class="math inline"><em>α</em></span></td>
<td align="left">0.005</td>
<td align="left">unitless</td>
<td align="left">Adult parasite density dependent fecundity parameter</td>
</tr>
<tr class="even">
<td align="left"><span class="math inline"><em>ξ</em></span></td>
<td align="left">0.0018</td>
<td align="left">unitless</td>
<td align="left">Density dependent parasite establishment acquired immunity parameter</td>
</tr>
</tbody>
</table>

*R*<sub>*e**f**f*</sub> Derivation
----------------------------------

### Dimensionality reduction: fast snail infection dynamics

We begin with the assumption that the rate of change of the intermediate host infection dynamics is fast compared to the adult parasites and therefore reaches an equilibirum, i.e. $\\frac{dS}{dt}=\\frac{dE}{dt}=\\frac{dI}{dt}=0$. We then solve for the equilibirum infected snail population, *I*<sup>\*</sup>:

#### Solve for *I*<sup>\*</sup>

![](http://latex.codecogs.com/gif.latex?%5Csigma%20E%5E*%20-%20%5Cmu_I%20I%5E*%3D0)

![](http://latex.codecogs.com/gif.latex?I%5E*%3D%5Cfrac%7B%5Csigma%20E%5E*%7D%7B%5Cmu_I%7D)

#### Solve for *E*<sup>\*</sup>

![](http://latex.codecogs.com/gif.latex?%5Cbeta%20M(W)S%5E*-(%5Cmu_N+%5Csigma)E%5E*%3D0)

![](http://latex.codecogs.com/gif.latex?E%5E*%3D%5Cfrac%7B%5Cbeta%20M(W)S%5E*%7D%7B%5Cmu_N+%5Csigma%7D)

#### Solve for *S*<sup>\*</sup>

First, we write the equilibrium total snail population, *N*<sup>\*</sup>, in terms of *S*<sup>\*</sup> by substituting for *E*<sup>\*</sup> and *I*<sup>\*</sup> from *N*<sup>\*</sup> = *S*<sup>\*</sup> + *E*<sup>\*</sup> + *I*<sup>\*</sup>: ![](http://latex.codecogs.com/gif.latex?N%5E*%3DS%5E*+%5Cfrac%7B%5Cbeta%20M(W)S%5E*%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%20E%5E*%7D%7B%5Cmu_I%7D)

![](http://latex.codecogs.com/gif.latex?N%5E*%3DS%5E*+%5Cfrac%7B%5Cbeta%20M(W)S%5E*%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5Cbeta%20M(W)S%5E*%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D)

![](http://latex.codecogs.com/gif.latex?N%5E*%3DS%5E*%5CBig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5CBig))

Then, with ![](http://latex.codecogs.com/gif.latex?0%3Df_N%5Cbig(1-%5Cfrac%7BN%5E*%7D%7BK%7D%5Cbig)%5Cbig(S%5E*+E%5E*%5Cbig)-%5Cmu_N%20S%5E*-%5Cbeta%20M(W)S%5E*)

we substitute for *E*<sup>\*</sup> and *N*<sup>\*</sup> and divide by *S*<sup>\*</sup> to arrive at:

![](http://latex.codecogs.com/gif.latex?f_N%5CBigg(1-%5Cfrac%7BS%5E*%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5Cbig)%7D%7BK%7D%5CBigg)%5CBigg(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBigg)-%5Cmu_N-%5Cbeta%20M(W)%3D0)

Solving for *S*<sup>\*</sup>:

![](http://latex.codecogs.com/gif.latex?1-%5Cfrac%7BS%5E*%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5Cbig)%7D%7BK%7D%3D%5Cfrac%7B%5Cmu_N+%5Cbeta%20M(W)%7D%7Bf_N%5CBig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%7D)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BS%5E*%5CBig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5CBig)%7D%7BK%7D%3D1-%5Cfrac%7B%5Cmu_N+%5Cbeta%20M(W)%7D%7Bf_N%5CBig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%7D)

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5CBigg(1-%5Cfrac%7B%5Cmu_N+%5Cbeta%20M(W)%7D%7Bf_N%5CBig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%7D%5CBigg)%5CBigg(%5Cfrac%7BK%7D%7B%5CBig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5CBig)%7D%5CBigg))

Now want to try and simplify this. To start, we'll assign ![](http://latex.codecogs.com/gif.latex?C_1%3D%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D) to give:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5CBigg(1-%5Cfrac%7B%5Cmu_N+%5Cbeta%20M(W)%7D%7Bf_N%5Cbig(1+C_1%5Cbig)%7D%5CBigg)%5CBigg(%5Cfrac%7BK%7D%7B%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D%5CBigg))

Then distribute:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%7D%7B%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D-%5Cfrac%7BK%5Cbig(%5Cmu_N+%5Cbeta%20M(W)%5Cbig)%7D%7Bf_N%5Cbig(1+C_1%5Cbig)%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D)

Then multiply the LHS by ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7Bf_N%5Cbig(1+C_1%5Cbig)%7D%7Bf_N%5Cbig(1+C_1%5Cbig)%7D)

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(f_N%5Cbig(1+C_1%5Cbig)%5Cbig)-K%5Cbig(%5Cmu_N+%5Cbeta%20M(W)%5Cbig)%7D%7Bf_N%5Cbig(1+C_1%5Cbig)%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D)

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(f_N+f_N%20C_1-%5Cmu_N-%5Cbeta%20M(W)%5Cbig)%7D%7Bf_N%5Cbig(1+C_1%5Cbig)%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D)

This is basically where we left it in the [Elimination Feasibility paper](https://doi.org/10.1371/journal.pntd.0006794), but I think it might be possible to simplify a bit more by first factoring out *f*<sub>*N*</sub> from the numerator and then canceling:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Ccancel%7Bf_N%7D%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Ccancel%7Bf_N%7D%5Cbig(1+C_1%5Cbig)%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D)

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cbig(1+C_1%5Cbig)%5Cbig(1+C_1+%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig)%7D)

Then distribute in the denominator:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cfrac%7B2%5Csigma%20C_1%5E2%7D%7B%5Cmu_I%7D+%5Cfrac%7B3%5Csigma%20C_1%7D%7B%5Cmu_I%7D+1%7D)

And then factor out ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D) from the denominator:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5Cbig(2C_1+3+%5Cfrac%7B%5Cmu_I%7D%7B%5Csigma%20C_1%7D%5Cbig)%7D)

Now with the rate of change of the mean worm burden:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5Clambda%5Crho(W,%20%5Ckappa,%20%5Cxi)%20I%5E*-(%5Cmu_W+%5Cmu_H)W)

And: ![](http://latex.codecogs.com/gif.latex?I%5E*%3D%5CBig(%5Cfrac%7B%5Csigma%7D%7B%5Cmu_I%7D%5CBig)%5CBig(%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)S%5E*) and ![](http://latex.codecogs.com/gif.latex?C_1%3D%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D)

We get:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5Clambda%5Crho(W,%20%5Ckappa,%20%5Cxi)%20%5CBig(%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%5CBig)S%5E*-(%5Cmu_W+%5Cmu_H)W)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5Clambda%5Crho(W,%20%5Ckappa,%20%5Cxi)%20%5CBig(%5Ccancel%7B%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%7D%5CBig)%5CBig(%5Cfrac%7BK%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Ccancel%7B%5Cfrac%7B%5Csigma%20C_1%7D%7B%5Cmu_I%7D%7D%5Cbig(2C_1+3+%5Cfrac%7B%5Cmu_I%7D%7B%5Csigma%20C_1%7D%5Cbig)%7D%5CBig)-(%5Cmu_W+%5Cmu_H)W)

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5CBig(%5Cfrac%7BK%5Clambda%5Crho%5Cbig(W,%20%5Ckappa,%20%5Cxi%5Cbig)%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cbig(2C_1+3+%5Cfrac%7B%5Cmu_I%7D%7B%5Csigma%20C_1%7D%5Cbig)%7D%5CBig)-(%5Cmu_W+%5Cmu_H)W)

Now factor out (*μ*<sub>*W*</sub> + *μ*<sub>*H*</sub>)*W* to get:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D%5CBig((%5Cmu_W+%5Cmu_H)W%5CBig)%5CBig(%5Cfrac%7BK%5Clambda%5Crho%5Cbig(W,%20%5Ckappa,%20%5Cxi%5Cbig)%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cbig(2C_1+3+%5Cfrac%7B%5Cmu_I%7D%7B%5Csigma%20C_1%7D%5Cbig)%5Cbig((%5Cmu_W+%5Cmu_H)W%5Cbig)%7D-1%5CBig))

Given the definition of *R*<sub>*e**f**f*</sub> as the number of adult worms produced by a single adult worm over its lifespan, we interpret (*μ*<sub>*W*</sub> + *μ*<sub>*H*</sub>) as the mean lifespan and therefore have:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdW%7D%7Bdt%7D%3D(%5Cmu_W+%5Cmu_H)W(R_%7Beff%7D-1)) and therefore:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D%3D%5Cfrac%7BK%5Clambda%5Crho%5Cbig(W,%20%5Ckappa,%20%5Cxi%5Cbig)%5Cbig(1+C_1-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cbig(2C_1+3+%5Cfrac%7B%5Cmu_I%7D%7B%5Csigma%20C_1%7D%5Cbig)%5Cbig((%5Cmu_W+%5Cmu_H)W%5Cbig)%7D)

and plugging back in for *C*<sub>1</sub>:

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D%3D%5Cfrac%7BK%5Clambda%5Crho%5Cbig(W,%20%5Ckappa,%20%5Cxi%5Cbig)%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D-%5Cfrac%7B%5Cmu_N%7D%7Bf_N%7D-%5Cfrac%7B%5Cbeta%20M(W)%7D%7Bf_N%7D%5Cbig)%7D%7B%5Cbig(%5Cfrac%7B2%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+3+%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Csigma(%5Cbeta%20M(W))%7D%5Cbig)%5Cbig((%5Cmu_W+%5Cmu_H)W%5Cbig)%7D)

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D%3D%5Cfrac%7BK%5Clambda%5Crho%5Cbig(W,%20%5Ckappa,%20%5Cxi%5Cbig)%5Cbig(f_N(%5Cmu_N+%5Csigma)+f_N%5Cbeta%20M(W)-%5Cmu_N(%5Cmu_N+%5Csigma)-%5Cbeta%20M(W)(%5Cmu_N+%5Csigma)%5Cbig)%7D%7B%5Cbig(f_N(%5Cmu_N+%5Csigma)%5Cbig)%5Cbig(%5Cfrac%7B2%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+3+%5Cfrac%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%7B%5Csigma(%5Cbeta%20M(W))%7D%5Cbig)%5Cbig((%5Cmu_W+%5Cmu_H)W%5Cbig)%7D)

Then distribute the denominator and factor out ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7Bf_N%7D%7B%5Csigma%5Cbeta%20M(W)%7D)

![](http://latex.codecogs.com/gif.latex?R_%7Beff%7D%3D%5Cfrac%7BK%5Clambda%5Crho%5Cbig(W,%20%5Ckappa,%20%5Cxi%5Cbig)%5Cbig(%5Csigma%5Cbeta%20M(W)%5Cbig)%5Cbig(f_N(%5Cmu_N+%5Csigma)+f_N%5Cbeta%20M(W)-%5Cmu_N(%5Cmu_N+%5Csigma)-%5Cbeta%20M(W)(%5Cmu_N+%5Csigma)%5Cbig)%7D%7Bf_N%5CBig(2%5Csigma%5Cbig(%5Cbeta%20M(W)%5Cbig)%5E2+3%5Csigma%5Cbig(%5Cmu_N+%5Csigma%5Cbig)%5Cbig(%5Cbeta%20M(W)%5Cbig)+%5Cmu_I%5Cbig(%5Cmu_N+%5Csigma%5Cbig)%5E2%5CBig)%5CBig((%5Cmu_W+%5Cmu_H)W%5CBig)%7D)

A separate attempt at simplifying *S*<sup>\*</sup> further that don't end up anywhere
-------------------------------------------------------------------------------------

### Attempt 2

after plugging ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D) back in for *C*<sub>1</sub>:

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(f_N+%5Cfrac%7Bf_N%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D-%5Cmu_N-%5Cbeta%20M(W)%5Cbig)%7D%7Bf_N%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%20%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5Cbig)%7D)

Multiply the three terms in the numerator by ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cmu_N+%5Csigma%7D) to factor out a *μ*<sub>*N*</sub> + *σ*

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(%5Cfrac%7Bf_N(%5Cmu_N+%5Csigma)+f_N%5Cbeta%20M(W)-%5Cmu_N(%5Cmu_N+%5Csigma)-%5Cbeta%20M(W)(%5Cmu_N+%5Csigma)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%7Bf_N%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%5Cbig(1+%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D+%5Cfrac%7B%5Csigma%20%5Cbeta%20M(W)%7D%7B%5Cmu_I(%5Cmu_N+%5Csigma)%7D%5Cbig)%7D)

Now factor out some ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D) from the denominator.

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(%5Cfrac%7Bf_N(%5Cmu_N+%5Csigma)+f_N%5Cbeta%20M(W)-%5Cmu_N(%5Cmu_N+%5Csigma)-%5Cbeta%20M(W)(%5Cmu_N+%5Csigma)%7D%7B%5Cmu_N+%5Csigma%7D%5Cbig)%7D%7Bf_N%5CBig(%5Cfrac%7B%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%5E2%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1+%5Cfrac%7B%5Csigma%7D%7B%5Cmu_I%7D%5CBig)%7D)

Then get rid of a *μ*<sub>*N*</sub> + *σ* in both the numerator and denominator

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(f_N(%5Cmu_N+%5Csigma)+f_N%5Cbeta%20M(W)-%5Cmu_N(%5Cmu_N+%5Csigma)-%5Cbeta%20M(W)(%5Cmu_N+%5Csigma)%5Cbig)%7D%7Bf_N%5CBig(%5Cfrac%7B%5Cbig(%5Cbeta%20M(W)%5Cbig)%5E2%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1+%5Cfrac%7B%5Csigma%7D%7B%5Cmu_I%7D%5CBig)%7D)

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(%5Cmu_N+%5Csigma%5Cbig)%5Cbig(f_N+%5Cfrac%7Bf_N%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D-%5Cmu_N-%5Cbeta%20M(W)%5Cbig)%7D%7Bf_N%5CBig(%5Cfrac%7B%5Cbig(%5Cbeta%20M(W)%5Cbig)%5E2%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1+%5Cfrac%7B%5Csigma%7D%7B%5Cmu_I%7D%5CBig)%7D)

Then redistribute the ![](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cbig(%5Cbeta%20M(W)%5Cbig)%5E2%7D%7B%5Cmu_N+%5Csigma%7D) in the denominator

![](http://latex.codecogs.com/gif.latex?S%5E*%3D%5Cfrac%7BK%5Cbig(%5Cmu_N+%5Csigma%5Cbig)%5Cbig(f_N+%5Cfrac%7Bf_N%5Cbeta%20M(W)%7D%7B%5Cmu_N+%5Csigma%7D-%5Cmu_N-%5Cbeta%20M(W)%5Cbig)%7D%7Bf_N%5CBig(%5Cfrac%7B%5Cbig(%5Cbeta%20M(W)%5Cbig)%5E2%7D%7B%5Cmu_N+%5Csigma%7D%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1%5CBig)%5CBig(%5Cfrac%7B%5Cmu_N+%5Csigma%7D%7B%5Cbeta%20M(W)%7D+1+%5Cfrac%7B%5Csigma%7D%7B%5Cmu_I%7D%5CBig)%7D)
