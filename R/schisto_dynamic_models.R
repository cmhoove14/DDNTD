#' Basic schistosomiasis model
#'
#' A dynamic schistosomiasis model with negative (crowding induced reductions in fecundity)
#' and positive (mating limitation) density dependencies and functionality
#' to simulate mass drug administration, snail control, and other interventions.
#' Note this is a function that is wrapped into sim_schisto_base_mod which
#' should be used to simulate the model, this function is fed into the ode solver
#' from deSolve
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param parameters Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_base_mod <- function(t, n, parameters) {
  with(as.list(parameters),{

    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I

    W=(cvrg*Wt) + ((1-cvrg)*Wu) #weighting treated and untreated populations


  #Density dependencies based on worm burden
    f = f_Wgk(W, gamma, k)
    R = R_Wv(W, v)
    phi = phi_Wk(W = W, k = k)  #Mating probability

  #Miracidial estimate assuming 1:1 sex ratio, mating probability, density dependence
    M=((0.5*W*H)*phi*f)#*m*u_H*(v*vq)


    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human
    dWtdt= (lambda*I*R) - ((mu_W+mu_H)*Wt)
    dWudt= (lambda*I*R) - ((mu_W+mu_H)*Wu)



    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  })
}

#' Simulate the base schistosomiasis model
#'
#' Uses the `schisto_base_mod` function to simulate the model through time
#' given parameter set, time frame, starting conditions, and potential events (e.g. MDA)
#'
#' @param nstart Named vector of starting values for state variables
#' @param time Numeric vector of times at which state variables should be estimated
#' @param model Name of the ode function to use, defaults to `schisto_base_mod`
#' @param parameters Named vector or list of parameter values
#' @param events_df Data frame of events such as MDA
#'
#' @return dataframe of state variable values at requested times
#' @export
#'
sim_schisto_base_mod <- function(nstart,
                                 time,
                                 model = schisto_base_mod,
                                 parameters,
                                 events_df){
  as.data.frame(ode(nstart, time, model, parameters,
                    events = list(data = events_df)))
}


#' Stochastic version of the base schistosomiasis model
#'
#' A stochastic schistosomiasis model with same structure as the deterministic model
#' but implemented with the adaptivetau package which uses tau leaping
#' to simulate stochastic transmission. Note this is a function that is wrapped into
#' sim_schisto_stoch_mod which should be used to simulate the model,
#' this function is fed into the tau leaping algorithm from adaptivetau
#'
#' @param x vector of state variables
#' @param p named vector of parameters
#' @param t numeric indicating length of simulation
#'
#' @return vector of transition probabilities
#' @export
#'
schisto_stoch_mod <- function(x, p, t){
  S = x['S']
  E = x['E']
  I = x['I']
  N = S + I + E
  Wt = x['Wt']
  Wu = x['Wu']
  W = p$cvrg*Wt+(1-p$cvrg)*Wu

  return(c(p$f_N * (1-N/p$C) * (S + E),   #Snail birth
           p$mu_N * S,        #Susceptible snail death
           p$beta * 0.5 * W * p$H * S * phi_Wk(W = W, k = p$k) * f_Wgk(W, p$gam, p$k),  #Snail exposure
           p$mu_N * E,       #Exposed snail dies
           p$sigma * E,      #Exposed snail becomes infected
           (p$mu_N + p$mu_I) * I,   #Infected snail dies
           p$lamda * I * R_Wv(W, p$v),        #infected snail produces adult worm
           (p$mu_W + p$mu_H) * Wt,    #Adult worm in treated population dies
           (p$mu_W + p$mu_H) * Wu))    #Adult worm in untreated population dies

}

#' Simulate the stochastic schistosomiasis model
#'
#' Uses the `schisto_stoch_mod` function to simulate the stochastic model through time
#'
#' @param nstart starting values of state variables, must be whole numbers
#' @param transitions list of all possible transitions between state variables, defaults to transitions possible in `schisto_stoch_mod`
#' @param sfx function to feed into adaptive tau, uses `schisto_stoch_mod` as default
#' @param params named vector of parameter values used in model
#' @param tf numeric of total time to run simulation
#'
#' @return matrix of state variables and values at time points up until tf
#' @export
#'
sim_schisto_stoch_mod <- function(nstart,
                                  transitions = list(c(S = 1),             #New snail born
                                                     c(S = -1),            #Susceptible snail dies
                                                     c(S = -1, E = 1),     #Susceptible snail becomes exposed
                                                     c(E = -1),            #Exposed snail dies
                                                     c(E = -1, I = 1),     #Exposed snail becomes Infected
                                                     c(I = -1),            #Infected snail dies
                                                     c(Wt = 1, Wu = 1),    #Infected snail emits cercaria that produces an adult worm
                                                     c(Wt = -1),           #Adult worm in the treated population dies
                                                     c(Wu = -1)),           #Adult worm in the untreated population dies
                                  sfx = schisto_stoch_mod,
                                  params,
                                  tf){
  adaptivetau::ssa.adaptivetau(nstart, transitions, sfx, params, tf)
}
