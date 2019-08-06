#' Barebones schistosomiasis model
#'
#' A dynamic schistosomiasis model with SEI snail infection dynamics, a single mean worm burden population,
#' negative (crowding induced reductions in fecundity)
#' and positive (mating limitation) density dependencies and functionality
#' to simulate mass drug administration, snail control, and other interventions.
#' Note this is a function that is wrapped into `sim_schisto_base_mod`` which
#' should be used to simulate the model, this function is fed into the ode solver
#' from deSolve
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param pars Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_mod <- function(t, n, pars) {
  with(as.list(pars),{

    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]

    N=S+E+I

  #Clumping parameters based on worm burden
    k_W = k_from_log_W(W)

  #Miracidial estimate from treated and untreated populations assuming 1:1 sex ratio, mating probability, density dependence
    M = (0.5*W*H)*phi_Wk(W, k_W)*rho_Wk(W, zeta, k_W)

  #Snail infection dynamics
    dSdt= f_N*(1-(N/K))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

  #worm burden in humans
    dWdt= (lambda*I*gam_Wxi(W, xi)) - ((mu_W+mu_H)*W)


    return(list(c(dSdt,dEdt,dIdt,dWdt)))
  })
}


#' Basic schistosomiasis model
#'
#' A dynamic schistosomiasis model with SEI snail infection dynamics, separate
#' compartments for the mean worm burden in treated and untreated segments of the
#' human population, negative (crowding induced reductions in fecundity)
#' and positive (mating limitation) density dependencies and functionality
#' to simulate mass drug administration, snail control, and other interventions.
#' Note this is a function that is wrapped into `sim_schisto_base_mod`` which
#' should be used to simulate the model, this function is fed into the ode solver
#' from deSolve
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param pars Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_base_mod <- function(t, n, pars) {

  f_N <- pars["f_N"]
  K <- pars["K"]
  mu_N <- pars["mu_N"]
  sigma <- pars["sigma"]
  mu_I <- pars["mu_I"]
  mu_W <- pars["mu_W"]
  H <- pars["H"]
  mu_H <- pars["mu_H"]
  lambda <- pars["lambda"]
  beta <- pars["beta"]
  cvrg <- pars["cvrg"]
  zeta <- pars["zeta"]
  xi <- pars["xi"]

    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]

    N=S+E+I

    W=(cvrg*Wt) + ((1-cvrg)*Wu) #weighting treated and untreated populations

  #Clumping parameters based on worm burden
    k_Wt = k_from_log_W(Wt)
    k_Wu = k_from_log_W(Wu)

  #Miracidial estimate from treated and untreated populations assuming 1:1 sex ratio, mating probability, density dependence
    M = (0.5*Wt*H*cvrg)*phi_Wk(Wt, k_Wt)*rho_Wk(Wt, zeta, k_Wt) +
        (0.5*Wu*H*(1-cvrg))*phi_Wk(Wu, k_Wu)*rho_Wk(Wu, zeta, k_Wu)


    dSdt= f_N*(1-(N/K))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human
    dWtdt= (lambda*I*gam_Wxi(Wt, xi)) - ((mu_W+mu_H)*Wt)
    dWudt= (lambda*I*gam_Wxi(Wu, xi)) - ((mu_W+mu_H)*Wu)



    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
}

#' Simulate a schistosomiasis model
#'
#' Simulate a schisto deterministic model through time using `ode` function
#' given the schisto model deSolve function, parameter set, time frame, starting conditions, and potential events (e.g. MDA)
#'
#' @param nstart Named vector of starting values for state variables
#' @param time Numeric vector of times at which state variables should be estimated
#' @param model Name of the ode function to use, defaults to `schisto_base_mod`
#' @param pars Named vector or list of parameter values
#' @param events_df Data frame of events such as MDA with columns "var", "time", "value", and "method"
#'
#' @return dataframe of state variable values at requested times
#' @export
#'
sim_schisto_mod <- function(nstart,
                            time,
                            model = schisto_base_mod,
                            pars,
                            events_df = NA){
  if(is.na(events_df)){
    as.data.frame(ode(nstart, time, model, pars))
  } else {
    as.data.frame(ode(nstart, time, model, pars,
                      events = list(data = events_df)))
  }
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

  phi_Wt = ifelse(Wt == 0, 0, phi_Wk(W = Wt, k = k_from_log_W(Wt)))
  phi_Wu = ifelse(Wu == 0, 0, phi_Wk(W = Wu, k = k_from_log_W(Wu)))

  return(c(p["f_N"] * (1-N/p["K"]) * (S + E),   #Snail birth
           p["mu_N"] * S,        #Susceptible snail death
           ((0.5 * Wt * p["H"] * p["cvrg"] * phi_Wt * rho_Wk(Wt, p["zeta"], k_from_log_W(Wt))) +
             (0.5 * Wu * p["H"] * (1-p["cvrg"]) * phi_Wu * rho_Wk(Wu, p["zeta"], k_from_log_W(Wu)))) *
              S * p["beta"],  #Snail exposure
           p["mu_N"] * E,       #Exposed snail dies
           p["sigma"] * E,      #Exposed snail becomes infected
           (p["mu_N"] + p["mu_I"]) * I,   #Infected snail dies
           p["lambda"] * I * gam_Wxi(Wt, p["xi"]), #infected snail produces adult worm in treated population
           p["lambda"] * I * gam_Wxi(Wu, p["xi"]), #infected snail produces adult worm in untreated population
           (p["mu_W"] + p["mu_H"]) * Wt,    #Adult worm in treated population dies
           (p["mu_W"] + p["mu_H"]) * Wu))    #Adult worm in untreated population dies

}

