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
  r <- pars["r"]
  K <- pars["K"]
  mu_N <- pars["mu_N"]
  sigma <- pars["sigma"]
  mu_I <- pars["mu_I"]
  mu_W <- pars["mu_W"]
  H <- pars["H"]
  mu_H <- pars["mu_H"]
  theta <- pars["theta"]
  U <- pars["U"]
  m <- pars["m"]
  v <- pars["v"]
  omega <- pars["omega"]
  alpha <- pars["alpha"]
  beta <- pars["beta"]
  cvrg <- pars["cvrg"]
  zeta <- pars["zeta"]
  xi <- pars["xi"]

    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]

    N=S+E+I

  #Clumping parameters based on worm burden
    k_W = k_w_fx(W)

  #Miracidial estimate from treated and untreated populations assuming 1:1 sex ratio, mating probability, density dependence
    M = 0.5*W*H*phi_Wk(W, k_W)*rho_Wk(W, zeta, k_W)*omega*U*m*v

    dSdt= r*(1-(N/K))*(S+E) - (mu_N + (beta*M)/N)*S #Susceptible snails

    dEdt= ((beta*M)/N)*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - mu_I*I #Infected snails

    #worm burden in human
    dWdt= (alpha*omega*theta*I) - ((mu_W+mu_H)*W)



    return(list(c(dSdt,dEdt,dIdt,dWdt)))
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

  r <- pars["r"]
  K <- pars["K"]
  mu_N <- pars["mu_N"]
  sigma <- pars["sigma"]
  mu_I <- pars["mu_I"]
  mu_W <- pars["mu_W"]
  H <- pars["H"]
  mu_H <- pars["mu_H"]
  theta <- pars["theta"]
  U <- pars["U"]
  m <- pars["m"]
  v <- pars["v"]
  omega <- pars["omega"]
  alpha <- pars["alpha"]
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
    k_Wt = k_w_fx(Wt)
    k_Wu = k_w_fx(Wu)

  #Miracidial estimate from treated and untreated populations assuming 1:1 sex ratio, mating probability, density dependence
    M = 0.5*Wt*H*cvrg*phi_Wk(Wt, k_Wt)*rho_Wk(Wt, zeta, k_Wt)*omega*U*m*v +
        0.5*Wu*H*(1-cvrg)*phi_Wk(Wu, k_Wu)*rho_Wk(Wu, zeta, k_Wu)*omega*U*m*v


    dSdt= r*(1-(N/K))*(S+E) - (mu_N + (beta*M)/N)*S #Susceptible snails

    dEdt= ((beta*M)/N)*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - mu_I*I #Infected snails

    #worm burden in human
    dWtdt= (alpha*omega*theta*I) - ((mu_W+mu_H)*Wt)
    dWudt= (alpha*omega*theta*I) - ((mu_W+mu_H)*Wu)



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

  phi_Wt = ifelse(Wt == 0, 0, phi_Wk(W = Wt, k = k_w_fx(Wt)))
  phi_Wu = ifelse(Wu == 0, 0, phi_Wk(W = Wu, k = k_w_fx(Wu)))

#Miracidia as sum of contribution from treated and untreated populations
  M = 0.5*Wt*p["H"]*p["cvrg"]*phi_Wt*rho_Wk(Wt, p["zeta"], k_w_fx(Wt))*p["omega"]*p["U"]*p["m"]*p["v"] +
      0.5*Wu*p["H"]*(1-p["cvrg"])*phi_Wu*rho_Wk(Wu, p["zeta"], k_w_fx(Wu))*p["omega"]*p["U"]*p["m"]*p["v"]

  return(c(p["r"] * (1-N/p["K"]) * (S + E),   #Snail birth
           p["mu_N"] * S,        #Susceptible snail death
           S*(p["beta"]*M/N),  #Snail exposure
           p["mu_N"] * E,       #Exposed snail dies
           p["sigma"] * E,      #Exposed snail becomes infected
           p["mu_I"] * I,   #Infected snail dies
           p["alpha"] * p["theta"] * p["omega"] * I * gam_Wxi(Wt, p["xi"]), #infected snail produces adult worm in treated population
           p["alpha"] * p["theta"] * p["omega"] * I * gam_Wxi(Wu, p["xi"]), #infected snail produces adult worm in untreated population
           (p["mu_W"] + p["mu_H"]) * Wt,    #Adult worm in treated population dies
           (p["mu_W"] + p["mu_H"]) * Wu))    #Adult worm in untreated population dies

}

#' Simulate the stochastic version of age-stratified schistosomiasis model
#'
#' Uses the `schisto_stoch_mod` function to simulate the stochastic model through time
#'
#' @param nstart starting values of state variables, must be whole numbers
#' @param transitions list of all possible transitions between state variables, defaults to transitions possible in `schisto_stoch_mod`
#' @param sfx function to feed into adaptive tau, uses `schisto_stoch_mod` as default
#' @param params named vector of parameter values used in model
#' @param tf numeric of total time to run simulation
#' @param events_df Data frame of events such as MDA with columns "var", "time", "value", and "method"
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
                                                     c(Wt = 1),            #Infected snail emits cercaria that produces an adult worm in treated population
                                                     c(Wu = 1),            #Infected snail emits cercaria that produces an adult worm in untreated population
                                                     c(Wt = -1),           #Adult worm in the treated population dies
                                                     c(Wu = -1)),           #Adult worm in the untreated population dies
                                  sfx = schisto_stoch_mod,
                                  params,
                                  tf,
                                  events_df = NA){
  if(class(events_df) == "data.frame"){
    n_parts <- length(unique(events_df$time)) + 2
    n_vars <- length(unique(events_df$var))

    stoch_sim <- list()

    #Start with 1 year run in
    stoch_sim[[1]] <- adaptivetau::ssa.adaptivetau(nstart, transitions, sfx, params, 365)

    for(i in 1:(n_parts-2)){
      # reset initial states from end of last sim section
      init = stoch_sim[[i]][dim(stoch_sim[[i]])[1],c(2:6)]

      # event occurrence on affected state variables
      for(j in 1:n_vars){
        init[[as.character(events_df[["var"]][((i-1)*n_vars)+j])]] =
          round(init[as.character(events_df[["var"]][((i-1)*n_vars)+j])]*events_df[["value"]][((i-1)*n_vars)+j])
      }

      #stochastic sim for allotted time (time between event occurrences in events_df)
      t_sim <- ifelse(i == 1,
                      events_df[["time"]][i],
                      events_df[["time"]][i*n_vars+1] - events_df[["time"]][i*n_vars])

      stoch_sim[[i+1]] = adaptivetau::ssa.adaptivetau(init, transitions, sfx, params, t_sim)

    #adjust time in section of full simulation
      #remove very last observation
      stoch_sim[[i+1]] <- stoch_sim[[i+1]][-nrow(stoch_sim[[i+1]]),]

      stoch_sim[[i+1]][,"time"] = stoch_sim[[i+1]][,"time"] + events_df[["time"]][i*n_vars]

    }

    #Correct time on one year run in simulation
    stoch_sim[[1]] <- stoch_sim[[1]][-nrow(stoch_sim[[1]]),]

    #Run for remaining time allotted
    stoch_sim[[n_parts]] <- adaptivetau::ssa.adaptivetau(init, transitions, sfx, params, tf - max(events_df[["time"]]))

    #adjust time for remaining sim
    stoch_sim[[n_parts]][,"time"] = stoch_sim[[n_parts]][,"time"] + max(events_df[["time"]])

    return(as.data.frame(do.call(rbind, stoch_sim)))

  } else {
    return(as.data.frame(adaptivetau::ssa.adaptivetau(nstart, transitions, sfx, params, tf)))
  }

}
