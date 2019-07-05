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
#' @param parameters Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_mod <- function(t, n, parameters) {
  with(as.list(parameters),{

    S=n[1]
    E=n[2]
    I=n[3]
    W=n[4]

    N=S+E+I

  #Clumping parameters based on worm burden
    k_W = k_from_log_W(W)

  #Miracidial estimate from treated and untreated populations assuming 1:1 sex ratio, mating probability, density dependence
    M = (0.5*W*H)*phi_Wk(W, k_W)*f_Wgk(W, gamma, k_W)

  #Snail infection dynamics
    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

  #worm burden in humans
    dWdt= (lambda*I*R_Wv(W, xi)) - ((mu_W+mu_H)*W)


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
#' @param parameters Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_base_mod <- function(t, n, parameters) {

  f_N <- parameters["f_N"]
  C <- parameters["C"]
  mu_N <- parameters["mu_N"]
  sigma <- parameters["sigma"]
  mu_I <- parameters["mu_I"]
  mu_W <- parameters["mu_W"]
  H <- parameters["H"]
  mu_H <- parameters["mu_H"]
  lambda <- parameters["lambda"]
  beta <- parameters["beta"]
  cvrg <- parameters["cvrg"]
  gamma <- parameters["gamma"]
  xi <- parameters["xi"]

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
    M = (0.5*Wt*H*cvrg)*phi_Wk(Wt, k_Wt)*f_Wgk(Wt, gamma, k_Wt) +
        (0.5*Wu*H*(1-cvrg))*phi_Wk(Wu, k_Wu)*f_Wgk(Wu, gamma, k_Wu)


    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human
    dWtdt= (lambda*I*R_Wv(Wt, xi)) - ((mu_W+mu_H)*Wt)
    dWudt= (lambda*I*R_Wv(Wu, xi)) - ((mu_W+mu_H)*Wu)



    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
}

#' Age-stratified schistosomiasis model
#'
#' An age-stratified dynamic schistosomiasis model with treated and untreated
#' school-aged children (SAC) and adult populations,
#' negative (crowding induced reductions in fecundity)
#' and positive (mating limitation) density dependencies and functionality
#' to simulate mass drug administration, snail control, and other interventions.
#' Also incorporates non-linear snail-to-man transmission dynamics as derived in Gurarie et al 2018
#' Note this is a function that is wrapped into `sim_schisto_age_strat_mod`` which
#' should be used to simulate the model, this function is fed into the ode solver
#' from deSolve
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param parameters Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_age_strat_mod <- function(t, n, parameters) {
  with(as.list(parameters),{

    S=n[1]
    E=n[2]
    I=n[3]
    Wt_SAC=n[4]
    Wu_SAC=n[5]
    Wt_adult=n[6]
    Wu_adult=n[7]

    #Total snail population
      N=S+E+I

    #weighting treated and untreated populations among SAC
      W_SAC = (cvrg_SAC*Wt_SAC) + ((1-cvrg_SAC)*Wu_SAC)

    #weighting treated and untreated populations among adults
      W_adult = (cvrg_adult*Wt_adult) + ((1-cvrg_adult)*Wu_adult)

    #Weighting SAC and adult populations
      W_tot = W_SAC*prop_SAC + W_adult*prop_adult


    #Update clumping parameter, k from estimate of eggs burden per 10mL estimate
      k_t_SAC = k_from_log_W(Wt_SAC)
      k_u_SAC = k_from_log_W(Wu_SAC)
      k_t_adult = k_from_log_W(Wt_adult)
      k_u_adult = k_from_log_W(Wu_adult)

    #Estimate mating probability within each strata
      phi_Wt_SAC = phi_Wk(W = Wt_SAC, k = k_t_SAC)  #Mating probability in treated SAC population
      phi_Wu_SAC = phi_Wk(W = Wu_SAC, k = k_u_SAC)  #Mating probability in untreated SAC population
      phi_Wt_adult = phi_Wk(W = Wt_adult, k = k_t_adult)  #Mating probability in treated adult population
      phi_Wu_adult = phi_Wk(W = Wu_adult, k = k_u_adult)  #Mating probability in untreated adult population


    #Estimate mean eggs produced per person in each strata as product of
    # mated female worms,
    # eggs produced per female worm per 10mL urine,
    # and reduction in fecundity due to crowding
      eggs_Wt_SAC = 0.5*(Wt_SAC*phi_Wt_SAC) * m * f_Wgk(Wt_SAC, gamma, k_t_SAC)

      eggs_Wu_SAC = 0.5*(Wu_SAC*phi_Wu_SAC) * m * f_Wgk(Wu_SAC, gamma, k_u_SAC)

      eggs_Wt_adult = 0.5*(Wt_adult*phi_Wt_adult) * m * f_Wgk(Wt_adult, gamma, k_t_adult)

      eggs_Wu_adult = 0.5*(Wu_adult*phi_Wu_adult) * m * f_Wgk(Wu_adult, gamma, k_u_adult)

    #Estimate miracidia produced by each strata as product of
    # mean eggs per 10 mL urine for individuals in each strata,
    # number of people in each strata,
    # egg viability
    # mean mL urine produced by an average individual in each group/10,
    # contamination coefficient for SAC/adults,

      M_Wt_SAC = eggs_Wt_SAC * ((H*prop_SAC)*cvrg_SAC) * v * u_SAC * rho_SAC

      M_Wu_SAC = eggs_Wu_SAC * ((H*prop_SAC)*(1-cvrg_SAC)) * v * u_SAC * rho_SAC

      M_Wt_adult = eggs_Wt_adult * ((H*prop_adult)*cvrg_adult) * v * u_adult * rho_adult

      M_Wu_adult = eggs_Wu_adult * ((H*prop_adult)*(1-cvrg_adult)) * v * u_adult * rho_adult

    # Estimate total miracidia entering snail habitat
      M_tot = M_Wt_SAC + M_Wu_SAC + M_Wt_adult + M_Wu_adult

    # Snail infection dynamics
      dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*(1-exp(-M_tot/N))*S #Susceptible snails

      dEdt= beta*(1-exp(-M_tot/N))*S - (mu_N+sigma)*E #Exposed snails

      dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human populations
      dWt_SACdt= (omega_SAC*lambda*I*R_Wv(Wt_SAC, xi)) - (mu_W*Wt_SAC)
      dWu_SACdt= (omega_SAC*lambda*I*R_Wv(Wu_SAC, xi)) - (mu_W*Wu_SAC)
      dWt_adultdt= (omega_adult*lambda*I*R_Wv(Wt_adult, xi)) - (mu_W*Wt_adult)
      dWu_adultdt= (omega_adult*lambda*I*R_Wv(Wu_adult, xi)) - (mu_W*Wu_adult)

    return(list(c(dSdt,dEdt,dIdt,
                  dWt_SACdt,dWu_SACdt,
                  dWt_adultdt, dWu_adultdt)))
  })
}

#' Simulate a schistosomiasis model
#'
#' Simulate a schisto deterministic model through time using `ode` function
#' given the schisto model deSolve function, parameter set, time frame, starting conditions, and potential events (e.g. MDA)
#'
#' @param nstart Named vector of starting values for state variables
#' @param time Numeric vector of times at which state variables should be estimated
#' @param model Name of the ode function to use, defaults to `schisto_base_mod`
#' @param parameters Named vector or list of parameter values
#' @param events_df Data frame of events such as MDA with columns "var", "time", "value", and "method"
#'
#' @return dataframe of state variable values at requested times
#' @export
#'
sim_schisto_mod <- function(nstart,
                            time,
                            model = schisto_base_mod,
                            parameters,
                            events_df = NA){
  if(is.na(events_df)){
    as.data.frame(ode(nstart, time, model, parameters))
  } else {
    as.data.frame(ode(nstart, time, model, parameters,
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

  return(c(p$f_N * (1-N/p$C) * (S + E),   #Snail birth
           p$mu_N * S,        #Susceptible snail death
           ((0.5 * Wt * p$H * p$cvrg * phi_Wt * f_Wgk(Wt, p$gam, k_from_log_W(Wt))) +
             (0.5 * Wu * p$H * (1-p$cvrg) * phi_Wu * f_Wgk(Wu, p$gam, k_from_log_W(Wu)))) *
              S * p$beta,  #Snail exposure
           p$mu_N * E,       #Exposed snail dies
           p$sigma * E,      #Exposed snail becomes infected
           (p$mu_N + p$mu_I) * I,   #Infected snail dies
           p$lambda * I * R_Wv(Wt, p$xi), #infected snail produces adult worm in treated population
           p$lambda * I * R_Wv(Wu, p$xi), #infected snail produces adult worm in untreated population
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
  if(is.na(events_df)){
    return(as.data.frame(adaptivetau::ssa.adaptivetau(nstart, transitions, sfx, params, tf)))
  } else {
    n_parts <- nrow(events_df) + 2

    stoch_sim <- list()

    #Start with 1 year run in
    stoch_sim[[1]] <- adaptivetau::ssa.adaptivetau(nstart, transitions, sfx, params, 365)

    for(i in 1:nrow(events_df)){
      # reset initial states
      init = setNames(as.numeric(stoch_sim[[i]][dim(stoch_sim[[i]])[1],c(2:6)]),
                      colnames(stoch_sim[[i]])[c(2:6)])

      # event occurrence
      init[[as.character(events_df[["var"]][i])]] =
        round(init[as.character(events_df[["var"]][i])]*events_df[["value"]][i])

      #stochastic sim for allotted time (time between event occurrences in events_df)
      t_sim <- ifelse(i == 1,
                      events_df[["time"]][i],
                      events_df[["time"]][i] - events_df[["time"]][i-1])

      stoch_sim[[i+1]] = adaptivetau::ssa.adaptivetau(init, transitions, sfx, params, t_sim)

      #adjust time in section of full simulation
      stoch_sim[[i+1]] <- stoch_sim[[i+1]][-nrow(stoch_sim[[i+1]]),]

      stoch_sim[[i+1]][,"time"] = stoch_sim[[i+1]][,"time"] + events_df[["time"]][i]

    }

    #Correct time on one year run in simulation
    stoch_sim[[1]] <- stoch_sim[[1]][-nrow(stoch_sim[[1]]),]

    #Run for remaining time allotted
    stoch_sim[[n_parts]] <- adaptivetau::ssa.adaptivetau(init, transitions, sfx, params, tf - max(events_df[["time"]]))

    #adjust time for remaining sim
    stoch_sim[[n_parts]][,"time"] = stoch_sim[[n_parts]][,"time"] + max(events_df[["time"]]) + 365

    return(as.data.frame(do.call(rbind, stoch_sim)))
  }

}
