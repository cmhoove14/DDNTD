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
  ##standard snail parameters
    r=parameters["r"]             # recruitment rate (from sokolow et al)
    K=parameters["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=parameters["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=parameters["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=parameters["mu_I"]          # Increased mortality rate of infected snails
    theta=parameters["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = parameters["mu_W"]   # death rate of adult worms
    mu_H_A = parameters["mu_H_A"] # death rate of adult humans
    mu_H_C = parameters["mu_H_C"] # death rate of children
    m = parameters["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = parameters["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = parameters["zeta"]       # parameter of fecundity reduction function
    xi = parameters["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H_TC = parameters["H_TC"]         # Total number of treated children
    H_UC = parameters["H_UC"]          # Total number of untreated children
    H_TA = parameters["H_TA"]           # Total number of treated adults
    H_UA = parameters["H_UA"]         # Total number of untreated adults
    U_C = parameters["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = parameters["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = parameters["omega"]          # relative infection risk/contamination of SAC vs adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=parameters["alpha"]       # Cercarial infection probability
    Lambda_0=parameters["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=parameters["beta"]         # second parameter of non-linear man-to-snail FOI

  #State variables #######
    S=n[1]
    E=n[2]
    I=n[3]
    W_TC=n[4]
    W_UC=n[5]
    W_TA=n[6]
    W_UA=n[7]

    #Total snail population
      N=S+E+I

    #Update clumping parameter, k from estimate of worm burden in each population
      k_TC = k_from_log_W(W_TC)
      k_UC = k_from_log_W(W_UC)
      k_TA = k_from_log_W(W_TA)
      k_UA = k_from_log_W(W_UA)

    #Estimate mating probability within each strata
      phi_W_TC = phi_Wk(W = W_TC, k = k_TC)  #Mating probability in treated SAC population
      phi_W_UC = phi_Wk(W = W_UC, k = k_UC)  #Mating probability in untreated SAC population
      phi_W_TA = phi_Wk(W = W_TA, k = k_TA)  #Mating probability in treated adult population
      phi_W_UA = phi_Wk(W = W_UA, k = k_UA)  #Mating probability in untreated adult population


    #Estimate mean eggs produced per person in each strata as product of
      # mated female worms (0.5*worm burden*mating probability),
      # eggs produced per female worm per 10mL urine,
      # urine produced per day /10mL
      # and reduction in fecundity due to crowding
        eggs_W_TC = 0.5*(W_TC*phi_W_TC) * m * rho_Wk(W_TC, zeta, k_TC) * U_C

        eggs_W_UC = 0.5*(W_UC*phi_W_UC) * m * rho_Wk(W_UC, zeta, k_UC) * U_C

        eggs_W_TA = 0.5*(W_TA*phi_W_TA) * m * rho_Wk(W_TA, zeta, k_TA) * U_A

        eggs_W_UA = 0.5*(W_UA*phi_W_UA) * m * rho_Wk(W_UA, zeta, k_UA) * U_A

    #Estimate miracidia produced by each strata as product of
      # mean eggs produced by individuals in each strata,
      # number of people in each strata,
      # egg viability
      # contamination coefficient for SAC/adults,

      M_W_TC = eggs_W_TC * v * H_TC

      M_W_UC = eggs_W_UC * v * H_UC

      M_W_TA = eggs_W_TA * v * H_TA / omega

      M_W_UA = eggs_W_UA * v * H_UA / omega

    # Estimate total miracidia entering snail habitat
      M_tot = M_W_TC + M_W_UC + M_W_TA + M_W_UA

    # Snail infection dynamics
      dSdt= r*(1-(N/K))*(S+E) - mu_N*S - Lambda_0*(1-exp(-beta*M_tot/N))*S #Susceptible snails

      dEdt= Lambda_0*(1-exp(-beta*M_tot/N))*S - (mu_N+sigma)*E #Exposed snails

      dIdt= sigma*E - mu_I*I #Infected snails

    # Estimate cercarial output from infected snail population
      C = theta*I

    #worm burden in human populations; worm acquisition a product of
      # contamination coefficient
      # probability infection per exposure
      # cercarial density/exposure
      # density dependent probability of establishment

      dW_TCdt= omega*alpha*C*gam_Wxi(W_TC, xi) - (mu_W+mu_H_C)*W_TC
      dW_UCdt= omega*alpha*C*gam_Wxi(W_UC, xi) - (mu_W+mu_H_C)*W_UC
      dW_TAdt= alpha*C*gam_Wxi(W_TA, xi) - (mu_W+mu_H_A)*W_TA
      dW_UAdt= alpha*C*gam_Wxi(W_UA, xi) - (mu_W+mu_H_A)*W_UA

    return(list(c(dSdt,dEdt,dIdt,
                  dW_TCdt,dW_UCdt,
                  dW_TAdt, dW_UAdt)))
}

#' Estimate exposed snail fraction from infected snail prevalence
#'
#' Function to estimate proportion of snail population in exposed, E, compartment from model parameters
#' and infected snail, I, prevalence
#'
#' @param I infected snail prevalence
#' @param pars parameter set
#'
#' @return Estimate of porportion of snail population in E compartment
#' @export

I_get_E <- function(I, pars){
  mu_I = pars["mu_I"]
  sigma = pars["sigma"]

  E = (I*mu_I)/sigma

  return(E)
}

#' Estimate peak snail infection from input snail infection prevalence
#'
#' Peak infection rate can be estimated as a function of observed infected snail prevalence as in
#' www.pnas.org/cgi/doi/10.1073/pnas.1708729114. This is a function to estimate this parameter
#' from input infected snail prevalence and model parameters. Depends on function `I_get_E` as well
#'
#' @param I infected snail prevalence
#' @param pars parameter set
#'
#' @return Estimate of porportion of snail population in E compartment
#' @export

I_get_lambda0 <- function(I, pars){
  mu_N = pars["mu_N"]
  mu_I = pars["mu_I"]
  E = I_get_E(I, pars)

  lambda0 = (mu_N*E+mu_I*I)/(1-E-I)

  return(lambda0)
}

#' Model to feed to `rootSolve::multiroot()` to estimate key parameters of the age stratified model based on observed inputs
#'
#' We have four unknown parameters and four solvable equations. This function
#' takes inputs and puts them into equations to estimate the value of these four unkown parameters,
#' alpha, beta, omega, and N_eq, from inputs: adult infection (prevalence and intensity),
#' child infection (prevalence and intensity), child and adult population sizes, infected snail prevalence
#'
#' @param x input estimates of parameters beta, N_eq, alpha, and omega (IN THAT ORDER)
#' @param I infected snail prevalence
#' @param W_A estimated mean worm burden in adult population
#' @param prev_A estimated prevalence in adult population
#' @param H_A adult population size
#' @param W_C estimated mean worm burden in child population
#' @param prev_C estimated prevalence in child population
#' @param H_C child population size
#' @param parameters parameter set of other key parameters
#'
#' @return Vector of four equation solutions
#' @export

age_strat_fit_mod <- function(x, I, W_A, prev_A, H_A, W_C, prev_C, H_C, parameters){

  Lambda_0 <- I_get_lambda0(I, parameters)
  kap_A <- prev_W_get_k(W_A, prev_A)
  kap_C <- prev_W_get_k(W_C, prev_C)

  ##standard snail parameters
    r=parameters["r"]             # recruitment rate (from sokolow et al)
    K=(H_A + H_C)*0.75          # carrying capacity as a function of human population size (assumption that snail area is proportional to community size)
    mu_N=parameters["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=parameters["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=parameters["mu_I"]          # Increased mortality rate of infected snails
    theta=parameters["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = parameters["mu_W"]   # death rate of adult worms
    mu_H_A = parameters["mu_H_A"] # death rate of adult humans
    mu_H_C = parameters["mu_H_C"] # death rate of children
    m = parameters["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = parameters["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = parameters["zeta"]       # parameter of fecundity reduction function
    xi = parameters["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    U_C = parameters["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = parameters["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4

  E_C <- 0.5*W_C*phi_Wk(W_C, kap_C)*m*rho_Wk(W_C,zeta,kap_C)*U_C*v*H_C

    F1 <- I*((mu_I*(mu_N+sigma))*(1+(Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_N+sigma)+(sigma*Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_I*(mu_N+sigma))))*(1/(sigma*Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))) - 1

    F2 <- x[3]*theta*gam_Wxi(W_A, xi)*((sigma*Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2])))*x[2])/(mu_I*(mu_N+sigma)*(1+(Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_N+sigma)+(sigma*Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_I*(mu_N+sigma))))) - (mu_H_A+mu_W)*W_A

    F3 <- x[3]*x[4]*theta*gam_Wxi(W_C, xi)*((sigma*Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2])))*x[2])/(mu_I*(mu_N+sigma)*(1+(Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_N+sigma)+(sigma*Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_I*(mu_N+sigma))))) - (mu_H_A+mu_W)*W_C

    F4 <- r*(1-x[2]/K)*(1+(Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))/(mu_N+sigma))-(mu_N + Lambda_0*(1-exp(-(x[1]*(0.5*W_A*phi_Wk(W_A, kap_A)*m*rho_Wk(W_A,zeta,kap_A)*U_A*v*H_A*(1/x[4])+E_C)/x[2]))))

  c(F1 = F1, F2 = F2, F3 = F3, F4 = F4)
}

#' Function taking `age_strat_fit_mod` and implementing within `rootSolve::multiroot()`
#' to return estimates of key parameters of the age stratified model based on observed inputs
#'
#' We have four unknown parameters and four solvable equations. This function
#' takes inputs and puts them into equations to estimate the value of these four unkown parameters,
#' alpha, beta, omega, and N_eq, from inputs: adult infection (prevalence and intensity),
#' child infection (prevalence and intensity), child and adult population sizes, infected snail prevalence
#'
#' @param f function with equations to feed to `multiroot` function. Defaults to `age_strat_fit_mod`
#' @param start input estimates of parameters beta, N_eq, alpha, and omega (IN THAT ORDER)
#' @param I infected snail prevalence
#' @param W_A estimated mean worm burden in adult population
#' @param prev_A estimated prevalence in adult population
#' @param H_A adult population size
#' @param W_C estimated mean worm burden in child population
#' @param prev_C estimated prevalence in child population
#' @param H_C child population size
#' @param parameters parameter set of other key parameters
#' @param tolerance minimum error to tolerate, else return error
#'
#' @return Estimate of porportion of snail population in E compartment
#' @export

age_strat_get_pars <- function(f = age_strat_fit_mod,
                               start = c(1e-2, 1e3, 1e-6, 4),
                               I, W_A, prev_A, H_A, W_C, prev_C, H_C, parameters,
                               tolerance = 1e-6){
  fit_ob <- multiroot(f = age_strat_fit_mod,
                      start = start,
                      I = I,
                      W_A = W_A,
                      prev_A = prev_A,
                      H_A = H_A,
                      W_C = W_C,
                      prev_C = prev_C,
                      H_C = H_C,
                      parameters = parameters)

  if(fit_ob$estim.precis > tolerance){
    stop("Model fit does not reach your tolerance threshold. Try different starting parameter estimates, lower your tolerance, or something is wrong")
  } else {
    fit_pars <- fit_ob$root
    names(fit_pars) <- c("beta", "N_eq", "alpha", "omega")

    fit_pars

  }

}

#' Augment parameter set for age stratified schisto model with fit transmission and observed demographic/infection parameters
#'
#' Takes as input a few observed variables, uses `age_strat_get_pars` to estimate transmission parameters, adds them to the
#' base parameters, and returns the resulting full parameter set
#'
#' @param I infected snail prevalence
#' @param W_A estimated mean worm burden in adult population
#' @param prev_A estimated prevalence in adult population
#' @param H_A adult population size
#' @param cvrg_A MDA coverage (proportion) among the adult population
#' @param W_C estimated mean worm burden in child population
#' @param prev_C estimated prevalence in child population
#' @param H_C child population size
#' @param cvrg_C MDA coverage (proportion) among the child population
#' @param parameters parameter set of other key parameters
#'
#' @return parameter sit with all parameters needed to run the age stratified model
#' @export
#'
augment_age_strat_parameters <- function(I, W_A, prev_A, H_A, cvrg_A, W_C, prev_C, H_C, cvrg_C, parameters){
  #Make copy to add new pars to
  new_pars <- parameters

  #Get transmisison parameters from fit
  fit_pars <- age_strat_get_pars(I = I,
                                 W_A = W_A,
                                 prev_A = prev_A,
                                 H_A = H_A,
                                 W_C = W_C,
                                 prev_C = prev_C,
                                 H_C = H_C,
                                 parameters = parameters)

  #Add parameters to new parameter set
  new_pars["Lambda_0"] <- as.numeric(I_get_lambda0(I, parameters))
  new_pars["alpha"] <- as.numeric(fit_pars["alpha"])
  new_pars["beta"] <- as.numeric(fit_pars["beta"])
  new_pars["omega"] <- as.numeric(fit_pars["omega"])
  new_pars["K"] <- as.numeric(H_A + H_C * 0.8)  # Snail environmental carrying capacity (related to area, habitat suitability, etc.) a function of human population size
  new_pars["H_TC"] <- round(H_C*cvrg_C)         # Total number of treated children
  new_pars["H_UC"] <- H_C - round(H_C*cvrg_C)   # Total number of untreated children
  new_pars["H_TA"] <- round(H_A*cvrg_A)         # Total number of treated adults
  new_pars["H_UA"] <- H_A - round(H_A*cvrg_A)   # Total number of untreated adults

}

#' Stochastic version of the age stratified schistosomiasis model
#'
#' A stochastic schistosomiasis model with same structure as the age stratified schisto model
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

  } else {
    return(as.data.frame(adaptivetau::ssa.adaptivetau(nstart, transitions, sfx, params, tf)))
  }

}
