#' Garchitorena et al model for environmentally transmitted infectious diseases
#'
#' Deterministic version of the simple, two-differential-equation model
#' presented in Garchitorena et al 2018
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param p Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export
#'
garch_mod = function(t, n, p){
  with(as.list(p),{
    I = n[1]
    W = n[2]

    dIdt = (beta_e*W + beta_d*I)*(1-I) - gamma*I
    dWdt = omega + V*sigma*lambda*I - rho*W

    return(list(c(dIdt, dWdt)))
  })
}

#' Garchitorena et al model for environmentally transmitted infectious diseases
#' augmented with density dependent function
#'
#' Deterministic version of the simple, two-differential-equation model
#' presented in Garchitorena et al 2018 with a function dependent on
#' prevalence (I) that represents a density dependence in transmission
#' and gives rise to a transmission breakpoint and endemic equilibrium
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param p Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export
#'
garch_mod_dd = function(t, n, p){
  with(as.list(p),{
    I = n[1]
    W = n[2]

    dIdt = (beta_e*W + beta_d*I)*(1-I) - gamma*I
    dWdt = omega + V*sigma*lambda*I*(4*I*(1-I)) - rho*W

    return(list(c(dIdt, dWdt)))
  })
}


#' Simulate the garchitorena model
#'
#' Simulate garchitorena deterministic model through time using `ode` function
#' given the model deSolve function, parameter set, time frame, starting conditions, and potential events (e.g. MDA)
#'
#' @param nstart Named vector of starting values for state variables
#' @param time Numeric vector of times at which state variables should be estimated
#' @param model Name of the ode function to use, defaults to `garch_mod`
#' @param parameters Named vector or list of parameter values
#' @param events_df Data frame of events such as MDA with columns "var", "time", "value", and "method"
#'
#' @return dataframe of state variable values at requested times
#' @export
#'
sim_garch_mod <- function(nstart,
                            time,
                            model = garch_mod,
                            parameters,
                            events_df = NA){
  if(is.data.frame(events_df)){
    as.data.frame(ode(nstart, time, model, parameters,
                      events = list(data = events_df)))
  } else {
    as.data.frame(ode(nstart, time, model, parameters))
  }
}


#' Garchitorena et al $R_0$ estimation
#'
#' For the simple, two-equation infectious disease model in
#' `garch_mod`, estimate basic reproduction number, R_0 for the given parameter set
#'
#' @param p parameter set
#'
#' @return estimate of $R_0$
#' @export
#'
garch_r0 <- function(p){
  with(as.list(p),{
  #person to person contribution
  r0d <- beta_d/gamma

  #environmental contribution
  r0e <- (sigma*V*lambda*beta_e) / (gamma*rho)

  r0d + r0e
  })
}

#' Simulate transmission and intervention through time with variable inputs
#'
#' For the simple, two-equation infectious disease model in
#' `garch_mod`, simulate transmission and intervention with choices
#' for transmission intensity (R0), relative cost of environmental to human
#' intervention (script_P), timeframe, capital available, allocation of capital,
#' parameters, and model formulation
#'
#' @param R0 transmission intensity
#' @param script_P relative cost of environmental to human intervention
#' @param T_frame total time frame to simulate
#' @param freq frequency of intervention
#' @param M total capital available per intervention
#' @param A numeric between 0 and 1 of proportion of capital to go towards human intervention (remainder goes towards environmental intervention)
#' @param pars parameter set with remainin transmission parameters
#' @param mod model formulation (e.g. with density dependence `garch_mod_pdd` or not `garch_mod`)
#'
#' @return estimate of utility associated with inputs. Measured as the sum over time of prevalence^1.5 (penalizes higher prevalence values)
#' @export
#'

sim_w_inputs <- function(R0, script_P, T_frame, freq, M, A, pars, mod){
  pars["beta_e"] <- (pars["gamma"] * pars["rho"]) * R0
  pars["mu"] <- pars["theta"]*script_P

  #Get reductions in prevalence based on capital and allocation
  M_A_W <- M*(1-A)/pars["mu"]*0.01
  M_A_I <- M*A/pars["theta"]*0.01

  #Create events dataframes based on capital allocation decisions
  W_events <- data.frame(var=rep('W', times = round(T_frame/freq)),
                         time = c(1:round(T_frame/freq))*freq,
                         value = rep(M_A_W, times = round(T_frame/freq)),
                         method = rep("mult", times = round(T_frame/freq)))

  I_events <- data.frame(var=rep('I', times = round(T_frame/freq)),
                         time = c(1:round(T_frame/freq))*freq,
                         value = rep(M_A_I, times = round(T_frame/freq)),
                         method = rep("mult", times = round(T_frame/freq)))

  eq_vals <- runsteady(y = c(I = 0.5, W = 20), func = mod,
                       parms = pars)[["y"]]

  sim <- as.data.frame(ode(eq_vals, c(1:T_frame), mod, pars,
                           events = list(data = rbind(W_events, I_events) %>% arrange(time))))

  U <- -sum((sim$I*100)^1.5)

  return(U)
}


#' Garchitorena discrete time simulator
#'
#' Same model as in `garch_mod`, but discrete time version for
#' ease of implementation with stochastic dynamic programing framework
#'
#' @param I value of the state variable I at time t
#' @param p parameter set
#'
#' @return value of I at t+1
#' @export
#'
garch_discrete_mod <- function(I, p){

  I + ((p["beta_e"] * p["V"] * p["sigma"] * p["lambda"] * I)/p["rho"])*(1-I) - p["gamma"]*I

}

#' Garchitorena discrete time model simulation
#'
#' Simulates discrete time version of the garchitorena model
#' builds on `garch_discrete_mod`
#'
#' @param I_0 Initial value of prevalence
#' @param time Total time (days) to simulate
#' @param model Name of the discrete sim function to use, defaults to `garch_discrete_mod`
#' @param parameters Named vector or list of parameter values
#' @param events_df Data frame of events such as MDA with columns "var", "time", "value", and "method"
#'
#' @return dataframe of state variable values at requested times
#' @export
#'
sim_discrete_mod <- function(I_0,
                             time,
                             model = garch_discrete_mod,
                             parameters,
                             events_df = NA){
  #initialize vector to fill
  fill <- vector("numeric", length = time)

  #fill first value with starting prevalence
  fill[1] <- I_0

  if(class(events_df) == "data.frame"){
    for(t in 2:time){
      if(t %in% events_df$time){
        fill[t] <- garch_discrete_mod(fill[t-1], parameters) - fill[t-1]*events_df$value[which(t == events_df$time)]
      } else {
        fill[t] <- garch_discrete_mod(fill[t-1], parameters)
      }
    }
  } else {
    for(t in 2:time){
        fill[t] <- garch_discrete_mod(fill[t-1], parameters)
      }
  }

  return(data.frame(t = c(1:time),
                    I = fill))
}

#' Garchitorena discrete time model simulation with intervention variable
#'
#' Simulates discrete time version of the garchitorena model
#' builds on `garch_discrete_mod` and incorporates decision variable A_t
#' that regulates where capital is committed to intervention
#'
#' @param A_t decision variable: proportion of capital to allocate towards MDA
#' @param I_0 Initial value of prevalence
#' @param int_par parameter for remaining capital intervention to be allocated towards
#' @param time Total time (days) to simulate
#' @param parameters Named vector or list of parameter values
#'
#' @return value of prevalence state variable at the end of the time period of intervention
#' @export
#'
sim_int_choice <- function(A_t,
                           I_0,
                           int_par,
                           time = 365,
                           parameters){
  #initialize vector to fill
    fill <- vector("numeric", length = time+2)

  #fill first value with starting prevalence
    fill[1] <- I_0

  #Adjust parameters based on intervention allocation
    use_pars <- parameters
    use_pars[int_par] <- parameters[int_par] * (1-(1-A_t)*(parameters["M"]/parameters["mu"])*0.01)

  #Implement MDA based on capital allocation in the first time step
   fill[2] <- garch_discrete_mod(fill[1], use_pars) - fill[1]*(A_t*use_pars["M"]/use_pars["theta"])*0.01

  #Run the model out for a year
    for(t in 3:(time+2)){
      fill[t] <- garch_discrete_mod(fill[t-1], use_pars)
    }

  return(fill)
}

#' Implement an optimal policy identified with MDPtoolbox through time
#'
#' Simulates optimal treatment allocation based on control variable A_t
#' in the discrete time version of the garchitorena model
#' builds on `garch_discrete_mod` and incorporates decision variable A_t
#' that regulates where capital is committed to intervention
#'
#' @param A_vec vector of possible decisions
#' @param states vector of possible states of the system
#' @param opt_list List from `MDPtoolbox::mdp_finite_horizon`
#' @param p_start index of which state to start at
#' @param t_per_step Total time (days) to simulate
#' @param int_par parameter environmental intervention acts on
#' @param parameters Named vector or list of parameter values
#'
#' @return Data frame of prevalence over course of intervention with optimal intervention implemented at each time step
#' @export
#'
sim_opt_choice <- function(A_vec,
                           states,
                           opt_list,
                           p_start,
                           t_per_step,
                           int_par,
                           parameters){

  #initialize vector to fill
    fill <- list()

  #fill first value with starting prevalence for a year run in
    fill[[1]] <- list(time = c(1:t_per_step),
                      I = rep(states[p_start], t_per_step))

  #Fill subsequent time based on time steps
    for(i in 1:ncol(opt_list[["policy"]])){
    #optimal allocation at time step
      A_t <- A_vec[opt_list[["policy"]][p_start, i]]

  #Adjust parameters based on intervention allocation
    use_pars <- parameters
    use_pars[int_par] <- as.numeric(parameters[int_par] * (1-(1-A_t)*(parameters["M"]/parameters["mu"])*0.01))

    #Implement MDA based on capital allocation in the time step
      fill_next <- numeric()
      fill_next[1] <- fill[[i]][["I"]][t_per_step]
      fill_next[1] <- garch_discrete_mod(fill_next[1], use_pars) - fill_next[1]*(A_t*use_pars["M"]/use_pars["theta"])*0.01

    #Run the model out for a year
      for(t in 2:t_per_step){
        fill_next[t] <- garch_discrete_mod(fill_next[t-1], use_pars)
      }

    fill[[i+1]] <- list(time = c((t_per_step*i+1):(t_per_step*(i+1))),
                        I = fill_next)

  }
  return(as.data.frame(bind_rows(fill)))
}
