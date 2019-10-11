#' Schistosomiasis mating probability
#'
#' Estimate the mating probability from the population mean worm burden and the
#' clumping parameter of the negative binomial distribution
#'
#' @param W mean population worm burden
#' @param k clumping parameter of the negative binomial distribution
#'
#' @return probability a female worm is succefully mated (ranges from 0-1)
#' @export

phi_Wk <- function(W, k) {
  if(W == 0){
    return(0)
  } else {
    b <- ((1-(W/(W + k)))^(1+k))/(2*pi)
    i = integrate(f = function(x, W, k){(1-cos(x))/((1 + (W/(W + k))*cos(x))^(1+k))},
    lower = 0,
    upper = 2*pi,
    stop.on.error = FALSE,
    W = W, k = k)$value

    return(1-b*i)
  }
}

#' Schistosomiasis density dependent fecundity
#'
#' Reductions in egg output due to crowding of adult female worms at high worm burdens
#'
#' @param W Mean worm burden in human population
#' @param zeta shape parameter regulating response of egg ouptut to worm burden
#' @param k clumping parameter of the negative binomial distribution
#'
#' @return Proportional reduction in female worm egg output due to crowding
#' @export

  rho_Wk <- function(W, zeta, k) {
    if(k == 0){
      rho = 1
    } else {
      rho = (1 + (1-exp(-zeta))*(W/k))^-(k+1)
    }
    return(rho)
  }


#' Schistosomiasis acquired immunity
#'
#' Reductions in transmission at high worm burdens due to acquired immunity in human population
#'
#' @param W mean worm burden in human population
#' @param xi density dependence shape parameter regulating response of immunity effect to mean worm burden
#'
#' @return reduction in snail-to-man infection as a result of acquired immunity in human population
#' @export

gam_Wxi <- function(W,xi){
    exp(1-xi*W-exp(-xi*W))
  }

#' Estimate prevalence as function of clumping parameter and mean infection intensity
#'
#' Prevalence, mean intensity and the clumping parameter are all related,
#' therefore estimation of 1 can be achieved if the other two are known
#'
#' @param W Mean worm burden or infection intensity
#' @param k clumping parameter of the negative binomial distribution
#'
#' @return Estimate of the prevalence
#' @export

  est_prev_W_k <- function(W, k){
    1-(1+W/k)^-k
  }

#' Estimate clumping parameter, kappa, as function of mean infection intensity and prevalence using uniroot
#'
#' Prevalence, mean intensity and the clumping parameter are all related,
#' therefore estimation of 1 can be achieved if the other two are known
#'
#' @param W Mean worm burden or infection intensity
#' @param prev prevalence of infection
#'
#' @return Estimate of the clumping parameter
#' @export

prev_W_get_k <- function(W, prev){
  uniroot(function(k){ 1-(1+W/k)^-k-prev},
          interval = c(0,10))$root
}

#' Estimate clumping parameter as function of prevalence and intensity
#'
#' Prevalence, mean intensity and the clumping parameter are all related,
#' therefore estimation of 1 can be achieved if the other two are known
#'
#' @param W Mean worm burden or infection intensity
#' @param Prev Prevalence in th epopulation
#'
#' @return Estimate of the clumping parameter
#' @export

  prev_intensity_to_clump <- function(W, Prev){
    rootSolve::uniroot.all(function(x) 1-Prev-(1+W/x)^-x,
                           interval = c(0,10))
  }

#' Estimate clumping parameter as a function of log worm burden
#'
#' Function fit to our data from 16 communities in Senegal with annual MDA
#' and egg burden assessed annually for three years
#'
#' @param W Mean worm burden or infection intensity
#'
#' @return Estimate of the clumping parameter of the negative binomial distribution
#' @export

  k_w_fx <- function(intensity){
    DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*intensity))
  }

#' Estimate mean worm burden and clumping parameter as function of prevalence
#'
#' This is based on the assumption that the clumping parameter is determined by mean worm burden
#' which reduces the number of unkowns to 2 therefore if prevalence is known, mean worm burden
#' and kappa as determined by the mean worm burden are known
#'
#' @param Prev Prevalence in th epopulation
#'
#' @return Estimate of the mean worm burden
#' @export

  prev_get_W <- function(Prev){
    rootSolve::uniroot.all(function(w) 1-Prev-(1+w/(exp(DDNTD::base_pars["a"] + DDNTD::base_pars["b"]*log(w))))^-(exp(DDNTD::base_pars["a"] + DDNTD::base_pars["b"]*log(w))),
                           interval = c(0,1000))
  }


#' Generate seasonal model forcing function
#'
#' Function to generate a seasonal forcing function based on the simulation time frame
#' the parameter affected, and the phase of the seasonality.
#' time_frame/freq should be a whole number indicating the number of seasonal cycles
#'
#' @param burn_in time period in the beginning to leave parameter unaffected if desired
#' @param t_int time period during which intervention occurs
#' @param t_max total time frame of simulation as numeric
#' @param freq frequency of the seasonality
#' @param par_set named vector of the parameter set that contains the parameter affected
#' @param par character of the parameter in the parameter set that is affected
#' @param par_change relative amplitude of the seasonal change in the parameter value
#'
#' @return Forcing function
#' @export
#'
gen_force_fx <- function(burn_in, t_int, t_max, freq, par_set, par, par_change){
  force_fx <- approxfun(as.data.frame(list(
    times = c(1:t_max),
    par_val = c(rep(par_set[par], times = burn_in),
                sapply(rep(seq(0,pi, length.out = freq), times = (t_int - burn_in)/freq),
                     function(t) par_set[par]*par_change + par_set[par]*(1-par_change) * cos(t+pi)),
                rep(par_set[par], times = (t_max - t_int)))
  )), rule = 2)

  return(force_fx)
}

#' Return 1
#'
#' Function to take any input and return 1, used as a placeholder in some functions where input density dependent functions are required
#'
#' @return 1
#' @export
#'
nil_1 <- function(...){
  return(1)
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

#' Function to estimate Lambda (snail-to-man FOI) as function of snail infection prevalence and snail parameters
#'
#' Snail FOI can be estimated as a function of observed infected snail prevalence and additional snail parameters
#' This is a function to estimate this parameter from input infected snail prevalence and model parameters.
#'
#' @param I_P infected snail prevalence
#' @param mu_N snail daily mortality rate
#' @param mu_I infected snail daily mortality rate
#' @param sigma inverse of snail pre-patency period
#'
#' @return Estimate of man-to-snail FOI
#' @export

I_get_Lambda <- function(I_P, mu_N, mu_I, sigma){
  as.numeric((mu_I*(mu_N+sigma))/((sigma/I_P)-mu_I-sigma))
}

#' Function to estimate equilibirum snail population size, N_eq, as function of Lambda and snail parameters
#'
#' Total equilibrium snail population can be estimated as a function of Lambda, the man-to-snail FOI and other parameters.
#' This function takes as input Lambda (which can be estimated as a function of infected snail prevalence using function `I_get_Lambda`)
#' and additional snail parameters to estimate N_eq
#'
#' @param I_P infected snail prevalence
#' @param K snail population environmental carrying capacity
#' @param mu_N snail daily mortality rate
#' @param r snail intrinsic reproduction rate
#' @param sigma inverse of snail pre-patency period
#'
#' @return Estimate of man-to-snail FOI
#' @export

Lambda_get_N_eq <- function(Lambda, K, mu_N, r, sigma){
  as.numeric(K*(1-(mu_N+Lambda)/(r*(1+(Lambda/(mu_N+sigma))))))
}

#' Function to estimate beta from egg output and other parameters
#'
#'
#'
#' @param egg_output mean eggs produced per 10mL per person (observed diagnostic data)
#' @param H number of individuals in the human population
#' @param Lambda snail force of infection
#' @param N_eq snail population size
#' @param U mean urine produced per day per individual in units of 10mL
#' @param v mean egg viability (probability an egg produced successfully hatches into miracidia)
#' @param omega human exposure/contamination parameter that regulates proportion of eggs produced that interact with snail population
#'
#' @return estimate of the per miracidial contact probability of infection, beta
#' @export
#'

beta_from_eggs <- function(egg_output, H, Lambda, N_eq, U, v, omega){
  as.numeric((Lambda*N_eq)/(egg_output*H*U*v*omega))
}

#' Function used in `multiroot`` to fit transmission parameters
#'
#' Fits transmission parameters, alpha and beta, based on equilibrium solutions to snail and worm burden equations
#' Given unkown beta, alpha, and equilibirum snail population, N_eq
#'
#' @param x starting values of unkowns, a vector of beta, alpha, and N_eq
#' @param parms additional model parameters
#' @param W equilibirum mean worm burden
#' @param Ip equlibirum infected snail prevalence
#'
#' @return vector of solutions to three equations from inputs
#' @export
#'

W_Ip_eq_solns <- function(x, parms, W, Ip){
  #standard snail parameters
    r=parms["r"]             # recruitment rate (from sokolow et al)
    K=parms["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=parms["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=parms["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=parms["mu_I"]          # Increased mortality rate of infected snails
    theta=parms["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = parms["mu_W"]   # death rate of adult worms
    mu_H = parms["mu_H"] # death rate of adult humans
    m = parms["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = parms["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = parms["zeta"]       # parameter of fecundity reduction function
    xi = parms["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = parms["H"]
    U = parms["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = parms["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = parms["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  M_tot = 0.5*W*H*omega*U*m*v*phi_Wk(W,k_w_fx(W))*rho_Wk(W,parms["zeta"],k_w_fx(W))

  # Equilibrium snail population size from unkown beta (x[1]) and snail population size (x[3])
  F1 <- K*(1-(mu_N+(x[1]*M_tot/x[3]))/(r*(1+(x[1]*M_tot/(x[3]*(mu_N+sigma))))))-x[3]

  # Snail force of infection (Lambda) as function of infected snail prevalence
  F2 <- (x[3]*mu_I*(mu_N+sigma))/((sigma/Ip-mu_I-sigma)*x[1]*M_tot) - 1

  #Equlibirum worm burden from unkown beta (x[1]), alpha (x[2]) and snail population size (x[3])
  F3 <- (x[2]*omega*theta*sigma*x[3])/((mu_W+mu_H)*((x[3]*mu_I*(mu_N+sigma))/(x[1]*M_tot)+mu_I+sigma)) - W

  return(c(F1 = F1, F2 = F2, F3 = F3))
}


#' Uses `W_Ip_eq_solns` within `multiroot` to solve for unkown transmission parameters
#'
#' Uses euqlibirum equations for snail population, mean worm burden, and snail FOI
#' and input equilibirum mean worm burden and infected snail prevalence to solve for unkown transmission
#' parameters and infected snail prevalence
#'
#' @param beta_guess initial estimate for snail FOI parameter, beta
#' @param alpha_guess initial estimate for human FOI parameter, alpha
#' @param Neq_guess initial estimate for equlibirum snail population size
#' @param W input equlibirum mean worm burden
#' @param Ip input equlibirum infected snail prevalence
#' @param pars additional model parameters
#'
#' @return vector of alpha, beta, and N_eq estimates
#' @export

fit_pars_from_eq_vals <- function(beta_guess, alpha_guess, Neq_guess, W, Ip, pars){

  multiroot(W_Ip_eq_solns, start = c(beta_guess, alpha_guess, Neq_guess),
            parms = pars, positive = TRUE,
            W = W, Ip = Ip)$root

}

#' Function used in multiroot to estimate worm burden and clumping parameter from egg burden, prevalence and parameters
#'
#' Function which uses equations representing 1) an estimate of prevalence given worm burden and clumping parameter and
#' 2) estimated mean egg output per 10mL urine given worm burden, mating probability, DD fecundity, and peak egg release per 10mL
#' to estimate the worm burden and clumping parameter from input egg burden and prevalence
#'
#' @param x initial estimates for the mean worm burden and clumping parameter, respectively
#' @param parms input egg burden (measured from surveys) measured in eggs/10mL for the population group
#'
#' @return Value of the two equations given inputs, but is mostly irrelevant other than its use within multiroot (see function `convert_burden_egg_to_worm`)
#' @export
#'

egg_to_worm_fx <- function(x, parms){

  egg_burden <- parms[1]
  prev <- parms[2]
  m <- parms[3]
  zeta <- parms[4]

  mate_integral <- function(t){
      (1-cos(t))/((1 + (x[1]/(x[1] + x[2]))*cos(t))^(1+x[2]))
  }

  mate_prob <- ifelse(x[2] == 0, 1,
                      (1-integrate(mate_integral, 0, 2*pi)$value*((1-(x[1]/(x[1] + x[2])))^(1+x[2]))/(2*pi)))


# From equation relating mean egg burden to mean worm burden and clumping parameter
  F1 <- 0.5*x[1]*(mate_prob*m*(1 + (1-exp(-zeta))*(x[1]/x[2]))^(-x[2]-1))- # Estimate of eggs produced per 10mL urine
    egg_burden

# From equation relating prevalence of mated pairs to mean worm burden and clumping parameter
  F2 <- 1 - 2*(1+x[1]/(2*x[2]))^(-x[2]) + (1+x[1]/x[2])^(-x[2])-prev

  return(c(F1 = F1, F2 = F2))
}

#' Function to get solutions for worm burden and clumping parameter from input egg burden and prevalence
#'
#' Uses `egg_to_worm_fx` within `rootSolve::multiroot` equation solver to come up with estimates of mean worm burden
#' and clumping parameter from input egg burden, prevalence and egg output parameters
#'
#' @param egg_burden observed mean egg burden in population fraction
#' @param prevalence observed prevalence in population fraction
#' @param m peak egg output per 10mL per mated female worm
#' @param zeta density dependent fecundity parameter
#'
#' @return estimate of the mean worm burden and clumping parameter
#' @export
#'

convert_burden_egg_to_worm <- function(egg_burden, prevalence, m, zeta){
  multiroot(egg_to_worm_fx, start = c(egg_burden/2, 0.001), # Use egg burden/2 and prevalence as guesses for starting W and kap estimates
            parms = c(egg_burden, prevalence, m, zeta), positive = TRUE)$root
}

#' Function to estimate worm burden given mean egg burden and dispersion of egg burden
#'
#' Uses `uniroot` to solve for W assuming mean egg burden = 0.5Wphi(W,k)rho(W,k)m
#'
#' @param eggs mean egg output
#' @param kap disperion parameter on egg output
#' @param zeta fecundity reduction parameter
#' @param m mean egg output per mated female worm with no fecundity reduction
#'
#' @return estimate of the mean worm burden
#' @export
#'

eggs_kap_get_W <- function(eggs, kap, zeta, m){
  uniroot(function(W){
    mate_integral <- integrate(f = function(x, W, kap){(1-cos(x))/((1 + (W/(W + kap))*cos(x))^(1+kap))},
                               lower = 0,
                               upper = 2*pi,
                               stop.on.error = FALSE,
                               W = W, k = kap)$value
    0.5*W*m*((1 + (1-exp(-zeta))*(W/kap))^-(kap+1))*(1-((1-(W/(W + kap)))^(1+kap))/(2*pi)*mate_integral)-eggs
  }, interval = c(0,1000))$root
}

#' Cheever eggs to females function
#'
#' Uses reported log-log relationship between eggs in urine and female worms to estimate worm burden from egg burden
#'
#' @param eggs observed egg burden in eggs/10mL urine
#'
#' @return estimate of fecund female worms from egg burden
#' @export
#'

cheever_eggs_to_females <- function(eggs){
  exp((log(1+eggs)-0.75)/1.467) - 1
}

#' Estimate R0 from starting worm burden and density dependence parameters
#'
#'
#'
#' @param w_start endemic equilibrium mean worm burden
#' @param kap endemic equilibrium dispersion parameter
#' @param zeta negative densiy dependent fecundity parameter
#'
#' @return estimate of basic reproduction number, R0
#' @export

w_start_get_r0 <- function(w_start, kap, zeta){
  1/(phi_Wk(w_start, kap)*rho_Wk(w_start, zeta, kap))
}
