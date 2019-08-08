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
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param parameters Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_age_strat_mod <- function(t, n, pars) {
  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    h_tc = pars["h_tc"]         # Total number of treated children
    h_uc = pars["h_uc"]          # Total number of untreated children
    h_ta = pars["h_ta"]           # Total number of treated adults
    h_ua = pars["h_ua"]         # Total number of untreated adults
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = pars["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = pars["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI

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

    #Estimate total miracidia produced by each strata as product of
      # mean eggs produced by individuals in each strata,
      # number of people in each strata,
      # egg viability
      # contamination coefficient for SAC/adults,

    # Estimate total miracidia entering snail habitat
      M_tot = H*omega_a*v*(eggs_W_TC*h_tc*Omega + eggs_W_UC*h_uc*Omega + eggs_W_TA*h_ta + eggs_W_UA*h_ua)

    # Snail infection dynamics
      dSdt= r*(1-(N/K))*(S+E) - mu_N*S - Lambda_0*(1-exp(-M_tot/N))*S #Susceptible snails

      dEdt= Lambda_0*(1-exp(-M_tot/N))*S - (mu_N+sigma)*E #Exposed snails

      dIdt= sigma*E - mu_I*I #Infected snails

    # Estimate cercarial output from infected snail population
      C = theta*I

    #worm burden in human populations; worm acquisition a product of
      # contamination coefficient
      # probability infection per exposure
      # cercarial density/exposure
      # density dependent probability of establishment

      dW_TCdt= omega_c*alpha*C - (mu_W+mu_H_C)*W_TC
      dW_UCdt= omega_c*alpha*C - (mu_W+mu_H_C)*W_UC
      dW_TAdt= omega_a*alpha*C - (mu_W+mu_H_A)*W_TA
      dW_UAdt= omega_a*alpha*C - (mu_W+mu_H_A)*W_UA

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


  F1 <- (2*egg_burden)/ #2* egg burden for 1:1 sex ratio
    ((1-integrate(mate_integral, 0, 2*pi)$value*((1-(x[1]/(x[1] + x[2])))^(1+x[2]))/(2*pi))* #mating probability
    m* #egg release per 10mL urine
    (1 + ((x[1]*(1-(exp(-zeta))))/x[2]))^(-x[2]-1))- #DD fecundity
    x[1]

  F2 <- 1-(1+x[1]/x[2])^-x[2]-prev

  #print(c(x[1], x[2], integrate(mate_integral, 0, 2*pi)$value))

  #print(c(x[1], x[2], F1, F2, (1 + ((x[1]*(1-(exp(-zeta))))/x[2]))^(-x[2]-1), (1-integrate(mate_integral, 0, 2*pi)$value*((1-(x[1]/(x[1] + x[2])))^(1+x[2]))/(2*pi))))

  return(c(F1 = F1, F2 = F2))
}

#' Function to get solutions for worm burden and clumping parameter from input egg burden and prevalence
#'
#' Uses `egg_to_worm_fx` within `rootSolve::multiroot` equation solver to come up with estimates of mean worm burden
#' and clumping parameter from input egg burden, prevalence and egg output parameters
#'
#' @param W_guess initial estimate fed to multiroot for mean worm burden
#' @param kap_guess initial estimate of clumping parameter fed to multiroot
#' @param egg_burden observed mean egg burden in population fraction
#' @param prevalence observed prevalence in population fraction
#' @param m peak egg output per 10mL per mated female worm
#' @param zeta density dependent fecundity parameter
#'
#' @return estimate of the mean worm burden and clumping parameter
#' @export
#'

convert_burden_egg_to_worm <- function(W_guess, kap_guess, egg_burden, prevalence, m, zeta){
  multiroot(egg_to_worm_fx, start = c(W_guess, kap_guess), parms = c(egg_burden, prevalence, m, zeta), positive = TRUE)$root
}

#' Function to estimate miracidial invasion rate and relative contamination coefficient from equilibirum egg burden and prevalence inputs
#'
#' Function which takes input infection and demographic information and estimates relative infection between adults and children
#' assuming equilibrium differences in worm burden are a result of this parameter. Then estimates peak miracidial invasion rate
#' from a resulting estimate of the equilibrium miracidial density as a function of worm burden, contamination coefficients and other parameters
#'
#' @param W_A equilibrium mean worm burden in adult population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param kap_A equilibrium clumping parameter in adult population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param H_A adult human population size
#' @param W_C equilibrium mean worm burden in child population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param kap_C equilibrium clumping parameter in child population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param H_C child human population size
#' @param K_ratio ratio of snail carrying capacity to human population size. defaults to 1, but should be increased for environments/communities that have particularly intense snail habitat or tranmission
#' @param I_P infected snail prevalence, observed or input
#' @param pars additional model parameters
#'
#' @return Estimates of the relative child to adult exposure/contamination ratio
#' @export
#'

eq_Ws_get_Omega_Lambda0 <- function(W_A,
                                    kap_A,
                                    H_A,
                                    W_C,
                                    kap_C,
                                    H_C,
                                    K_ratio = 1,
                                    I_P,
                                    pars){
  #parameters as functions of inputs
    H = H_A + H_C
    h_a = H_A/H
    h_c = H_C/H

  # Snail parameters
    r=pars["r"]           # recruitment rate (from sokolow et al)
    K=H * K_ratio         # carrying capacity as a function of human population size (assumption that snail area is proportional to community size)
    mu_N=pars["mu_N"]     # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]   # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]     # Increased mortality rate of infected snails
    theta=pars["theta"]   # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]  # amount of infectious material reaching environment

  # Parameters determined from inputs
    Lambda <- I_get_Lambda(I_P, mu_N, mu_I, sigma)
    N_eq <- Lambda_get_N_eq(Lambda, K, mu_N,r, sigma)
    Omega = W_C/W_A
    omega_a = omega_c/Omega

  #Miracidal concentration from inputs
    M = 0.5*H*omega_a*m*v*(W_C*phi_Wk(W_C, kap_C)*rho_Wk(W_C,zeta,kap_C)*U_C*h_c*Omega + W_A*phi_Wk(W_A, kap_A)*rho_Wk(W_A,zeta,kap_A)*U_A*h_a)

  #Lambda_0 estimate for non-linear snail FOI
    Lambda_0 <- Lambda/(1-exp(-M/N_eq))

  # Beta estimate for linear snail FOI
    beta <- Lambda*N_eq/M

    return(c(as.numeric(Omega), as.numeric(Lambda_0), as.numeric(beta)))

}

#' Function to estimate probability of worm establishment per cercarial exposure from equilibrium worm burden
#'
#' Assuming observed worm burden prior to any intervention is at euilibirum, we can estimate the worm acquisition rate
#' as a function of this equilibirum worm burden and host and worm turnover (mortality) rates. The only component
#' of the worm acquisition rate that is unkown is the probability of worm establishment per cercarial exposure, alpha, which we can estimate
#' from equilibirum mean worm burden and known parameters
#'
#' @param W_C equilibrium mean worm burden in the child population
#' @param omega_c contamination/exposure coefficient for children
#' @param mu_W daily mortality rate of adult worms
#' @param mu_H_C mortality rate of children
#' @param theta daily cercarial shedding rate of infected snails
#' @param I_P infected snail prevalence
#' @param N_eq equilibrium snail population size

eq_W_get_alpha <- function(W_C, omega_c, mu_W, mu_H_C, theta, I_P, N_eq){
   as.numeric((W_C*(mu_W + mu_H_C))/(omega_c*I_P*N_eq*theta))
}

#' Function to input observed equilibrium infection values and return fully "fit" parameter set based on observed infection values
#'
#' Uses `eq_W_get_alpha`, `eq_Ws_get_Omega_Lambda0`, `convert_burden_egg_to_worm`, `Lambda_get_N_eq`, and `I_get_Lambda` to convert input
#' infection values (mean egg burden and infection prevalence in adult and child populations, infected snail prevalence) into estimates of model parameters
#'
#' @param W_A equilibrium mean worm burden in adult population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param kap_A equilibrium clumping parameter in adult population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param H_A adult human population size
#' @param cvrg_A MDA coverage of the adult population
#' @param W_C equilibrium mean worm burden in child population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param kap_C equilibrium clumping parameter in child population, can be estimated from egg burden and prevalence using `convert_burden_egg_to_worm`
#' @param H_C child human population size
#' @param cvrg_C MDA coverage of the child population
#' @param K_ratio ratio of snail carrying capacity to human population size. defaults to 1, but should be increased for environments/communities that have particularly intense snail habitat or tranmission
#' @param I_P infected snail prevalence, observed or input
#' @param pars additional model parameters
#'
#' @return parameter set with fit parameters based on equilibirum values
#' @export

infection_inputs_get_pars <- function(W_A, kap_A, H_A, cvrg_A, W_C, kap_C, H_C, cvrg_C, K_ratio = 1, I_P, pars){
  # Snail parameters
    r=pars["r"]           # recruitment rate (from sokolow et al)
    mu_N=pars["mu_N"]     # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]   # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]     # Increased mortality rate of infected snails
    theta=pars["theta"]   # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]  # amount of infectious material reaching environment

  #Make copy to fill with new parameters
  new_pars <- pars

  #parameters as functions of inputs
    H = H_A+H_C
    h_a = H_A/H
    h_c = H_C/H
    h_ta = round(H_A*cvrg_A) / H
    h_ua = h_a - h_ta
    h_tc = round(H_C*cvrg_C) / H
    h_uc = h_c - h_tc

    K = H * K_ratio

    new_pars["H"] = H_A + H_C
    new_pars["h_a"] = h_a
    new_pars["h_c"] = h_c
    new_pars["h_t"] = (cvrg_C*H_C + cvrg_A*H_A)/H
    new_pars["h_u"] = ((1-cvrg_C)*H_C + (1-cvrg_A)*H_A)/H
    new_pars["h_tc"] = h_tc
    new_pars["h_uc"] = h_uc
    new_pars["h_ta"] = h_ta
    new_pars["h_ua"] = h_ua
    new_pars["cvrg_C"] = cvrg_C
    new_pars["cvrg_A"] = cvrg_A

    new_pars["K"] = K

  # Parameters determined from inputs
    Lambda = I_get_Lambda(I_P, mu_N, mu_I, sigma)
    N_eq = Lambda_get_N_eq(Lambda, K, mu_N, r, sigma)
    Omega = W_C/W_A
    omega_a = omega_c/Omega

    new_pars["Lambda"] <- Lambda
    new_pars["N_eq"] <- N_eq
    new_pars["I_eq"] <- new_pars["N_eq"]*I_P
    new_pars["E_eq"] <- new_pars["I_eq"]*mu_I/sigma
    new_pars["S_eq"] <- (new_pars["I_eq"]*mu_I*(mu_N+sigma))/(sigma*Lambda)
    new_pars["Omega"] = Omega
    new_pars["omega_a"] = omega_c/Omega
    new_pars["alpha"] = eq_W_get_alpha(W_C = W_C,
                                       omega_c = new_pars["omega_c"],
                                       mu_W = new_pars["mu_W"],
                                       mu_H_C = new_pars["mu_H_C"],
                                       theta = new_pars["theta"],
                                       I_P = I_P,
                                       N_eq = new_pars["N_eq"])

  #Miricidia concentration from inputs
      M = 0.5*H*omega_a*m*v*
          (W_C*phi_Wk(W_C, kap_C)*rho_Wk(W_C,zeta,kap_C)*U_C*h_c*Omega +
           W_A*phi_Wk(W_A, kap_A)*rho_Wk(W_A,zeta,kap_A)*U_A*h_a)

      new_pars["Lambda_0"] <- Lambda/(1-exp(-M/N_eq))
      new_pars["beta"] <- Lambda*N_eq/M

    return(new_pars)

}


#' Estimate saturating form of man-to-snail FOI
#'
#' Estimates the man-to-snail FOI as function of
#' model parameters and mean worm burden in each age and treatment population
#'
#' @param pars parameter set
#' @param W_TC mean worm burden in treated SAC group
#' @param W_UC mean worm burden in untreated SAC group
#' @param W_TA mean worm burden in treated adult group
#' @param W_UA mean worm burden in untreated adult group
#'
#' @return Estimate of the man-to-snail FOI from the saturating function
#' @export

Lambda_Wij <- function(pars, W_TC, W_UC, W_TA, W_UA){
  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    h_tc = pars["h_tc"]         # Proportion of treated children
    h_uc = pars["h_uc"]          # Proportion of untreated children
    h_ta = pars["h_ta"]           # Proportion of treated adults
    h_ua = pars["h_ua"]         # Proportion of untreated adults
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = pars["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = pars["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI

  # Get miracidial density as function of worm burdens
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

    # Estimate total miracidia entering snail habitat
      M_tot = 0.5*H*omega_a*v*m*((W_TC*phi_W_TC) * rho_Wk(W_TC, zeta, k_TC) * U_C*h_tc*Omega +
                                   (W_UC*phi_W_UC) * rho_Wk(W_UC, zeta, k_UC) * U_C*h_uc*Omega +
                                   (W_TA*phi_W_TA) * rho_Wk(W_TA, zeta, k_TA) * U_A*h_ta +
                                   (W_UA*phi_W_UA) * rho_Wk(W_UA, zeta, k_UA) * U_A*h_ua)

  # Get man-to-snail FOI as solution given M_tot and other parameters
    Lambda <- uniroot(function(L) Lambda_0*(1-exp(-M_tot/(K*(1-(mu_N+L)/(r*(1+L/(mu_N+sigma)))))))-L,
                      interval = c(1e-8,10))$root

  return(Lambda)
}

#' Estimate linear form of man-to-snail FOI
#'
#' Estimates the linear version of man-to-snail FOI as function of
#' model parameters and mean worm burden in each age and treatment population
#'
#' @param pars parameter set
#' @param W_TC mean worm burden in treated SAC group
#' @param W_UC mean worm burden in untreated SAC group
#' @param W_TA mean worm burden in treated adult group
#' @param W_UA mean worm burden in untreated adult group
#'
#' @return Estimate of the man-to-snail FOI from the linear function
#' @export

Lambda_Wij_linear <- function(pars, W_TC, W_UC, W_TA, W_UA){
  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    h_tc = pars["h_tc"]         # Proportion of treated children
    h_uc = pars["h_uc"]          # Proportion of untreated children
    h_ta = pars["h_ta"]           # Proportion of treated adults
    h_ua = pars["h_ua"]         # Proportion of untreated adults
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = pars["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = pars["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI

  # Get miracidial density as function of worm burdens
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

    # Estimate total miracidia entering snail habitat
      M_tot = 0.5*H*omega_a*v*m*((W_TC*phi_W_TC) * rho_Wk(W_TC, zeta, k_TC) * U_C*h_tc*Omega +
                                   (W_UC*phi_W_UC) * rho_Wk(W_UC, zeta, k_UC) * U_C*h_uc*Omega +
                                   (W_TA*phi_W_TA) * rho_Wk(W_TA, zeta, k_TA) * U_A*h_ta +
                                   (W_UA*phi_W_UA) * rho_Wk(W_UA, zeta, k_UA) * U_A*h_ua)

  # Get man-to-snail FOI as solution given M_tot and other parameters
      Lambda <- uniroot(function(L) beta*M/(K*(1-(mu_N+L)/(r*(1+L/(mu_N+sigma)))))-L,
                        interval = c(0, 10))$root

  return(Lambda)
}

#' Estimate endemic equilibrium given infected snail prevalence
#'
#' Using Reff expression and fact that Reff=1 at the endemic equilibrium, we can express the endemic equilibirum in terms of
#' the man-to-snail FOI, Lambda. We can also express Lambda in terms of the infected snail prevalence, I_P,
#' therefore we can estimate the endemic equilibrium as a function of model parameters and input snail prevalence
#'
#' @param I_P infected snail prevalence
#' @param pars named vector or list of other parameters
#'
#' @return estimate of the breakpoint
#' @export

I_P_get_W_eq <- function(I_P, pars){
  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    h_tc = pars["h_tc"]         # Proportion of treated children
    h_uc = pars["h_uc"]          # Proportion of untreated children
    h_ta = pars["h_ta"]           # Proportion of treated adults
    h_ua = pars["h_ua"]         # Proportion of untreated adults
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = pars["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = pars["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

    munsig = mu_N+sigma
    sigmui = mu_I+sigma

  W_eq <- (alpha*omega_c*theta*K*I_P*mu_I*munsig*(1+(I_P*mu_I/(sigma-I_P*sigmui))-mu_N/r-I_P*mu_I*munsig/(r*(sigma-I_P*sigmui)))) /
      ((1-I_P*mu_I/sigma-I_P)*(mu_W+mu_H_C)*(mu_I*(mu_N+sigma+I_P*mu_I*munsig/(sigma-I_P*sigmui))+I_P*mu_I*munsig/(1-I_P*mu_I/sigma-I_P))*(1+(I_P*mu_I/(sigma-I_P*sigmui))))

  return(W_eq)
}

#' Estimate Snail infection class distribution given mean worm burden in different infection/treatment classes
#'
#' Takes mean worm burden estimates in different worm burden classes and estimates miraicidial concentration
#' Then uses snail population dynamics equation as a function of LINEAR snail FOI to estimate snail population size
#' Then uses parameters to estimate infection class distribution from total snail population size
#'
#' @param W_TC mean worm burden in treated SAC group
#' @param W_UC mean worm burden in untreated SAC group
#' @param W_TA mean worm burden in treated adult group
#' @param W_UA mean worm burden in untreated adult group
#' @param pars named vector or list of other parameters
#'
#' @return estimate of N, S, E, and I
#' @export

W_ij_get_N_S_E_I <- function(W_TC, W_UC, W_TA, W_UA, pars){
  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    h_tc = pars["h_tc"]         # Proportion of treated children
    h_uc = pars["h_uc"]          # Proportion of untreated children
    h_ta = pars["h_ta"]           # Proportion of treated adults
    h_ua = pars["h_ua"]         # Proportion of untreated adults
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = pars["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = pars["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = pars["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

  # Get miracidial density as function of worm burdens
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

    # Estimate total miracidia entering snail habitat
      M = 0.5*H*omega_a*v*m*((W_TC*phi_W_TC) * rho_Wk(W_TC, zeta, k_TC) * U_C*h_tc*Omega +
                               (W_UC*phi_W_UC) * rho_Wk(W_UC, zeta, k_UC) * U_C*h_uc*Omega +
                               (W_TA*phi_W_TA) * rho_Wk(W_TA, zeta, k_TA) * U_A*h_ta +
                               (W_UA*phi_W_UA) * rho_Wk(W_UA, zeta, k_UA) * U_A*h_ua)

  #Estimate N
    N <- uniroot.all(function(N) r*(1-N/K)*(1+(beta*M)/(N*(mu_N+sigma)))-mu_N-beta*M/N,
                     interval = c(0,K))

    Lambda <- beta*M/N

    S <- as.numeric(N/(1+Lambda/(mu_N+sigma)+(sigma*Lambda)/(mu_I*(mu_N+sigma))))

    E <- as.numeric(S*Lambda/(mu_N+sigma))

    I <- as.numeric(S*Lambda*sigma/(mu_I*(mu_N+sigma)))

    return(c(N, S, E, I))

}











#NEED UPDATING
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



# Deprecated FUNCTIONS #############

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
#'

age_strat_fit_mod <- function(x, I, W_A, prev_A, H_A, W_C, prev_C, H_C, pars){

  Lambda_0 <- I_get_lambda0(I, pars)
  kap_A <- prev_W_get_k(W_A, prev_A)
  kap_C <- prev_W_get_k(W_C, prev_C)

  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=(H_A + H_C)*0.8          # carrying capacity as a function of human population size (assumption that snail area is proportional to community size)
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4

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
#' @param pars parameter set of other key parameters
#' @param tolerance minimum error to tolerate, else return error
#'
#' @return Estimate of porportion of snail population in E compartment
#'

age_strat_get_pars <- function(f = age_strat_fit_mod,
                               start = c(1e-2, 1e3, 1e-6, 4),
                               I, W_A, prev_A, H_A, W_C, prev_C, H_C, pars,
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
                      pars = pars)

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
#' @param pars parameter set of other key parameters
#'
#' @return parameter set with all parameters needed to run the age stratified model
#'
#'

augment_age_strat_parameters <- function(I, W_A, prev_A, H_A, cvrg_A, W_C, prev_C, H_C, cvrg_C, pars){
  #Make copy to add new pars to
  new_pars <- pars

  #Get transmisison parameters from fit
  fit_pars <- age_strat_get_pars(I = I,
                                 W_A = W_A,
                                 prev_A = prev_A,
                                 H_A = H_A,
                                 W_C = W_C,
                                 prev_C = prev_C,
                                 H_C = H_C,
                                 pars = pars)

  #Add parameters to new parameter set
  new_pars["Lambda_0"] <- as.numeric(I_get_lambda0(I, pars))
  new_pars["alpha"] <- as.numeric(fit_pars["alpha"])
  new_pars["beta"] <- as.numeric(fit_pars["beta"])
  new_pars["omega"] <- as.numeric(fit_pars["omega"])
  new_pars["N_eq"] <- as.numeric(fit_pars["N_eq"])
  new_pars["K"] <- as.numeric(H_A + H_C * 0.8)  # Snail environmental carrying capacity (related to area, habitat suitability, etc.) a function of human population size
  new_pars["H_TC"] <- round(H_C*cvrg_C)         # Total number of treated children
  new_pars["H_UC"] <- H_C - round(H_C*cvrg_C)   # Total number of untreated children
  new_pars["H_TA"] <- round(H_A*cvrg_A)         # Total number of treated adults
  new_pars["H_UA"] <- H_A - round(H_A*cvrg_A)   # Total number of untreated adults

  return(new_pars)
}

#' Solve for equilibrium susceptible snail and total snail population
#'
#' Takes as input the worm burden in each age and treatment group and
#' vector of model parameters and uses `rootSolve::multiroot` to determine the
#' total, N, and susceptible, S, snail population sizes. Returns these values along
#' with resulting estimates of exposed, E, and infected, I, snail populations.
#' These values returned as ratios, i.e. S/N; E/N; I/N
#'
#' @param x input estimates of eqbm population sizes, N and S (IN THAT ORDER) for initial estimates
#' @param pars parameter set
#' @param W_TC mean worm burden in treated SAC group
#' @param W_UC mean worm burden in untreated SAC group
#' @param W_TA mean worm burden in treated adult group
#' @param W_UA mean worm burden in untreated adult group
#'
#' @return vector of estimates of total equlibirum snail population size and proportion susceptible, exposed, and infected
#'

eqbm_snails_fit_mod <- function(x, pars, W_TC, W_UC, W_TA, W_UA){
  # Get total human population size
  H_tot <- pars["H_TC"]+pars["H_UC"]+pars["H_TA"]+pars["H_UA"]

  # Get mean worm burden of total population as weighted sum of worm burden in each age/treatment group
  W_bar <- W_TC*pars["H_TC"]/H_tot+
           W_UC*pars["H_UC"]/H_tot+
           W_TA*pars["H_TA"]/H_tot+
           W_UA*pars["H_UA"]/H_tot


    # Get all parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=H_tot*0.75          # carrying capacity as a function of human population size (assumption that snail area is proportional to community size)
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H_A = pars["mu_H_A"] # death rate of adult humans
    mu_H_C = pars["mu_H_C"] # death rate of children
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    U_C = pars["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4

  # parameters as function of Ws
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

     #Miracidial contribution of each group to estimate M
        M_W_TC = 0.5*(W_TC*phi_W_TC) * m * rho_Wk(W_TC, zeta, k_TC) * U_C * v * H_TC

        M_W_UC = 0.5*(W_UC*phi_W_UC) * m * rho_Wk(W_UC, zeta, k_UC) * U_C * v * H_UC

        M_W_TA = 0.5*(W_TA*phi_W_TA) * m * rho_Wk(W_TA, zeta, k_TA) * U_A * v * H_TA / omega

        M_W_UA = 0.5*(W_UA*phi_W_UA) * m * rho_Wk(W_UA, zeta, k_UA) * U_A * v * H_UA / omega

        M_tot = M_W_TC + M_W_UC + M_W_TA + M_W_UA

    N_eq <- uniroot.all(f = function(N){
      r*(1-N/K)*(1+Lambda_0*(1-exp(-beta*M_tot/N)))- (mu_N + Lambda_0*(1-exp(-beta*M_tot/N)))
    },
    interval = c(0,K))

}

