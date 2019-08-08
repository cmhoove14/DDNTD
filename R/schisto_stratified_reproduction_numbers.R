#' Estimate effective reproduction number
#'
#' Estimates the effective reproduction number, $R_eff$ as a function of
#' model parameters and mean worm burden in each age and treatment population
#'
#' @param pars parameter set
#' @param W_TC mean worm burden in treated SAC group
#' @param W_UC mean worm burden in untreated SAC group
#' @param W_TA mean worm burden in treated adult group
#' @param W_UA mean worm burden in untreated adult group
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param S_FOI form of man-to-snail FOI either "linear" or "saturating"
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of te effective reproduction number in each population class, the population-weighted average, and the population average mean worm burden
#' @export

Reff_Wij <- function(pars, W_TC, W_UC, W_TA, W_UA,
                     PDD, DDF, S_FOI, DDI){
  # Get mean worm burden of total population as weighted sum of worm burden in each age/treatment group
  W_bar <- W_TC*pars["h_tc"]+
           W_UC*pars["h_uc"]+
           W_TA*pars["h_ta"]+
           W_UA*pars["h_ua"]

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

  N_eq = pars["N_eq"]

  # Get miracidial density as function of worm burdens
    #Update clumping parameter, k from estimate of worm burden in each population
      k_TC = k_from_log_W(W_TC)
      k_UC = k_from_log_W(W_UC)
      k_TA = k_from_log_W(W_TA)
      k_UA = k_from_log_W(W_UA)

    #Estimate mating probability within each strata
      phi_W_TC = PDD(W = W_TC, k = k_TC)  #Mating probability in treated SAC population
      phi_W_UC = PDD(W = W_UC, k = k_UC)  #Mating probability in untreated SAC population
      phi_W_TA = PDD(W = W_TA, k = k_TA)  #Mating probability in treated adult population
      phi_W_UA = PDD(W = W_UA, k = k_UA)  #Mating probability in untreated adult population

    #Estimate density dependent fecundity
      rho_W_TC = DDF(W_TC, zeta, k_TC)
      rho_W_UC = DDF(W_UC, zeta, k_UC)
      rho_W_TA = DDF(W_TA, zeta, k_TA)
      rho_W_UA = DDF(W_UA, zeta, k_UA)

    # Estimate total miracidia entering snail habitat
      M_tot = 0.5*H*omega_a*v*m*((W_TC*phi_W_TC) * rho_W_TC * U_C*h_tc*Omega +
                                   (W_UC*phi_W_UC) * rho_W_UC * U_C*h_uc*Omega +
                                   (W_TA*phi_W_TA) * rho_W_TA * U_A*h_ta +
                                   (W_UA*phi_W_UA) * rho_W_UA * U_A*h_ua)

  # Get man-to-snail FOI as solution given M_tot and other parameters and choice of specification
    if(S_FOI == "linear"){

      #Estimate N
        N <- uniroot.all(function(N) r*(1-N/K)*(1+(beta*M_tot)/(N*(mu_N+sigma)))-mu_N-beta*M_tot/N,
                         interval = c(0,K))

      Lambda <- beta*M_tot/N

    } else if(S_FOI == "saturating"){

      N <- uniroot.all(function(N) r*(1-N/K)*(1+(Lambda_0*(1-exp(-M_tot/N)))/(mu_N+sigma))-mu_N-Lambda_0*(1-exp(-M_tot/N)),
                       interval = c(0,K))

      Lambda <- Lambda_0*(1-exp(-M_tot/N))

    } else if(!(S_FOI %in% c("linear", "saturating"))){

      stop("S_FOI must be either linear or saturating")

    }

  #Density dependent acquired immunity
    gam_W_TC = DDI(W_TC, xi)
    gam_W_UC = DDI(W_UC, xi)
    gam_W_TA = DDI(W_TA, xi)
    gam_W_UA = DDI(W_UA, xi)

  # Net Reff
    sum_bit <- h_tc*omega_c^2*phi_W_TC*rho_W_TC*gam_W_TC*U_C/(mu_W+mu_H_C) +
               h_uc*omega_c^2*phi_W_UC*rho_W_UC*gam_W_UC*U_C/(mu_W+mu_H_C) +
               h_ta*omega_a^2*phi_W_TA*rho_W_TA*gam_W_TA*U_A/(mu_W+mu_H_A) +
               h_ua*omega_a^2*phi_W_UA*rho_W_UA*gam_W_UA*U_A/(mu_W+mu_H_A)

    Reff <- (0.5*alpha*theta*beta*H*m*v*(sigma/(mu_I*(mu_N+sigma)/Lambda+mu_I+sigma))/(2*mu_I*(mu_N+sigma)))*sum_bit

  return(c("W_bar" = as.numeric(W_bar),
           "Reff" = as.numeric(Reff),
           "Reff_W_TC" = Reff_W_TC,
           "Reff_W_UC" = Reff_W_UC,
           "Reff_W_TA" = Reff_W_TA,
           "Reff_W_UA" = Reff_W_UA))
}


#' Estimate effective reproduction number as function of mean worm burden in different strata, infected snail prevalence and other parameters
#'
#' Estimates the effective reproduction number, $R_eff$ as a function of
#' model parameters, mean worm burden in each age and treatment population, infected snail prevalence, and
#' different options for density dependence parameters
#'
#' @param pars parameter set
#' @param W_TC mean worm burden in treated SAC group
#' @param W_UC mean worm burden in untreated SAC group
#' @param W_TA mean worm burden in treated adult group
#' @param W_UA mean worm burden in untreated adult group
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of the effective reproduction number and the population average mean worm burden
#' @export

Reff_Wij_I_P <- function(pars, W_TC, W_UC, W_TA, W_UA, I_P,
                         PDD, DDF, DDI){
  # Get mean worm burden of total population as weighted sum of worm burden in each age/treatment group
  W_bar <- W_TC*pars["h_tc"]+
           W_UC*pars["h_uc"]+
           W_TA*pars["h_ta"]+
           W_UA*pars["h_ua"]

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
      phi_W_TC = PDD(W = W_TC, k = k_TC)  #Mating probability in treated SAC population
      phi_W_UC = PDD(W = W_UC, k = k_UC)  #Mating probability in untreated SAC population
      phi_W_TA = PDD(W = W_TA, k = k_TA)  #Mating probability in treated adult population
      phi_W_UA = PDD(W = W_UA, k = k_UA)  #Mating probability in untreated adult population

    #Estimate density dependent fecundity
      rho_W_TC = DDF(W_TC, zeta, k_TC)
      rho_W_UC = DDF(W_UC, zeta, k_UC)
      rho_W_TA = DDF(W_TA, zeta, k_TA)
      rho_W_UA = DDF(W_UA, zeta, k_UA)

  #Density dependent acquired immunity
    gam_W_TC = DDI(W_TC, xi)
    gam_W_UC = DDI(W_UC, xi)
    gam_W_TA = DDI(W_TA, xi)
    gam_W_UA = DDI(W_UA, xi)

  # Net Reff
    sum_bit <- h_tc*omega_c^2*phi_W_TC*rho_W_TC*gam_W_TC*U_C/(mu_W+mu_H_C) +
               h_uc*omega_c^2*phi_W_UC*rho_W_UC*gam_W_UC*U_C/(mu_W+mu_H_C) +
               h_ta*omega_a^2*phi_W_TA*rho_W_TA*gam_W_TA*U_A/(mu_W+mu_H_A) +
               h_ua*omega_a^2*phi_W_UA*rho_W_UA*gam_W_UA*U_A/(mu_W+mu_H_A)

    Reff <- (0.5*alpha*theta*beta*H*m*v*I_P/(2*mu_I*(mu_N+sigma)))*sum_bit


  return(c("W_bar" = as.numeric(W_bar),
           "Reff" = as.numeric(Reff)))
}


#' Model to feed to `rootSolve::multiroot()` to estimate breakpoint from input parameters
#'
#' Uses the equilibrium snail population and reff equations to estimate equilibirum snail population size
#' and population breakpoint. This function takes estimates of these two values and returns values of these equations.
#' Best use is within `multiroot` to find roots which represent the snail population sizr and mean worm burden at the breakpoint
#'
#' @param x input estimates of values N_eq and W_bp (IN THAT ORDER)
#' @param parms parameter set of other key parameters
#'
#' @return Vector of two solutions
#' @export

W_bp_N_solver <- function(x, parms){
  ##standard snail parameters
    r=parms["r"]             # recruitment rate (from sokolow et al)
    K=parms["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=parms["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=parms["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=parms["mu_I"]          # Increased mortality rate of infected snails
    theta=parms["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = parms["mu_W"]   # death rate of adult worms
    mu_H_A = parms["mu_H_A"] # death rate of adult humans
    mu_H_C = parms["mu_H_C"] # death rate of children
    m = parms["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = parms["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = parms["zeta"]       # parameter of fecundity reduction function
    xi = parms["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = parms["H"]
    h_tc = parms["h_tc"]         # Proportion of treated children
    h_uc = parms["h_uc"]          # Proportion of untreated children
    h_ta = parms["h_ta"]           # Proportion of treated adults
    h_ua = parms["h_ua"]         # Proportion of untreated adults
    U_C = parms["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = parms["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = parms["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = parms["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = parms["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    alpha=parms["alpha"]       # Cercarial infection probability
    Lambda_0=parms["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=parms["beta"]       # Man to snail trnamission probability for linear FOI

    #Lambda equation
      F1 <- as.numeric(r*(1-x[1]/K)*(1+(beta*(0.5*x[2]*H*phi_Wk(x[2], k_from_log_W(x[2]))*rho_Wk(x[2], zeta, k_from_log_W(x[2]))*m*v*U_C*omega_c))/(x[1]*(mu_N+sigma)))-mu_N-beta*(0.5*x[2]*H*phi_Wk(x[2], k_from_log_W(x[2]))*rho_Wk(x[2], zeta, k_from_log_W(x[2]))*m*v*U_C*omega_c)/x[1])

    #W_bp equation
      F2 <- as.numeric((alpha*omega_c*theta*sigma*(0.5*H*phi_Wk(x[2], k_from_log_W(x[2]))*rho_Wk(x[2], zeta, k_from_log_W(x[2]))*m*v*U_C*omega_c))/
                         ((mu_I*(mu_N+sigma))*(1+(beta*(0.5*x[2]*H*phi_Wk(x[2], k_from_log_W(x[2]))*rho_Wk(x[2], zeta, k_from_log_W(x[2]))*m*v*U_C*omega_c))/(x[1]*(mu_N+sigma))+(sigma*beta*(0.5*x[2]*H*phi_Wk(x[2], k_from_log_W(x[2]))*rho_Wk(x[2], zeta, k_from_log_W(x[2]))*m*v*U_C*omega_c))/(x[1]*mu_I*(mu_N+sigma)))) - (mu_W+mu_H_C))

  #print(c(x[1], x[2], F1, F2))

  return(c(F1 = F1, F2 = F2))

}

#' Get breakpoint snail infection prevalence and mean worm burden
#'
#' Uses `multiroot` to get estimates of snail population size and
#' breakpoint mean worm burden from equations for snail population dynamics and worm burden estimation
#' then uses these estimates to determine breakpoint snail prevalence
#'
#' @param N_guess initial guess for breakpoint snail population size
#' @param W_bp_guess initial guess for breakpoint mean worm burden
#' @param parms parameter values
#'
#'

get_Ip_W_bp <- function(N_guess, W_bp_guess, parms){
  soln <- multiroot(f = W_bp_N_solver,
                    start = c(N_guess, W_bp_guess),
                    parms = parms,
                    positive = TRUE)

  N_eq = soln$root[1]
  W_bp = soln$root[2]

    ##standard snail parameters
    mu_N=parms["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=parms["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=parms["mu_I"]          # Increased mortality rate of infected snails

  #Adult Worm, Miracidia and Cercariae Parameters
    m = parms["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = parms["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = parms["zeta"]       # parameter of fecundity reduction function
    xi = parms["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = parms["H"]
    h_tc = parms["h_tc"]         # Proportion of treated children
    h_uc = parms["h_uc"]          # Proportion of untreated children
    h_ta = parms["h_ta"]           # Proportion of treated adults
    h_ua = parms["h_ua"]         # Proportion of untreated adults
    U_C = parms["U_C"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = parms["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega_c = parms["omega_c"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    omega_a = parms["omega_a"]          #  infection risk/contamination of adults (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4
    Omega = parms["Omega"]          # relative infection risk/contamination of SAC vs adults

  #Transmission parameters
    Lambda_0=parms["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=parms["beta"]       # Man to snail trnamission probability for linear FOI

  I_P <- sigma/((mu_I*(mu_N+sigma))/(beta*0.5*W_bp*phi_Wk(W_bp, k_from_log_W(W_bp))*rho_Wk(W_bp, zeta, k_from_log_W(W_bp))*m*v*U_C*omega_c)/N_eq+mu_I+sigma)

  return(c(I_P = as.numeric(I_P),
           W_bp = as.numeric(W_bp)))
}
