#' Estimate R_0 from model parameters
#'
#'
#'
#' @param pars parameters
#'
#' @return estimate of R0
#' @export
#'
get_R0 <- function(pars){
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

    R_0 = 0.5*alpha*theta*beta*H*m*v*I_P*omega_c^2*U_C /
      ((mu_W+mu_H_C)*(mu_I*(mu_N+sigma)))

}

#' Get Reff as function of W, model parameters, and DDs
#'
#'
#'
#' @param W mean worm burden
#' @param pars other model parameters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of Reff
#' @export

Reff_W <- function(W, pars,
                   PDD = phi_Wk, DDF = rho_Wk, DDI = gam_Wxi){

  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

  # Get miracidial density as function of worm burden
    #Update clumping parameter
      kap = k_from_log_W(W)

    #Estimate mating probability
      phi = PDD(W = W, k = kap)

    #Estimate density dependent fecundity
      rho = DDF(W, zeta, kap)

    # Estimate total miracidia entering snail habitat
      M_tot = 0.5*H*omega*v*m*W*phi*rho*U

  # Get man-to-snail FOI and snail population sizeas solution given M_tot and other parameters
      #Estimate N
        N <- uniroot.all(function(N) K*(1-(mu_N+(beta*M_tot/N))/(r*(1+(beta*M_tot/(N*(mu_N+sigma))))))-N,
                         interval = c(0,K))

        Lambda <- beta*M_tot/N

  #Density dependent acquired immunity
    gam = DDI(W, xi)

    I_P <- sigma/(((mu_I*(mu_N+sigma))/Lambda)+mu_I+sigma)

  # Reff
    Reff <- (alpha*gam*omega*theta*N*I_P)/
               (W*(mu_W+mu_H))


  return(c("Reff" = as.numeric(Reff),
           "I_P" = as.numeric(I_P),
           "Lambda" = as.numeric(Lambda)))
}

#' Estimate effective reproduction number as product of W, I_P and parameters
#'
#' Takes mean worm burden from unstratified population, infected snail prevalence
#' and additional parameters and returns both R_eff and R_0 estimates
#'
#' @param W mean worm burden
#' @param I_P infected snail prevalence
#' @param pars additional parmaters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return vector with Reff and R0 estimates
#' @export
#'
#'

R_eff_W_I_P <- function(W, I_P, pars,
                        PDD = phi_Wk, DDF = rho_Wk, DDI = gam_Wxi){
  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

    R_0 = alpha*omega^2*theta*beta*H*m*v*U*sigma /
      (2*(mu_W+mu_H)*(mu_I*(mu_N+sigma)))

    R_eff <- alpha*DDI(W, xi)*omega^2*theta*beta*PDD(W, k_from_log_W(W))*DDF(W, zeta, k_from_log_W(W))*H*m*v*U*I_P*(sigma/I_P-mu_I-sigma) /
      (2*(mu_W+mu_H)*(mu_I*(mu_N+sigma)))

    return(c(as.numeric(R_eff), as.numeric(R_0)))
}


#' Get rate of change of mean worm burden as function of W, model parameters, and DDs
#'
#'
#'
#' @param W mean worm burden
#' @param pars other model parameters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of Reff
#' @export

dWdt_W <- function(W, pars,
                   PDD = phi_Wk, DDF = rho_Wk, DDI = gam_Wxi){

  ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

  # Get miracidial density as function of worm burden
    #Update clumping parameter
      kap = k_from_log_W(W)

    #Estimate mating probability
      phi = PDD(W = W, k = kap)

    #Estimate density dependent fecundity
      rho = DDF(W, zeta, kap)

    # Estimate total miracidia entering snail habitat
      M_tot = 0.5*H*omega*v*m*W*phi*rho*U

  # Get man-to-snail FOI and snail population sizeas solution given M_tot and other parameters
      #Estimate N
        N <- uniroot.all(function(N) K*(1-(mu_N+(beta*M_tot/N))/(r*(1+(beta*M_tot/(N*(mu_N+sigma))))))-N,
                         interval = c(0,K))

        Lambda <- beta*M_tot/N

  #Density dependent acquired immunity
    gam = DDI(W, xi)

  # Reff
    dWdt <- (alpha*gam*omega*theta*sigma*N)/(((mu_I*(mu_N+sigma)/Lambda)+mu_I+sigma))-W*(mu_W+mu_H)


  return(as.numeric(dWdt))
}

#' Estimate breakpoint as function of I_P and parameters
#'
#'
#'
#' @param I_P infected snail prevalence
#' @param pars additional parmaters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return vector with Reff and R0 estimates
#' @export
#'
#'

W_bp_I_P <- function(I_P, pars, PDD, DDF, DDI){
    ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

    R_0 = alpha*omega^2*theta*beta*H*m*v*U*sigma /
      (2*(mu_W+mu_H)*(mu_I*(mu_N+sigma)))


  #Couldn't get uniroot to work with density dependence functions, so this is a bit of a hacky way to find the breakpoint
    #based on determining when the product of the DDs is equal to the constants
    if(R_0 < 1){
      W_bp = NA
      print("R_0 is less than 1, there is no breakpoint!")
    } else {
      test_Ws <- exp_seq(1e-6, 100, 10000)
      DDs_product <- sapply(test_Ws, function(W) DDI(W, xi)*PDD(W, k_from_log_W(W))*DDF(W, zeta, k_from_log_W(W)))

      constants_product <- (2*(mu_W+mu_H)*(mu_I*(mu_N+sigma))) /
        (alpha*omega^2*theta*beta*H*m*v*U*I_P*(sigma/I_P-mu_I-sigma))

    if(length(which(abs(DDs_product - constants_product) <= 1e-4) < 5)){
      W_eqs <- test_Ws[which(abs(DDs_product - constants_product) <= 1e-4)]

      W_bp = min(W_eqs)
    } else {
      W_eqs <- test_Ws[which(abs(DDs_product - constants_product) <= 1e-3)]

      W_bp = min(W_eqs)

    }

  }

    return(W_bp)
}

#' Estimate breakpoint as function of I_P and parameters
#'
#'
#'
#' @param I_P infected snail prevalence
#' @param pars additional parmaters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return vector with Reff and R0 estimates
#' @export
#'
#'

W_bp <- function(pars, PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1){
    ##standard snail parameters
    r=pars["r"]             # recruitment rate (from sokolow et al)
    K=pars["K"]          # carrying capacity corresponding to 50 snails per square meter
    mu_N=pars["mu_N"]          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=pars["sigma"]         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=pars["mu_I"]          # Increased mortality rate of infected snails
    theta=pars["theta"]          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = pars["mu_W"]   # death rate of adult worms
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = pars["zeta"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = pars["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    Lambda_0=pars["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

  #Couldn't get uniroot to work with density dependence functions, so this is a bit of a hacky way to find the breakpoint
    #based on determining when the product of the DDs is equal to the constants
      test_Ws <- exp_seq(1e-6, 30, 10000)

      Reff_Ws <- sapply(test_Ws, Reff_W, pars = pars,
                        PDD = PDD, DDF = DDF, DDI = DDI)[1,]
    if(max(Reff_Ws) < 1){
      print("R_0 less than 1, no breakpoint exists")
      return(NA)
    } else if(length(which(abs(Reff_Ws - 1) <= 1e-4) < 5)){
      W_eqs <- test_Ws[which(abs(Reff_Ws - 1) <= 1e-4)]

      W_bp = min(W_eqs)
    } else {
      W_eqs <- test_Ws[which(abs(Reff_Ws - 1) <= 1e-3)]

      W_bp = min(W_eqs)

    }

    return(W_bp)
}



#' Estimate W_bp from R_0
#'
#'
#'
#' @param R_0 basic reproduction number
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#' @param zeta DDF parameter
#' @param DDI parameter
#'
#' @return estimate of the breakpoint mean worm burden
#' @export
#'

W_bp_R_0 <- function(R_0, PDD, DDF, DDI, zeta, xi){
  test_Ws <- exp_seq(1e-6, 20, 10000)

  candidates <- sapply(test_Ws,
                       function(W) R_0*W*PDD(W, k_from_log_W(W))*DDF(W, zeta, k_from_log_W(W))*DDI(W, xi)-1)

  W_bp = min(test_Ws[which(abs(0-candidates) <= 1e-3)])

  return(as.numeric(W_bp))
}

#' Estimate MDA coverage necessary to reach breakpoint as a function of infected snail prevalence and parameters
#'
#'
#' @param W_t mean worm burden prior to MDA
#' @param I_P infected snail prevalence
#' @param epsilon MDA efficacy (proportion of worms cleared from those treated)
#' @param pars additional parmaters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return vector with Reff and R0 estimates
#' @export
#'
#'

cvrg_from_W_bp<- function(W_t, W_bp, epsilon){

    cvrg <- ((W_bp/W_t)-1)/-epsilon

    return(cvrg)

  }

#' Estimate MDA coverage necessary to reach breakpoint as a function of parameters
#' used to esimate breakpoint, density dependencies,
#'
#' @param W_t mean worm burden prior to MDA in treated population
#' @param W_u mean worm burden in untreated population
#' @param epsilon MDA efficacy (proportion of worms cleared from those treated)
#' @param pars additional parameters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of coverage required to reach breakpoint
#' @export
#'
#'

cvrg_from_pars<- function(W_t, W_u, epsilon = 0.94, K, omega, pars,
                          PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1){

  pars["K"] <- K
  pars["omega"] <- omega

    Wbp <- DDNTD::W_bp(pars = pars, PDD = PDD, DDF = DDF, DDI = DDI)

    cvrg <- (Wbp - W_u)/(W_t-epsilon*W_t-W_u)

    return(cvrg)

  }
