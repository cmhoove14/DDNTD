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
    mu_H = pars["mu_H"] # death rate of adult humans
    m = pars["m"]             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = pars["v"]           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    gamma = pars["gamma"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

  # R_0
    R_0 = (0.5*beta*m*v*U*H*omega^2*sigma*alpha*theta)/((mu_N+sigma)*mu_I*(mu_W+mu_H))

  return(R_0)
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
    gamma = pars["gamma"]       # parameter of fecundity reduction function
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
      kap = k_w_fx(W)

    #Estimate mating probability
      phi = PDD(W = W, k = kap)

    #Estimate density dependent fecundity
      rho = DDF(W, gamma, kap)

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
    Reff <- (alpha*omega*theta*N*I_P)/(W*(mu_H+mu_W))
    #  (K*alpha*omega*theta*sigma*(1+Lambda/(mu_N+sigma)-mu_I/r-Lambda/r)) / (W*(mu_W+mu_H)*(((mu_I*(mu_N+sigma))/Lambda)+mu_I+sigma)*(1+Lambda/(mu_N+sigma)))


  return(c("Reff" = as.numeric(Reff),
           "N" = as.numeric(N),
           "I_P" = as.numeric(I_P),
           "Lambda" = as.numeric(Lambda)))
}

#' Get Reff as function of W, model parameters, kappa, and DDs
#'
#'
#'
#' @param W mean worm burden
#' @param kap clumping parameter
#' @param pars other model parameters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of Reff
#' @export

Reff_W_kap <- function(W, kap, pars,
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
    gamma = pars["gamma"]       # parameter of fecundity reduction function
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
      kap = kap

    #Estimate mating probability
      phi = PDD(W = W, k = kap)

    #Estimate density dependent fecundity
      rho = DDF(W, gamma, kap)

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
    Reff <- (alpha*omega*theta*N*I_P)/(W*(mu_H+mu_W))
    #  (K*alpha*omega*theta*sigma*(1+Lambda/(mu_N+sigma)-mu_I/r-Lambda/r)) / (W*(mu_W+mu_H)*(((mu_I*(mu_N+sigma))/Lambda)+mu_I+sigma)*(1+Lambda/(mu_N+sigma)))


  return(c("Reff" = as.numeric(Reff),
           "N" = as.numeric(N),
           "I_P" = as.numeric(I_P),
           "Lambda" = as.numeric(Lambda)))
}

#' Estimate effective reproduction number as product of W, I and parameters
#'
#' Takes mean worm burden from unstratified population, infected snail population size
#' and additional parameters and returns R_eff estimate
#'
#' @param W mean worm burden
#' @param I infected snail population size
#' @param pars additional parmaters
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return vector with Reff and R0 estimates
#' @export
#'
#'

R_eff_W_I <- function(W, I, pars,
                        PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1){
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
    gamma = pars["gamma"]       # parameter of fecundity reduction function
    xi = pars["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = pars["H"]
    U = pars["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = pars["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=pars["alpha"]       # Cercarial infection probability
    beta=pars["beta"]       # Man to snail trnamission probability for linear FOI

    R_eff <- (alpha*DDI(W, xi)*omega*theta*I) / (W*(mu_W+mu_H))

    return(as.numeric(R_eff))
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
    gamma = pars["gamma"]       # parameter of fecundity reduction function
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
      kap = k_w_fx(W)

    #Estimate mating probability
      phi = PDD(W = W, k = kap)

    #Estimate density dependent fecundity
      rho = DDF(W, gamma, kap)

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
    gamma = pars["gamma"]       # parameter of fecundity reduction function
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
      DDs_product <- sapply(test_Ws, function(W) DDI(W, xi)*PDD(W, k_w_fx(W))*DDF(W, gamma, k_w_fx(W)))

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

#' Worm burden breakpoint equations for use in multiroot
#'
#' Equations for worm burden breakpoint and snail population at the breakpoint used in `multiroot` to identify roots
#' which represent the breakpoint and snail population at the breakpoint
#'
#' @param x vector of the two state variables, Wbp and N(Wbp)
#' @param parms additional parmaters
#'
#' @return vector with equation estimates from Wbp and N(Wbp)
#' @export
#'
#'

wbp_fx <- function(x, parms){
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
    gamma = parms["gamma"]       # parameter of fecundity reduction function
    xi = parms["xi"]        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    H = parms["H"]
    U = parms["U"]          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = parms["U_A"]          # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    omega = parms["omega"]          #  infection risk/contamination of SAC  (related to sanitation/education/water contact) 10.1186/s13071-016-1681-4

  #Transmission parameters
    alpha=parms["alpha"]       # Cercarial infection probability
    Lambda_0=parms["Lambda_0"]         # first parameter of non-linear man-to-snail FOI
    beta=parms["beta"]       # Man to snail trnamission probability for linear FOI

  mate_integral <- function(t){
    (1-cos(t))/((1 + (x[1]/(x[1] + (DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))*cos(t))^(1+(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))
  }


  F1 <- (alpha*omega*theta*sigma*x[2]) /
    ((mu_W+mu_H)*((mu_I*(mu_N+sigma)) /
                    (beta*0.5*H*omega*v*m*x[1]*U*
                       (1-integrate(mate_integral, 0, 2*pi)$value*((1-(x[1]/(x[1] + (DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1]))))))^(1+(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))/(2*pi))* #Phi
                       ((1 + ((x[1]*(1-(exp(-gamma))))/(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))^(-(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))-1))/x[2])+mu_I+sigma)) - x[1]

  F2 <- K*(1-(mu_N + (beta*0.5*H*omega*v*m*x[1]*U*
                       (1-integrate(mate_integral, 0, 2*pi)$value*((1-(x[1]/(x[1] + (DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1]))))))^(1+(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))/(2*pi))* #Phi
                       ((1 + ((x[1]*(1-(exp(-gamma))))/(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))^(-(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))-1))/x[2]))/(r*(1+(beta*0.5*H*omega*v*m*x[1]*U*
                       (1-integrate(mate_integral, 0, 2*pi)$value*((1-(x[1]/(x[1] + (DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1]))))))^(1+(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))/(2*pi))* #Phi
                       ((1 + ((x[1]*(1-(exp(-gamma))))/(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))))^(-(DDNTD::base_pars["a"]*(1-exp(DDNTD::base_pars["b"]*x[1])))-1))/(x[2]*(mu_N+sigma)))))) - x[2]

  return(c(F1 = F1, F2 = F2))

}

#' Worm burden breakpoint estimation
#'
#' Uses `multiroot` and `wbp_fx` functions to estimate Wb and N(Wbp)
#'
#' @param Wbp_guess initial value of Wbp to use
#' @param Nbp_guess initial value of N(Wbp) to use
#' @param pars additional parmaters
#'
#' @return vector with Wbp and N(Wbp)
#' @export
#'
#'

get_breakpoint_W_N <- function(Wbp_guess, Nbp_guess, pars){
  r0_init <- get_R0(pars)

  if(r0_init < 1){
    rtrn <- NA
    print("R0 less than 1, no breakpoint")
  } else {
    rtrn <- multiroot(wbp_fx, start = c(Wbp_guess, Nbp_guess), parms = pars, positive = TRUE)$root
  }
  return(rtrn)
}


#' Estimate W_bp from R_0
#'
#'
#'
#' @param R_0 basic reproduction number
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#' @param gamma DDF parameter
#' @param DDI parameter
#'
#' @return estimate of the breakpoint mean worm burden
#' @export
#'

W_bp_R_0 <- function(R_0, PDD, DDF, DDI, gamma, xi){
  test_Ws <- exp_seq(1e-6, 20, 10000)

  candidates <- sapply(test_Ws,
                       function(W) R_0*W*PDD(W, k_w_fx(W))*DDF(W, gamma, k_w_fx(W))*DDI(W, xi)-1)

  W_bp = min(test_Ws[which(abs(0-candidates) <= 1e-3)])

  return(as.numeric(W_bp))
}

#' Estimate MDA coverage necessary to reach breakpoint as a function of infected snail prevalence and parameters
#'
#'
#' @param W_t mean worm burden prior to MDA
#' @param W_bp breakpoint mean worm burden
#' @param epsilon MDA efficacy (proportion of worms cleared from those treated)
#'
#' @return estimate of MDA coverage necessary to reach breakpoint
#' @export
#'
#'

cvrg_from_W_bp<- function(W_t, W_bp, epsilon){

    cvrg <- (W_bp-W_t)/(-epsilon*W_t)

    return(cvrg)

  }

#' Estimate MDA coverage necessary to reach breakpoint as a function of parameters
#' used to esimate breakpoint, density dependencies,
#'
#' @param use_pars list of parameters necessary to estimate breakpoint and coverage including
#' mean pre-treatment worm burden in treated population, W_t, mean pre-treatment worm burden in the untreated population, W_u, drug efficacy, epsilon, snail population environmental carrying capacity, K, and sanitation/exposure parameter, omega, and other model parameters, pars
#' @param PDD positive density dependence function
#' @param DDF density dependent fecundity function
#' @param DDI density dependent acquired immunity function
#'
#' @return estimate of coverage required to reach breakpoint
#' @export
#'
#'

cvrg_from_pars<- function(use_pars,
                          PDD = phi_Wk, DDF = rho_Wk, DDI = nil_1){
  fin_pars <- use_pars$pars
  W_t = use_pars$W_t
  W_u = use_pars$W_u
  epsilon = use_pars$epsilon

  fin_pars["K"] <- use_pars$K
  fin_pars["omega"] <- use_pars$omega

    Wbp <- DDNTD::W_bp(pars = fin_pars, PDD = PDD, DDF = DDF, DDI = DDI)

    cvrg <- (Wbp - W_u)/(W_t-epsilon*W_t-W_u)

    return(c(Wbp, cvrg))

  }

#' Estimate R_0 from input egg burden and density dependence parameters
#'
#' @param eggs mean population egg burden
#' @param kap dispersion parameter
#' @param gamma negative density dependence parameter
#' @param m peak egg output
#'
#' @return estimate of basic reproduction number, R0
#' @export

R0_egg_kap_gamma_m <- function(eggs, kap, gamma, m){
  W_est <- eggs_kap_get_W(eggs, kap, gamma, m)

  R0 <- 1/(phi_Wk(W = W_est, k = kap)*rho_Wk(W = W_est, gamma = gamma, k = kap))

  return(R0)
}
