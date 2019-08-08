#' Estimate effective reproduction number as product of W, I_P and parameters
#'
#' Takes mean worm burden from unstratified population, infected snail prevalence
#' and additional parameters and returns both R_eff and R_0 estimates
#'
#' @param W mean worm burden
#' @param I_P infected snail prevalence
#' @param pars additional parmaters
#'
#' @return vector with Reff and R0 estimates
#' @export
#'
#'
R_eff_W_I_P <- function(W, I_P, pars){
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
      ((mu_W+mu_H_A)*(mu_I*(mu_N+sigma)))

    R_eff <- 0.5*alpha*theta*beta*H*m*v*I_P*omega_c^2*phi_Wk(W, k_from_log_W(W))*rho_Wk(W, zeta, k_from_log_W(W))*U_C /
      ((mu_W+mu_H_A)*(mu_I*(mu_N+sigma)))

    return(c(as.numeric(R_eff), as.numeric(R_0)))
}

#' Get estimate of effective reproduction number, $R_eff$ absen influence of positive density dependence
#'
#' Estimates the effective reproduction number, $R_eff$ as a function of
#' model parameters, population mean worm burden, and clumping parameter of negative binomial worm burden distribution
#' in the absence of density dependence expressed as the mating probability
#'
#' @param parameters Model parameters
#' @param W mean worm burden in human population
#' @param kap clumping parameter of the worm burden among human population
#'
#' @return estimate of the effective reproduction number, $R_eff$
#' @export

getReff_noPDD<-function(parameters, W, kap){
  with(as.list(parameters),{

  Num_1<-lambda*gam_Wxi(W, xi)*(sigma/(mu_N+mu_I))*(beta*0.5*W*H*rho_Wk(W, zeta, kap)/(mu_N+sigma))*C

  Num_2<- ( (1+(beta*0.5*W*H*rho_Wk(W, zeta, kap)/(mu_N+sigma)))*f_N ) - ( mu_N + (beta*0.5*H*W*rho_Wk(W, zeta, kap) ) )

  Den<-( 1+(beta*0.5*W*H*rho_Wk(W, zeta, kap)/(mu_N+sigma))+((sigma/(mu_N+mu_I))*(beta*0.5*W*H*rho_Wk(W, zeta, kap)/(mu_N+sigma))) ) * ( 1+(beta*0.5*W*H*rho_Wk(W, zeta, kap)/(mu_N+sigma)) ) * ( mu_W+mu_H )* W* f_N

  Reff<-as.numeric(Num_1*Num_2/Den)

  return(Reff)
  })
}
