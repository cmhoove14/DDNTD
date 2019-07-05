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

    I = n[1]
    W = n[2]

    dIdt = (p$beta_e*W + p$beta_d*I)*(1-I) - p$gamma*I
    dWdt = p$omega + p$V*p$sigma*p$lambda*I - p$rho*W

    return(list(c(dIdt, dWdt)))
}
