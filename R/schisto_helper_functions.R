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
    (1 + ((W*(1-(exp(-zeta))))/k))^(-k-1)
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

  k_from_log_W <- function(W){
    exp(-2.5 + 0.3*log(W))
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
    rootSolve::uniroot.all(function(x) 1-Prev-(1+x/(exp(-2.59748 + 0.41378*log(x))))^-(exp(-2.59748 + 0.41378*log(x))),
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
