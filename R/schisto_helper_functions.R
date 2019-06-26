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
  val = integrate(function(x, W, k){
    a <- ( W / (W + k) )
    b <- ((1-a)^(1+k))/(2*pi)
    return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
  },
  0, 2*pi, W = W, k = k)$value

  #}
    return(1-val)
}

#' Schistosomiasis fecundity reduction due to crowding at high densities
#'
#' Reductions in egg output due to crowding of adult female worms at high worm burdens
#'
#' @param W Mean worm burden in human population
#' @param gamma shape parameter regulating response of egg ouptut to worm burden
#' @param k clumping parameter of the negative binomial distribution
#'
#' @return Proportional reduction in female worm egg output due to crowding
#' @export

  f_Wgk <- function(W, gamma, k) {
    if(k <= 0){
      return(1)
    } else {
      return((1 + ((W*(1-(exp(-gamma))))/k))^(-k-1))
    }

  }


#' Schistosomiasis acquired immunity
#'
#' Reductions in transmission at high worm burdens due to acquired immunity in human population
#'
#' @param W mean worm burden in human population
#' @param v density dependence shape parameter regulating response of immunity effect to mean worm burden
#'
#' @return reduction in snail-to-man infection as a result of acquired immunity in human population
#' @export

  R_Wv <- function(W,v){
    exp(1-v*W-exp(-v*W))
  }

#' Estimate prevalence as function of clumping parameter and mean infection intensity
#'
#'   Prevalence, mean intensity and the clumping parameter are all related,
#'   therefore estimation of 1 can be achieved if the other two are known
#'
#'  @param W Mean worm burden or infection intensity
#'  @param k clumping parameter of the negative binomial distribution
#'
#'  @return Estimate of the prevalence
#'  @export

  est_prev_W_k <- function(W, k){
    1-(1+W/k)^-k
  }



#' Estimate clumping parameter as function of prevalence and intensity
#'
#'   Prevalence, mean intensity and the clumping parameter are all related,
#'   therefore estimation of 1 can be achieved if the other two are known
#'
#'  @param W Mean worm burden or infection intensity
#'  @param Prev Prevalence in th epopulation
#'
#'  @return Estimate of the clumping parameter
#'  @export

  prev_intensity_to_clump <- function(W, Prev){
    rootSolve::uniroot.all(function(x) 1-Prev-(1+W/x)^-x,
                           interval = c(0,10))
  }

#' Clumping parameter as a function of worm burden
#'
#'   Relaxes the assumption of constant worm burden through time and
#'   instead models the clumping parameter as a function of worm burden
#'
#'   @param
