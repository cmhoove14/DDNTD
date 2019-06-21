#' Get estimate of effective reproduction number, $R_eff$
#'
#' Estimates the effective reproduction number, $R_eff$ as a function of
#' model parameters, population mean worm burden, and clumping parameter of negative binomial worm burden distribution
#'
#' @param parameters Model parameters
#' @param W mean worm burden in human population
#' @param k clumping parameter of the worm burden among human population
#'
#' @return estimate of the effective reproduction number, $R_eff$
#' @export

getReff_addNDD<-function(parameters, W, k){
  with(as.list(parameters),{

  k1<-sigma/(mu_N+mu_I)
  k2<-beta*0.5*W*H*phi_Wk(W,k)*f_Wgk(W, gam, k)/(mu_N+sigma)

  Num_1<-lambda*R_Wv(W, v)*k1*k2*C
  Num_2<- ( (1+k2)*f_N ) - ( mu_N + (beta*0.5*H*W*phi_Wk(W,k)*f_Wgk(W, gam, k) ) )
  Den<-( 1+k2+(k1*k2) ) * ( 1+k2 ) * ( mu_W+mu_H )* W* f_N

  Reff<-as.numeric(Num_1*Num_2/Den)

  return(Reff)
  })
}
