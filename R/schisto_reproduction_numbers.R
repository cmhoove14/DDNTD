#' Get estimate of effective reproduction number, $R_eff$
#'
#' Estimates the effective reproduction number, $R_eff$ as a function of
#' model parameters, population mean worm burden, and clumping parameter of negative binomial worm burden distribution
#'
#' @param parameters Model parameters
#' @param W mean worm burden in human population
#' @param kap clumping parameter of the worm burden among human population
#'
#' @return estimate of the effective reproduction number, $R_eff$
#' @export

getReff<-function(parameters, W, kap){
  with(as.list(parameters),{

  Num_1<-lambda*R_Wv(W, v)*(sigma/(mu_N+mu_I))*(beta*0.5*W*H*phi_Wk(W,kap)*f_Wgk(W, gamma, kap)/(mu_N+sigma))*C

  Num_2<- ( (1+(beta*0.5*W*H*phi_Wk(W,kap)*f_Wgk(W, gamma, kap)/(mu_N+sigma)))*f_N ) - ( mu_N + (beta*0.5*H*W*phi_Wk(W,kap)*f_Wgk(W, gamma, kap) ) )

  Den<-( 1+(beta*0.5*W*H*phi_Wk(W,kap)*f_Wgk(W, gamma, kap)/(mu_N+sigma))+((sigma/(mu_N+mu_I))*(beta*0.5*W*H*phi_Wk(W,kap)*f_Wgk(W, gamma, kap)/(mu_N+sigma))) ) * ( 1+(beta*0.5*W*H*phi_Wk(W,kap)*f_Wgk(W, gamma, kap)/(mu_N+sigma)) ) * ( mu_W+mu_H )* W* f_N

  Reff<-as.numeric(Num_1*Num_2/Den)

  return(Reff)
  })
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

  Num_1<-lambda*R_Wv(W, v)*(sigma/(mu_N+mu_I))*(beta*0.5*W*H*f_Wgk(W, gamma, kap)/(mu_N+sigma))*C

  Num_2<- ( (1+(beta*0.5*W*H*f_Wgk(W, gamma, kap)/(mu_N+sigma)))*f_N ) - ( mu_N + (beta*0.5*H*W*f_Wgk(W, gamma, kap) ) )

  Den<-( 1+(beta*0.5*W*H*f_Wgk(W, gamma, kap)/(mu_N+sigma))+((sigma/(mu_N+mu_I))*(beta*0.5*W*H*f_Wgk(W, gamma, kap)/(mu_N+sigma))) ) * ( 1+(beta*0.5*W*H*f_Wgk(W, gamma, kap)/(mu_N+sigma)) ) * ( mu_W+mu_H )* W* f_N

  Reff<-as.numeric(Num_1*Num_2/Den)

  return(Reff)
  })
}
