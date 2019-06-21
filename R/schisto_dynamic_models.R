#' Basic schistosomiasis model
#'
#' A dynamic schistosomiasis model with negative (crowding induced reductions in fecundity)
#' and positive (mating limitation) density dependencies and functionality
#' to simulate mass drug administration, snail control, and other interventions.
#'
#' @param t Vector of timepoints to return state variable estiamtes
#' @param n Vector of state variable initial conditions
#' @param parameters Named vector of model parameters
#'
#' @return A matrix of the state variables at all requested time points
#' @export

schisto_base_mod = function(t, n, parameters) {
  with(as.list(parameters),{

    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I

    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations


  #Density dependencies based on worm burden
    f = f_Wgk(W, gam, k)
    R = R_Wv(W, v)
    phi = phi_Wk(W = W, k = k)  #Mating probability

  #Miracidial estimate assuming 1:1 sex ratio, mating probability, density dependence
    M=((0.5*W*H)*phi*f)#*m*u_H*(v*vq)


    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human
    dWtdt= (lambda*I*R) - ((mu_W+mu_H)*Wt)
    dWudt= (lambda*I*R) - ((mu_W+mu_H)*Wu)



    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  })
}

#' Simulate the base schistosomiasis model
#'
#' Uses the `schisto_base_mod` function to simulate the model through time
#' given parameter set, time frame, starting conditions, and potential events (e.g. MDA)
#'
#' @param nstart Named vector of starting values for state variables
#' @param time Numeric vector of times at which state variables should be estimated
#' @param model Name of the ode function to use, defaults to `schisto_base_mod`
#' @param parameters Named vector or list of parameter values
#'
#' @return dataframe of state variable values at requested times
#' @export
#'
sim_schisto_base_mod <- function(nstart,
                                 time,
                                 model = schisto_base_mod,
                                 parameters){
  as.data.frame(ode(nstart, time, model, parameters))
}
