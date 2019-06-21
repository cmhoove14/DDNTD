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

base_mod = function(t, n, parameters) {
  with(as.list(parameters),{

    S=n[1]
    E=n[2]
    I=n[3]
    Wt=n[4]
    Wu=n[5]
    N=S+E+I

    W=(cov*Wt) + ((1-cov)*Wu) #weighting treated and untreated populations
    #Miracidia production; function of adult female worms alive in the system (W*H*0.5) assuming 1:1 sex ratio,
  #DD functions

    f = f_Wgk(W, gam, k)
    R = R_Wv(W, v)
    phi = phi_Wk(W = W, k = k)  #Mating probability

    M=((0.5*W*H)*phi*f)#*m*u_H*(v*vq)


    dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*M*S #Susceptible snails

    dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails

    dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails

    #worm burden in human
    dWtdt= (lamda*I*R) - ((mu_W+mu_H)*Wt)
    dWudt= (lamda*I*R) - ((mu_W+mu_H)*Wu)



    return(list(c(dSdt,dEdt,dIdt,dWtdt,dWudt)))
  })
}
