area = 200
H = 300

base_pars <- c("f_N" = 0.1,
               "C" = 50*area,
               "mu_N" = 1/60,
               "sigma" = 1/40,
               "mu_I" = 1/10 - 1/60,
               "mu_W" = 1/(365*3.3),
               "H" = H,
               "mu_H" = 1/(60*365),
               "lambda" = 2.074609e-04,
               "beta" = 1.6e-6,
               "cvrg" = 0.43,
               "zeta" = 5e-3,
               "xi" = 2.8e-3)

usethis::use_data(base_pars, overwrite = TRUE)

age_strat_pars <- c(
  ##standard snail parameters
    r=0.10,             # recruitment rate (from sokolow et al)
    mu_N=1/60,          # Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=1/30,         # Transition rate from exposed to infected (assuming pre-patency period of ~4 weeks) doi:10.4269/ajtmh.16-0614
    mu_I=1/10,          # Increased mortality rate of infected snails
    theta=500,          # mean cercarial shedding rate per adult snail doi:10.4269/ajtmh.16-0614

  #Adult Worm, Miracidia and Cercariae Parameters
    mu_W = 1/(4*365),   # death rate of adult worms
    mu_H_A = 1/(50*365),# death rate of adult humans
    mu_H_C = 1/(70*365),# death rate of children
    m = 10,             # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = 0.08,           # mean egg viability (miracidia per egg)

  #Density dependence parameters
    zeta = 5e-3,        # parameter of fecundity reduction function
    xi = 2.8e-3,        # parameter for acquired immunity function http://doi.wiley.com/10.1111/j.1365-3024.1992.tb00029.x

  #Human parameters
    U_C = 110,          # mL urine produced per child per day /10mL https://doi.org/10.1186/s13071-016-1681-4
    U_A = 130           # mL urine produced per adult per day /10mL https://doi.org/10.1186/s13071-016-1681-4
)

usethis::use_data(age_strat_pars, overwrite = TRUE)


garch_pars <- c(
  beta_e = 4.06e-5, #environmental reservoir to human transmission rate
  beta_d = 0, #human to human transmission rate
  gamma = 1/(365*3), #rate of recover from infected to susceptible
  omega = 0, #recruitment of infectious agents not dependent on human infection (exogenous propogules)
  V = 1, #abundance of vectors/intermediate hosts/suitable environment for transmission
  lambda = 1, #recruitment rate of infectious agents by infectious individuals
  sigma = 1, #fraction of infectious agents produced that reach the environment
  rho = 1/90 #mortality rate of infectious agents in the environment
)

usethis::use_data(garch_pars, overwrite = TRUE)
