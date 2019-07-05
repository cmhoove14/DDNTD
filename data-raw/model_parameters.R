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
               "gamma" = 5e-3,
               "xi" = 2.8e-3)

usethis::use_data(base_pars, overwrite = TRUE)

age_strat_pars <- c(
  ##standard snail parameters
    f_N=0.10, # recruitment rate (from sokolow et al)
    C=50*area, # carrying capacity corresponding to 50 snails per square meter
    mu_N=1/60, #Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=1/40, #Transition rate from exposed to infected (assuming pre-patency period of 40 days)
    mu_I=1/10 - 1/60, # Increased mortality rate of infected snails

  #Adult Worm, Miracidia and Circariae Parameters
    mu_W = 1/(3.3*365), # death rate of adult worms
    m = 5.2,            # mean eggs shed per female worm per 10mL urine (truscott et al)
    v = 0.08,           # mean egg viability of eggs shed into environment

  #Density dependence parameters
    gamma = 5e-3,       # parameter of fecundity reduction function
    xi = 2.8e-3,        # parameter for acquired immunity funtion

  #Human parameters
    H = H,          # Total number of people
    prop_SAC = 0.5,       #Percent of people that are school age children (SAC)
    prop_adult = 0.5, #percent of people that are not school age children (assumed here to be adults)
    cvrg_SAC = 0.9,  #MDA coverage in SAC population
    cvrg_adult = 0,  #MDA coverage in adult population
    u_SAC = 50, # mL urine per SAC/day/10mL assumed to be half of adult (approximate, ranges from 20 - 100)
    u_adult = 100, # mL urine per adult/day/10mL (approximate, ranges from 80 - 200)
    rho_SAC = 0.3, # relative number of eggs shed by SAC that make it to snail habitat (~sanitation)
    rho_adult = 0.015, # relative number of eggs shed by adults that make it to snail habitat (~sanitation); 5% of SAC (truscott et al)
    omega_SAC = 1,     # relative infection risk of SAC (related to clean water access/education/water contact)
    omega_adult = 0.1,   # relative infection risk of adults (related to clean water access/education/water contact)

  #Transmission parameters
    lambda=1.2e-4, # snail-to-man transmission
    beta=1.6e-6  # man-to-snail transmission
)

usethis::use_data(age_strat_pars, overwrite = TRUE)
