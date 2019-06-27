base_pars <- c("f_N" = 0.1,
               "C" = 1e4,
               "mu_N" = 1/60,
               "sigma" = 1/40,
               "mu_I" = 1/10 - 1/60,
               "mu_W" = 1/(365*3.3),
               "H" = 300,
               "mu_H" = 1/(60*365),
               "k" = 0.08,
               "lambda" = 2.074609e-04,
               "beta" = 1.6e-6,
               "cvrg" = 0.43,
               "gamma" = 1e-3,
               "v" = 2.8e-3)

usethis::use_data(base_pars, overwrite = TRUE)
