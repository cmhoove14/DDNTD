require(readr)

schisto_dat <- readr::read_csv("data-raw/schistosoma_spp_infection_data.csv")

usethis::use_data(schisto_dat)
