% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/schisto_helper_functions.R
\name{convert_burden_egg_to_worm}
\alias{convert_burden_egg_to_worm}
\title{Function to get solutions for worm burden and clumping parameter from input egg burden and prevalence}
\usage{
convert_burden_egg_to_worm(egg_burden, prevalence, m, gamma)
}
\arguments{
\item{egg_burden}{observed mean egg burden in population fraction}

\item{prevalence}{observed prevalence in population fraction}

\item{m}{peak egg output per 10mL per mated female worm}

\item{gamma}{density dependent fecundity parameter}
}
\value{
estimate of the mean worm burden and clumping parameter
}
\description{
Uses `egg_to_worm_fx` within `rootSolve::multiroot` equation solver to come up with estimates of mean worm burden
and clumping parameter from input egg burden, prevalence and egg output parameters
}
