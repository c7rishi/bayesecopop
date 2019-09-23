% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_cylmix.R
\name{fit_cyl_mix}
\alias{fit_cyl_mix}
\title{Bayesian mixture modeling for cylendrical data}
\usage{
fit_cyl_mix(data, ncomp = 2, overfit = "super-impose", seed = 1,
  sampling_method = "mcmc", init_opt = NULL, angle_col = "angle",
  line_col = "line", rand_init = (sampling_method == "mcmc"),
  show_allocation = FALSE, ...)
}
\arguments{
\item{data}{a data frame, with one column providing data on the
angular variable, and one column on the linear variable.}
}
\description{
Bayesian mixture modeling for cylendrical data
}