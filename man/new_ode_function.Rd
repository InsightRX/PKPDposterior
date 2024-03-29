% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/new_ode_function.R
\name{new_ode_function}
\alias{new_ode_function}
\title{Generates the ODE function that is used in a Stan/Torsten model}
\usage{
new_ode_function(ode, parameter_names, n_cmt)
}
\arguments{
\item{ode}{Stan code for ODE block}

\item{parameter_names}{vector of parameter names. Needs to be in order defined by the
chosen Torsten analytic solver.}

\item{n_cmt}{number of compartments in the ODE}
}
\value{
Character vector containing valid Stan code with ODE function
}
\description{
Generates the ODE function that is used in a Stan/Torsten model
}
