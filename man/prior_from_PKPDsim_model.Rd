% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_from_PKPDsim_model.R
\name{prior_from_PKPDsim_model}
\alias{prior_from_PKPDsim_model}
\title{Get prior for Torsten from a PKPDsim model package}
\usage{
prior_from_PKPDsim_model(model, map = list(), drop = c())
}
\arguments{
\item{model}{package name}

\item{map}{remap model parameters from Stan to PKPDsim syntax, specified as a
list. E.g. \code{map = list(V1 = V)}}

\item{drop}{vector of parameter names that we are not interested in or don't
need to provide.}
}
\description{
Get prior for Torsten from a PKPDsim model package
}
