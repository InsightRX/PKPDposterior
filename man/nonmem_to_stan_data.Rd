% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonmem_to_stan_data.R
\name{nonmem_to_stan_data}
\alias{nonmem_to_stan_data}
\title{Reshape NONMEM data into PKPDposterior-ready or Stan-ready data}
\usage{
nonmem_to_stan_data(
  nm,
  covariate_cols,
  parameters,
  iiv,
  ruv,
  covariates_implementation = list(),
  dose_cmt = 1,
  ltbs = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{nm}{NONMEM data frame for a single patient, with columns for TIME,
EVID, DV, and AMT. Observations in DV are assumed to be taken from
compartment 2 unless otherwise specified in an (optional) \code{"CMT"}
column.}

\item{covariate_cols}{names of columns in \code{nm} that should be used as
covariates}

\item{parameters}{list of population parameters, e.g. \code{list(CL = 5, V = 50)}}

\item{iiv}{list of inter-individual variability for parameters. Should have
exact same list elements as \code{parameters}, and magnitude supplied on SD scale.}

\item{ruv}{magnitude of residual unexplained variability (RUV). Should be a
list specifying proportional and/or additive error magnitude on standard
deviation scale, e.g. \code{list("prop" = 0.1, "add" = 1)}. If \code{ltbs} is TRUE,
should specify only an \code{add} part, which applies an additive error on the
log-scale (which then becomes an approximate proportional error).}

\item{covariates_implementation}{named list indicating which implementation
method should be used for each covariate. See \link[PKPDsim:new_covariate]{PKPDsim::new_covariate} for
options.}

\item{dose_cmt}{Specify what dose compartment. Observation compartment in
dataset is irrelevant, handled in model.}

\item{ltbs}{use log-transform-both-sides approach for observations? Default
is \code{FALSE}.}

\item{verbose}{verbosity}
}
\value{
Named list suitable for passing on to Torsten.
}
\description{
Reshape NONMEM data into PKPDposterior-ready or Stan-ready data
}
