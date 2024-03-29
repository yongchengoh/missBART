% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypers.R
\name{hypers_list}
\alias{hypers_list}
\title{BART prior parameters}
\usage{
hypers_list(
  mu0 = 0,
  kappa = NULL,
  alpha = NULL,
  V = NULL,
  df = 10,
  q = 0.75,
  ...
)
}
\arguments{
\item{mu0}{hyperparameter}

\item{kappa}{hyperparameter}

\item{alpha}{hyperparameter}

\item{V}{hyperparameter}

\item{df}{hyperparameter}

\item{q}{hyperparameter}

\item{...}{Catches unused arguments}
}
\value{
list of hyperparameters
}
\description{
BART prior parameters
}
\examples{
hypers_list(mu0 = 5, kappa = 2, alpha = 1, V = 2)
}
