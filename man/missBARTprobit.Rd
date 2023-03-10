% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missBARTprobit.R
\name{missBARTprobit}
\alias{missBARTprobit}
\title{BART and probit regression for data with non-ignorable missingness in the response.}
\usage{
missBARTprobit(
  x,
  y,
  x_predict = NA,
  n_trees = 20,
  burn = 100,
  iters = 100,
  thin = 2,
  predict = TRUE,
  tree_prior_params = tree_list(),
  hypers = hypers_list(),
  scale = TRUE,
  include_x = TRUE,
  include_y = TRUE,
  show_progress = TRUE,
  progress_every = 10,
  ...
)
}
\arguments{
\item{x}{BART covariates}

\item{y}{data}

\item{x_predict}{out-of-sample covariates. If not specificied, the default is set to NA and no out-of-sample predictions will be made.}

\item{n_trees}{number of BART trees. Default is set to 90.}

\item{burn}{burn-in samples. Default is set to 1000.}

\item{iters}{post-burn-in samples (after thinning). Default is set to 1000.}

\item{thin}{thinning. Default is set to 3.}

\item{predict}{whether or not to make out-of-sample predictions}

\item{tree_prior_params}{prior parameters for BART trees}

\item{hypers}{prior parameters for BART parameters}

\item{scale}{logical. Whether to scale data to range (-0.5, 0.5).}

\item{include_x}{logical. Include x in probit model?}

\item{include_y}{logical. Include y in probit model?}

\item{show_progress}{logical.}

\item{progress_every}{integer value stating how often to update the progress bar.}

\item{...}{Catches unused arguments}
}
\value{
a list containing BART predictions and imputed values
}
\description{
BART and probit regression for data with non-ignorable missingness in the response.
}
\examples{
x <- matrix(runif(6), ncol = 2)
y <- matrix(runif(6), ncol = 2) \%*\% matrix(rnorm(4), ncol=2)
missBARTprobit(x, y, n_trees = 2, burn = 2, iters = 2, thin = 1, scale = FALSE)
}
