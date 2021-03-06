% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{pval_vT_cov}
\alias{pval_vT_cov}
\title{Calculate selective p-value for given covariance}
\usage{
pval_vT_cov(vT, VCOV, this_y, nrSamples, checkFun, twosided = TRUE,
  bayesian = FALSE, alpha = 0.05, maxiter = 10, trace = TRUE, ...)
}
\arguments{
\item{vT}{test vector of function}

\item{VCOV}{covariance used for distribution of test statistic}

\item{this_y}{original response vector}

\item{nrSamples}{number of samples to be used}

\item{checkFun}{function; function to congruency with initial selection}

\item{twosided}{logical; compute two-sided p-values?}

\item{bayesian}{see \code{\link{mocasin}}}

\item{alpha}{see \code{\link{mocasin}}}

\item{maxiter}{maximum number of iteratoins to perform the linesearch used 
in the sampling procedure}

\item{trace}{see \code{\link{mocasin}}}

\item{...}{further arguments passed to vT if specified as function}
}
\description{
Calculate selective p-value for given covariance
}
