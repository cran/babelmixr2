% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pknca.R
\name{nlmixr2Est.pknca}
\alias{nlmixr2Est.pknca}
\title{Estimate starting parameters using PKNCA}
\usage{
\method{nlmixr2Est}{pknca}(env, ...)
}
\arguments{
\item{env}{Environment for the nlmixr2 estimation routines.

This needs to have:

- rxode2 ui object in `$ui`

- data to fit in the estimation routine in `$data`

- control for the estimation routine's control options in `$ui`}

\item{...}{Other arguments provided to `nlmixr2Est()` provided for
flexibility but not currently used inside nlmixr}
}
\value{
A model with updated starting parameters.  In the model a new element
named "nca" will be available which includes the PKNCA results used for the
calculation.
}
\description{
Estimate starting parameters using PKNCA
}
\details{
Parameters are estimated as follows:
\itemize{
\item \code{ka} 4 half-lives to Tmax but not higher than 3:  \code{log(2)/(tmax/4)}
\item \code{vc} Inverse of dose-normalized Cmax
\item \code{cl} Estimated as the median clearance
\item \code{vp,vp2}2- and 4-fold the \code{vc}, respectively by default,
controlled by the \code{vpMult} and \code{vp2Mult} arguments to
\code{pkncaControl}
\item \code{q,q2} 0.5- and 0.25-fold the \code{cl}, respectively by default,
controlled by the \code{qMult} and \code{q2Mult} arguments to
\code{pkncaControl}
}

The bounds for the parameter estimates are set to 10\% of the first percentile
and 10 times the 99th percentile.  (For ka, the lower bound is set to the
lower of 10\% of the first percentile or 0.03 and the upper bound is not
modified from 10 times the 99th percentile.)

Parameter estimation methods may be changed in a future version.
}
