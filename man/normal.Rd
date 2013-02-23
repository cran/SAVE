\name{normal,uniform}
\alias{normal}
\alias{uniform}
\title{ Auxiliary functions. }
\description{ Auxiliary functions. }
\usage{
uniform(var.name, lower, upper)

normal(var.name, mean, sd, lower, upper)
}
\arguments{
  \item{var.name}{The name of the variable whose prior distribution is being specified.}
  \item{lower}{The lower prior bound for \code{var.name}.}
  \item{upper}{The upper prior bound for \code{var.name}.}
  \item{mean}{The prior mean for \code{var.name}.}
  \item{sd}{The prior standard deviation for \code{var.name}.}
}  
     
\details{ This functions are intended as user-friendly way of specifying the type of priors assumed for the calibration parameters with the function \code{bayesfit}.

 For example, for a problem with calibration parameters named ("delta1", "shift") a possible prior could be\cr

 \code{prior=c(uniform(var.name="delta1", lower=0, upper=1), normal(var.name="shift", mean=2, sd=1, lower=0, upper=8))}
 }

