\name{bayesfit.SAVE}
\alias{bayesfit.SAVE}
\alias{bayesfit,SAVE-method}
\alias{bayesfit}

\title{Bayesian methods for the analysis of a computer model}
\description{ Fitting into a Bayesian framework an object of the class \code{SAVE}. }
\usage{
\S4method{bayesfit}{SAVE}(object, prior, mcmcMultmle=1, prob.prop=0.5, method=2, n.iter, nMH=20, n.burnin=0, n.thin=1, ...)
}
\arguments{
  \item{object}{A \code{SAVE} object normally produced by a call to the \code{\link{SAVE}} function.}
  \item{prior}{The prior distribution for the calibration parameters. This should be specified concatenating, \code{c}, calls to functions \code{\link{uniform}} and/or \code{\link{normal}}. See details below.}
  \item{mcmcMultmle}{A factor that is used to specify the prior for lambdaF and lambdaM (the precision parameters in the field and the bias, respectively).}
  \item{prob.prop}{The probability of proposing from the prior in the Metropolis-Hastings algorithm. See details below.}
  \item{method}{Method implemented in the Gibbs sampling. See details below.}
  \item{n.iter}{Number of total iterations. See details below.}
  \item{nMH}{Number of Metropolis-Hastings steps. See details below.} 
  \item{n.burnin}{How many of the first simlations are discharged. See details below.} 
  \item{n.thin}{The thinin to be applied to the resulting MCMC sample. See details below.}
  \item{...}{Extra arguments to be passed to the function (still not implemented).}
}  
     
\details{This function runs an MCMC to obtain samples from the posterior of (lambdaM, lambdaF, calibration parameters). 
	
About the prior: the prior for lambdaM and lambdaF is the product of two independent exponential densities with mean the corresponding maximum likelihood estimators times the factor defined in \code{mcmcMultmle}. The prior for each calibration parameter can be either a uniform distribution or a truncated normal and should be specified concatenating calls to \code{\link{uniform}} and/or \code{\link{normal}}. For example, for a problem with calibration parameters named ("delta1", "shift") a possible prior could be\cr

 \code{prior=c(uniform(var.name="delta1", lower=0, upper=1), normal(var.name="shift", mean=2, sd=1, lower=0, upper=8))}

About the MCMC: the algorithm implemented is based on a Gibbs sampling scheme. If \code{method=2} then computer model and bias are integrated out (analytically) and only the full conditionals for lambdaM, lambdaF and calibration parameters are sampled from. When \code{method=1} the computer model and the bias form part of the sampling scheme. In any case, the corresponding full conditional for the calibration parameters is sampled via a Metropolis-Hastings algorithm with proposal being a mixture of the prior specified and a local move. In this mixture, the probability of a proposal coming from the prior is specified by \code{prob.prop}. This Metropolis-Hastings algorithm is run \code{nMH} times before the sample is accepted. The MCMC is run with \code{n.iter} iterations to which the first \code{n.burnin} are discharged and to the rest a thinin of \code{n.thin} is applied.  
 }

\value{
  \code{SAVE} returns a copy of the \code{SAVE} object used as argument to the function, but with the following slots filled (or replaced if they where no empty)
 \describe{ 
    \item{\code{method}:}{the value given to \code{method}.}
    \item{\code{mcmcMultmle}:}{the value given to \code{mcmcMultmle}.}
    \item{\code{n.iter}:}{the value given to \code{n.iter}.}
    \item{\code{nMH}:}{The value given to \code{nMH}.}
    \item{\code{mcmcsample}:}{A \code{matrix} with the  result of the MCMC sampling after the burnin and thinin has been applied.}
    \item{\code{bayesfitcall}:}{The \code{call} to \code{bayesfit}.}
  }
  
  }

\references{
Bayarri MJ, Berger JO, Paulo R, Sacks J, Cafeo JA, Cavendish J, Lin CH, Tu J (2007). A Framework for Validation of Computer Models. Technometrics, 49, 138-154.

Craig P, Goldstein M, Seheult A, Smith J (1996). Bayes linear strategies for history matching of hydrocarbon reservoirs. In JM Bernardo, JO Berger, AP Dawid, D Heckerman, AFM Smith (eds.), Bayesian Statistics 5. Oxford University Press: London. (with discussion).

Higdon D, Kennedy MC, Cavendish J, Cafeo J, Ryne RD (2004). Combining field data and computer simulations for calibration and prediction. SIAM Journal on Scientific Computing, 26, 448-466.

Kennedy MC, O Hagan A (2001). Bayesian calibration of computer models (with discussion). Journal of the Royal Statistical Society B, 63, 425-464.

Roustant O., Ginsbourger D. and Deville Y. (2012).
  DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer
  Experiments by Kriging-Based Metamodeling and Optimization. Journal
  of Statistical Software, 51(1), 1-55.
  
}


\author{Jesus Palomo, Rui Paulo and Gonzalo Garcia-Donato}
  
\seealso{ \code{"\link[=plot,SAVE-method]{plot}"},
	  	  \code{\link[=predictreality.SAVE]{predictreality}},
	      \code{\link[=validate,SAVE-method]{validate}}
}
	
\examples{
\dontrun{
library(SAVE)

#############
# load data
#############

data(spotweldfield,package='SAVE')
data(spotweldmodel,package='SAVE')

##############
# create the SAVE object which describes the problem and
# compute the corresponding mle estimates
##############

gfsw <- SAVE(response.name="N", controllable.names=c("C", "L", "G"), calibration.names=c("t"), field.data=spotweldfield, model.data=spotweldmodel, mean.formula=as.formula("~1"), bestguess=list(t=4.0))

##############
# obtain the posterior distribution of the unknown parameters 
##############

gfsw <- bayesfit(object=gfsw, prior=c(uniform("t", upper=8, lower=0.8)), n.iter=20000, n.burnin=100, n.thin=2)

# summary of the results
summary(gfsw)
	
	}
}

