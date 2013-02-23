\name{predictreality}
\alias{predictreality}
\alias{predictreality.SAVE}
\alias{predictreality.SAVE-class}
\alias{predictreality,SAVE-method}
\alias{summary.predictreality.SAVE}
\alias{summary,predictreality.SAVE-method}
\alias{show.summary.predictreality.SAVE}
\alias{show,predictreality.SAVE-method}
\alias{show,summary.predictreality.SAVE-method}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predicting values of reality at new input points }
\description{
The emulator of the computer model and the Bayesian fit are used to produce samples from the posterior distribution of the computer model and from the bias function at associated new entry points. Then, bias-corrected predictions of the response under study can be obtained through addition (model+bias). We call these, predictions of reality.
}
\usage{
\S4method{predictreality}{SAVE}(object, newdesign, n.burnin=0, n.thin=1, tol=1E-10, ...)

\S4method{summary}{predictreality.SAVE}(object)

\S4method{show}{summary.predictreality.SAVE}(object)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of the corresponding signature. }
  \item{newdesign}{ A named matrix containing the points (controllable inputs) where predictions are to be performed. Column names should contain the \code{object@controllablenames}.}
  \item{n.burnin}{ The burnin to be applied (see details below).} 
  \item{n.thin}{ The thinin to be applied (see details below).}   
  \item{tol}{ The tolerance in the Cholesky decomposition.}
  \item{...}{ Extra arguments to be passed to the function (still not implemented).}
}

\value{
	Returns an S4 object of class \code{predictreality.SAVE} with the following slots:
  \item{modelpred}{A matrix with the simulations from the posterior distribution of the computer model output evaluated at the new design}
  \item{biaspred }{ A matrix with the simulations from the posterior distribution of the bias function evaluated at the new design. }
  \item{newdesign }{ A copy of the design given as argument. }
  \item{predictrealitycall }{ The call to the function. }

}

\details{
Posterior predictions of the computer model and bias over a given set of controllable inputs are simulated based on a Monte Carlo Markov Chain sampling mechanism. If \code{n.thin=1} and \code{n.burnin=0} the simulations have the same size as the sample size in \code{object@mcmcsample}. Otherwise, the first \code{n.burnin} are discharged and only each of \code{n.burnin} are kept. 

A first analysis of the resulting sample can be performed with \code{summary} which provides punctual estimations and tolerance bounds of the resulting predictions.
}

\author{ Jesus Palomo, Rui Paulo and Gonzalo Garcia-Donato. }

\seealso{  \code{"\link[=plot,predictreality.SAVE-method]{predictreality}"}, 
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

# summary of the results

summary(gfsw)

##############
# obtain the posterior distribution of the unknown parameters
##############

gfsw <- bayesfit(object=gfsw, prior=c(uniform("t", upper=8, lower=0.8)), n.iter=20000, n.burnin=100, n.thin=2)

#########
# bias-corrected prediction at a set of inputs
# using predictreality
##########

load <- c(4.0,5.3)
curr <- seq(from=20,to=30,length=20)
g <- c(1,2)

xnew <- as.data.frame(expand.grid(curr,load,g))
names(xnew)<-c("C","L","G")

# Obtain samples
prsw <- predictreality(object=gfsw, newdesign=xnew, tol=1.E-12)

#Summarize the results:
summary(prsw)

}	
	
	}