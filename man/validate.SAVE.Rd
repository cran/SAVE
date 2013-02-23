\name{validate.SAVE}
%\alias{validate.SAVE}
\alias{validate.SAVE-class}
\alias{validate,SAVE-method}
\alias{validate}
\alias{summary.validate.SAVE}
\alias{summary,validate.SAVE-method}
\alias{show.summary.validate.SAVE}
%\alias{plot.validate.SAVE}
\alias{plot,validate.SAVE-method}
\alias{show,summary.validate.SAVE-method}
\alias{show,validate.SAVE-method}


\title{Validation of a computer model}
\description{Checking the validity of a computer model at a given set of controllable inputs}
\usage{
\S4method{validate}{SAVE}(object, newdesign, calibration.value="mean", prob=0.90, n.burnin=0, n.thin=1, tol=1E-10, ...)

\S4method{summary}{validate.SAVE}(object)

\S4method{show}{summary.validate.SAVE}(object)

\S4method{plot}{validate.SAVE}(x, ...)

}
\arguments{
  \item{object}{ An object of corresponding signature. }
  \item{x}{ An object of class \code{validate.SAVE}. }
  \item{newdesign}{ A named matrix containing the points (controllable inputs) where predictions are to be performed. Column names should contain the \code{object@controllablenames}. }
  \item{calibration.value}{ The summary of the posterior distribution of the calibration parameters which is taken as the true value for the computer model. Possible values are "mean" and "median". }
  \item{prob}{ The probability used for the representation of the credible intervals. }
  \item{n.burnin}{ The burnin to be applied (see details below).} 
  \item{n.thin}{ The thinin to be applied (see details below).}   
  \item{tol}{ The tolerance in the Cholesky decomposition.}
  \item{...}{Additional arguments to be passed, still not implemented. }
}  
     
\value{
Returns an S4 object of class \code{validate.SAVE} with the following slots:
 \describe{ 
 	\item{\code{bayesfitcall}:}{ The \code{call} to \code{bayesfit}. }
	\item{\code{call}:}{ The original \code{call} to \code{SAVE} function. }
	\item{\code{newdesign}:}{ A copy of the design given. }
	\item{\code{validate}:}{ A matrix with the results of the validation analysis (see details below). }
	\item{\code{validatecall}:}{ The call to the \code{validate} function. }
	}
  
  }

\details{ Following the framework for the analysis of computer models by Bayarri et al (2007), validation of a computer model translates to the question: is the computer model producing results that are accurate enough for its intended use?

Answering this question implies the comparison of the responses of the computer model at the `true' value of the calibration parameters (\code{calibration.value}) with reality. This comparison should be performed at the set of controllable inputs that are of interest for the researcher (which in turn represent the `intended use' for the computer model) and that are passed as the argument \code{newdesign} to the function.  

For this comparison to be performed, \code{validation} returns a matrix (\code{@validate}) containing (for each input point on rows) the prediction of reality (column called bias.corrected) jointly with the estimation of the model response, both with corresponding accompanying confidence measures (columns called tau.bc and tau.pm respectively). These measures should be interpreted as:

Prob(|estimate-real value|<tau)=\code{prob}

Alternatively, the discrepancy between computer model and reality can be done through the estimated bias (column called bias) with the associated 100\code{prob}\% credible interval (columns called bias.Lower and bias.Upper).

In the calculations, the simulated sample from the posterior distribution contained in \code{object@mcmsample} is used. It is possible to discharge some of the samples in this matrix and not consider the first \code{n.burnin} or only keeping one each \code{n.thin}.

The results can be conveniently visualized with the functions \code{summary} and \code{plot}.

 }



\references{
Bayarri MJ, Berger JO, Paulo R, Sacks J, Cafeo JA, Cavendish J, Lin CH, Tu J (2007). A Framework for Validation of Computer Models. Technometrics, 49, 138-154.
}

\author{Jesus Palomo, Rui Paulo and Gonzalo Garcia-Donato}
  
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

##############
# validate the computer model at chosen set of controllable
# inputs
###############

load <- c(4.0,5.3)
curr <- seq(from=20,to=30,length=20)
g <- c(1,2)

xnew <- as.data.frame(expand.grid(curr,load,g))
names(xnew)<-c("C","L","G")

valsw <- validate(object=gfsw,newdesign=xnew,n.burnin=100)

# summary of results
summary(valsw)
# plot results
plot(valsw)


	}
}
