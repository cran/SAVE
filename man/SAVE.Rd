\name{SAVE}
\alias{SAVE}
\alias{show.SAVE}
\alias{show,SAVE-method}
\alias{show,summary.SAVE-method}
\alias{summary.SAVE}
\alias{summary.SAVE-class}
\alias{summary,SAVE-method}
\alias{show.summary.SAVE}
\title{ Functions to start the Statistical Analysis and Validation Engine of a computer model and corresponding exploration of results }
\description{\code{SAVE} starts the main object for subsequent analyses.}
\usage{
SAVE(response.name, controllable.names, calibration.names, field.data, model.data,
             mean.formula, bestguess)

\S4method{show}{SAVE}(object)

\S4method{summary}{SAVE}(object)

\S4method{show}{summary.SAVE}(object)
}
\arguments{
  \item{response.name}{A \code{character} object with the name of the response variable.}
  \item{controllable.names}{A \code{character} object with the names of the controllable inputs.}
  \item{calibration.names}{A \code{character} object with the names of the calibration inputs. Use NULL if the model does not depend on any calibration inputs.}
  \item{field.data}{A \code{data.frame} with the observations coming from the field experiment. The column names of this \code{data.frame} should contain the \code{controllable.names} and the \code{response.name}.}
  \item{model.data}{A \code{data.frame} with the observations coming from runs of the computer model. The column names of this \code{data.frame} should contain the \code{controllable.names}, the \code{calibration.names} and the \code{response.name}.}
  \item{mean.formula}{A \code{formula} specifying the mean function of the Gaussian Process approximation to the output of the computer model (the emulator). Can only involve terms in \code{controllable.names}.}
  \item{bestguess}{A named \code{list} specifying the best guess for the calibration parameters.} 
  \item{object}{An object of the corresponding signature.}
}  
     
\details{This function fits the emulator of the computer model and an estimation of the bias function. Calculations are partially based on the package \code{DiceKriging}}

\value{
  \code{SAVE} returns an S4 object of class \code{SAVE} (see \code{SAVE-class}). 
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
  
\seealso{ \code{"\link[=plot,predictcode.SAVE-method]{predictcode}"},
	  	  \code{\link[=bayesfit.SAVE]{bayesfit}}  
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
}

}


