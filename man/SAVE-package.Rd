% File man/SAVE-package.Rd
\name{SAVE-package}
\alias{SAVE-package}
\docType{package}
\title{
Simulator Analysis and Validation Engine}
\description{
This package provides a statistical framework for the analysis of complex computer models, with special emphasis on aspects related with emulation, calibration and validation. The methodology implemented in SAVE is directly based on Bayarri et al (2007) but has roots also in the papers by Craig et al (1996), Kennedy and O'Hagan (2001) and Higdon et al (2004).
}
\details{
\tabular{ll}{
Package: \tab SAVE\cr
Type: \tab Package\cr
Version: \tab 0.9.3.3\cr
Date: \tab 2013-04-05\cr
License: \tab GPL-2\cr
}

Index:
\preformatted{
}
}
\author{
Jesus Palomo, Rui Paulo and Gonzalo Garcia-Donato

Maintainer: <jesus.palomo@urjc.es>
}
\references{
Bayarri MJ, Berger JO, Paulo R, Sacks J, Cafeo JA, Cavendish J, Lin CH, Tu J (2007). A Framework for Validation of Computer Models. Technometrics, 49, 138-154.

Craig P, Goldstein M, Seheult A, Smith J (1996). Bayes linear strategies for history matching of hydrocarbon reservoirs. In JM Bernardo, JO Berger, AP Dawid, D Heckerman, AFM Smith (eds.), Bayesian Statistics 5. Oxford University Press: London. (with discussion).

Higdon D, Kennedy MC, Cavendish J, Cafeo J, Ryne RD (2004). Combining field data and computer simulations for calibration and prediction. SIAM Journal on Scientific Computing, 26, 448-466.

Kennedy MC, O Hagan A (2001). Bayesian calibration of computer models (with discussion). Journal of the Royal Statistical Society B, 63, 425-464.

}
\keyword{package}
\keyword{computer models}
\keyword{simulation}
\keyword{calibration}
\keyword{validation}

\seealso{
\code{\link[SAVE]{SAVE}},
\code{\link[SAVE]{SAVE-class}}
}
\examples{
\dontrun{
	
#Validate the computer model in the Spotweld example (see Bayarri et al 2007 for details)
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

# summary of validation exercise:
summary(valsw)
# plot results
plot(valsw)
}
}
