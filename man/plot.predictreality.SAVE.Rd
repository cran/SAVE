\name{plot.predictreality.SAVE}
\alias{plot,predictreality.SAVE-method}
\title{
A function for plotting summaries of an object of class \code{predictreality.SAVE}.
}
\description{
Two different plots to summarize graphically the results in an object of class \code{predicreality.SAVE}.
}
\usage{
\S4method{plot}{predictreality.SAVE}(x, option = "trace", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{predictreality.SAVE}
  	
}
  \item{option}{One of "biascorr" or "biasfun" (see details)
}
  \item{...}{Additional graphical parameters to be passed}
}
\details{
If \code{option}="biascorr" this function returns a plot with punctual predictions and 95\% tolerance bounds of reality at the given set of controllable inputs. If \code{option}="biasfun" the plot represents the estimated bias and 95\% tolerance bounds.
}

\author{ Jesus Palomo, Rui Paulo and Gonzalo Garcia-Donato. }

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

#Plot the results:
#Represent reality and tolerance bounds:
plot(prsw, option="biascorr")
#Represent bias and tolerance bounds:
plot(prsw, option="biasfun")


}
}

