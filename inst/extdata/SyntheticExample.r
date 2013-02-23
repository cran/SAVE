# Synthetic Example
library(SAVE)

#######
# load data
#######

data(synthfield,package='SAVE')
data(synthmodel,package='SAVE')

##############
# create the SAVE object which describes the problem and
# compute the corresponding mle estimates
##############

synth <- SAVE(response.name="y",controllable.names="x1",
calibration.names="x2",field.data=synthfield,model.data=synthmodel,
mean.formula=as.formula("~1+x1"),bestguess=list(x2=1.5))

##############
# Bayesian fit
##############

synth <- bayesfit(object=synth,
     prior=c(uniform("x2", upper=3., lower=0.)),
     n.iter=20000)
     
summary(synth)
plot(synth, option="trace")
plot(synth, option="calibration")
plot(synth, option="precision")

##############
# validate at a grid of x1 points
##############

xnew <- data.frame(seq(from=0.05,to=3.05,length=25))    
colnames(xnew)<-"x1"

valsynth <- validate(object=synth,newdesign=xnew,n.burnin=100)
summary(valsynth)
plot(valsynth)

#############
# using predictreality
#############

prsynth <- predictreality(object=synth,
     newdesign=xnew, tol=1.E-12, n.burnin=100)
     
summary(prsynth)
plot(prsynth,option="biascorr")
plot(prsynth,option="biasfun")

real <- prsynth@modelpred + prsynth@biaspred
av.real <- apply(real,2,mean)
delta.u <- apply(real,2,quantile,probs=0.975)
delta.l <- apply(real,2,quantile,probs=0.025)

par(mfrow=c(1,1))
a <- min(delta.l)
b <- max(delta.u)
plot(xnew$x1,av.real,ty="n",ylim=c(0.95*a,1.05*b),
     main="Bias-corrected Prediction",
     xlab="x1",ylab="y")
lines(xnew$x1,av.real,lty=1)
lines(xnew$x1,delta.l,lty=2)     
lines(xnew$x1,delta.u,lty=2)
points(field$x1,field$y,pch="*")

###########
# using predictcode to predict code at posterior mode
###########

ustar <- mean(synth@mcmcsample[-(1:100),"x2"])
xnewpure <- data.frame(xnew,"x2"=rep(ustar,length(xnew$x1)))

pcsynth <- predictcode(synth, newdesign=xnewpure, n.iter=20000, tol=1.E-12)

# compare simulation-based and exact estimates 
samples <- pcsynth@samples
av <- apply(samples,2,mean)
u <- apply(samples,2,quantile,probs=0.975)
l <- apply(samples,2,quantile,probs=0.025)

mm <- pcsynth@modelmean
v <- diag(pcsynth@covmat)
d <- 1.96*sqrt(v)

par(mfrow=c(1,1))
a <- min(u)
b <- max(l)
plot(xnew$x1,av,ty="n",ylim=c(0.95*a,1.05*b),
     main="Emulation",
     xlab="x1",ylab="y")
lines(xnew$x1,av,lty=1)
lines(xnew$x1,l,lty=2)     
lines(xnew$x1,u,lty=2)

points(xnew$x1,mm,pch="*",col=2)
points(xnew$x1,mm+d,pch="+",col=2)     
points(xnew$x1,mm-d,pch="+",col=2)

# tolerance bounds for this pure-model prediction
tmpdata <- matrix(mm,ncol=dim(real)[2],nrow=dim(real)[1],
                  byrow=T)
tmpdata <- real - tmpdata
tmpdata <- apply(tmpdata,2,abs)
tau.pure <- apply(tmpdata,2,quantile,0.90)

par(mfrow=c(1,1))
a <- min(mm-tau.pure)
b <- max(mm+tau.pure)
plot(xnew$x1,mm,ty="n",ylim=c(0.95*a,1.05*b),
     main="Pure-model prediction",
     xlab="x1",ylab="y")
lines(xnew$x1,mm,lty=1)
lines(xnew$x1,mm-tau.pure,lty=2)
lines(xnew$x1,mm+tau.pure,lty=2)
points(field$x1,field$y,pch="*")
