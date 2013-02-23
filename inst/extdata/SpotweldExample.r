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

# summary of the results
summary(gfsw)
# traceplots
plot(gfsw, option="trace", col=2)
# posterior of the calibration parameter
plot(gfsw, option="calibration", col=4, lty=3)
# posterior of the measurement error precision and
# bias function precision
plot(gfsw, option="precision", col=2, lty=4)

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


##########
# emulate the output of the model using predictcode
##########

# construct design at which to emulate the model
u <- 3.2
load <- c(4.0,5.3)
curr <- seq(from=20,to=30,length=20)
g <- c(1,2)

xnewpure <- expand.grid(curr,load,g)
xnewpure <- cbind(xnewpure,rep(u,dim(xnewpure)[1]))
names(xnewpure) <- c("C","L","G","t")
xnewpure <- as.data.frame(xnewpure)

pcsw<- predictcode(object=gfsw, newdesign=xnewpure, n.iter=20000, tol=1.E-12)

# Plot results
paths <- pcsw@samples
# compare estimates using simulation and the explicit formulae
avpure <- apply(paths,2,mean)
meanpure <- pcsw@modelmean

qts <- apply(paths,2,quantile,probs=c(0.025,0.975))
qtspure <- 1.96*sqrt(diag(pcsw@covmat))

par(oma=c(1,1,2,0),mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
    
    # explicit formulae
    plot(curr,meanpure[v],type="l",ylim=c(3,9),
    xlab="current",ylab="weld diameter",col=2)
    lines(curr,meanpure[v]+qtspure[v],lty=1,col=2)
    lines(curr,meanpure[v]-qtspure[v],lty=1,col=2)
    text(22,9,paste("gauge= ",g[i],", 
    load=",load[j],sep=""),cex=0.8,pos=1)
    # simulation-based
    lines(curr,avpure[v],type="l",col=1,lty=3)
    lines(curr,qts[1,v],lty=3,col=1)
    lines(curr,qts[2,v],lty=3,col=1)
    
  }
}
mtext("Emulate output of code at u = 3.20 - Spotweld Data",side=3,outer=T,cex=1.2)

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

# Plot results
# reality = model + bias
real <- prsw@modelpred+prsw@biaspred

# estimate
av <- apply(real,2,mean)

# tolerance bounds
tmpdata <- matrix(av,ncol=dim(real)[2],nrow=dim(real)[1],
                  byrow=T)
tmpdata <- real - tmpdata
tmpdata <- apply(tmpdata,2,abs)
tau.real <- apply(tmpdata,2,quantile,0.90)

# plot

par(oma=c(1,1,2,0),mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
    plot(curr,av[v],type="l",ylim=c(3,9),
    xlab="current",ylab="weld diameter")
    lines(curr,av[v]+tau.real[v],lty=3)
    lines(curr,av[v]-tau.real[v],lty=3)
    text(22,9,paste("gauge= ",g[i],", 
    load=",load[j],sep=""),cex=0.8,pos=1)
    # field data that correspond to this situation
    v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
    data <- spotweldfield$N[v]
    inputs <- spotweldfield$C[v]
    points(inputs,data)
    #dev.off()
  }
}
mtext("Bias-corrected predictions - Spotweld Data",side=3,
      outer=T,cex=1.2)

# tolerance bounds for the pure model predictions
# at u = 3.2

# Note: avpure was computed above

tmpdata <- matrix(avpure,ncol=dim(real)[2],nrow=dim(real)[1], byrow=T)
tmpdata <- real - tmpdata
tmpdata <- apply(tmpdata,2,abs)
tau.pure <- apply(tmpdata,2,quantile,0.90)

# plot

par(oma=c(1,1,2,0),mfrow=c(2,2))
for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avpure[v],type="l",ylim=c(3,9),
			 xlab="current",ylab="weld diameter")
        lines(curr,avpure[v]+tau.pure[v],lty=3)
        lines(curr,avpure[v]-tau.pure[v],lty=3)
		text(22,9,paste("gauge= ",g[i],", 
						load=",load[j],sep=""),cex=0.8,pos=1)
# field data that correspond to this situation
		v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
		data <- spotweldfield$N[v]
		inputs <- spotweldfield$C[v]
		points(inputs,data)
#dev.off()
	}
}
mtext("Pure-model predictions - Spotweld Data",side=3,
outer=T,cex=1.2)

# plots using the output of validate()
# Note: valsw was computed above

avpure <- valsw@validate[,"pure.model"]
tau.pure <- valsw@validate[,"tau.pm"]

avbc <- valsw@validate[,"bias.corrected"]
taubc <- valsw@validate[,"tau.bc"]

avbias <- valsw@validate[,"bias"]
biasL <- valsw@validate[,"bias.Lower"]
biasU <- valsw@validate[,"bias.Upper"]

par(mfrow=c(3,4))
for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avpure[v],type="l",ylim=c(3.5,8.5),
			 xlab="current",ylab="weld diameter")
        lines(curr,avpure[v]+tau.pure[v],lty=3)
        lines(curr,avpure[v]-tau.pure[v],lty=3)
        title(main=paste("G=",g[i],", L=",load[j],sep=""),cex=0.8)
		
		# field data that correspond to this situation
		v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
		data <- spotweldfield$N[v]
		inputs <- spotweldfield$C[v]
		points(inputs,data)
#dev.off()
	}
}

for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avbias[v],type="l",ylim=c(-2,2),
			 xlab="current",ylab="weld diameter")
        lines(curr,biasL[v],lty=3)
        lines(curr,biasU[v],lty=3)
	}
}

for(i in 1:2){
	for(j in 1:2){
		v <- ((i-1)*40+(j-1)*20+1):((i-1)*40+j*20)
		plot(curr,avbc[v],type="l",ylim=c(3.5,8.5),
			 xlab="current",ylab="weld diameter")
        lines(curr,avbc[v]+taubc[v],lty=3)
        lines(curr,avbc[v]-taubc[v],lty=3)
		
# field data that correspond to this situation
		v <- ((i-1)*60+(j-1)*30+1):((i-1)*60+j*30)
		data <- spotweldfield$N[v]
		inputs <- spotweldfield$C[v]
		points(inputs,data)
	}
}