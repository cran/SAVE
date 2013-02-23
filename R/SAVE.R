##########################################################################
## SAVE main Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

`SAVE` <-
function (response.name=NULL, controllable.names=NULL,
			calibration.names=NULL,field.data=NULL,
            model.data=NULL,mean.formula=as.formula("~1"),bestguess=NULL){
             	
	model <- new ("SAVE")
	dprct <-.deprecate.parameters(call=sys.call(sys.parent(1)))
  	model@call <- as.call(dprct)

	#####
	#Field data:
  	if (is.null(field.data)){
  		bestguess <- NULL
  		model@yf <- NULL
  	  	if (!is.null(bestguess))
  			cat("The parameter bestguess unused since field data is not provided\n")
  		cat("The result wont have Stage II parameters.\n")
  		stop("Field data cannot be NULL\n")
  		}
	else {
		#Select those columns corresponding to the controllable inputs:
        yf <- as.vector(t(field.data[,response.name]))
        model@yf <- yf
	}	

	if (is.null(controllable.names))
		stop ("Controllable names cannot be NULL\n")
	else{
#df <- data.frame(field.data[,controllable.names])
		df <- as.matrix(field.data[,controllable.names])
		
    	colnames(df) <- controllable.names
		model@df <- df
	}
  	if (is.null(calibration.names)){
  		bestguess <- NULL
  		if (!is.null(bestguess))
  			cat("The parameter bestguess unused since calibartion.names parameter has not been provided\n")
	}
	else {
		if (is.null(bestguess))
			stop("bestguess cannot be NULL when there are calibration parameters.")
	}	
	
	model@controllablenames <- controllable.names
	model@responsename <- response.name
	model@calibrationnames <- calibration.names
	
    ####
	#Model data:
	dm <- as.matrix(model.data[,c(model@controllablenames,model@calibrationnames)])
	colnames(dm) <- c(model@controllablenames,model@calibrationnames)
	model@dm <- dm
    ym <- as.vector(t(model.data[,model@responsename]))
    model@ym <- ym

 	model@meanformula <- mean.formula 
    ####
    # Estimate stage I parameters using DiceKriging
    m <- km(model@meanformula,design=model@dm,response=model@ym,covtype="powexp")
    alphaM <- (m@covariance)@shape.val
    betaM <- ((m@covariance)@range.val)^(-alphaM)
    lambdaM <- 1/((m@covariance)@sd2)
    thtaM <- c(lambdaM,betaM,alphaM)
    names (thtaM) <- c("lambdaM",paste("betaM",c(model@controllablenames,model@calibrationnames),sep='.'),paste("alphaM",c(model@controllablenames,model@calibrationnames),sep='.'))
    thtaL <- m@trend.coef

    ####
    # Estimate stage II parameters

    # predict code at field design augmented with guess of u

    # produce design
		#takes the values on the list u.guess and puts on the right order (the one given by calibration.names)
		if (!is.null(bestguess)){
			model@bestguess<- vapply(bestguess, FUN=function(x){x[1]}, FUN.VALUE=c(0))[model@calibrationnames]
		}
		else model@bestguess<- NULL
				
				
		if (is.null(model@calibrationnames)){
			xnew<- model@df
		}
		else{
		    aux <- matrix(rep(model@bestguess,length(model@yf)),ncol=length(model@bestguess),byrow=T)
 			aux <- as.data.frame(aux)
			names(aux) <- model@calibrationnames
			xnew <- data.frame(model@df,aux)
			names(xnew) <- c(model@controllablenames,model@calibrationnames)
        }
    # predict code
    ymnew <- predict(m,newdata=xnew,type="UK",se.compute=FALSE,
                         light.return=TRUE)

    # resulting bias
    biashat <- model@yf-ymnew$mean
    # print(df)

    # estimate parameters
    # trend <- 0.0
    maux <- km(formula=~1,design=model@df,response=biashat,covtype="gauss",
                   # coef.trend=trend,
                   nugget.estim=TRUE)

    betab <- ((maux@covariance)@range.val)^(-2)/2
    lambdab <- 1/((maux@covariance)@sd2)
    lambdaF <- 1/((maux@covariance)@nugget)
    alphab <- rep(2.0,length(betab))

    thtaF <- c(lambdab,betab,alphab,lambdaF)
    names (thtaF) <- c("lambdaB",paste("betaB",model@controllablenames,sep='.'),paste("alphab",model@controllablenames,sep='.'),"lambdaF")

	model@mle <-list(thetaL=thtaL,thetaM=thtaM,thetaF=thtaF)
				
	
	## End of gaspfit
	## Beginning of bayesfit

	####################	
	#load some auxiliary functions
	duplicates <- function(dat)
	{
		s <- do.call("order", as.data.frame(dat))
        if(dim(as.matrix(dat))[2]==1){
          non.dup <- !duplicated(as.matrix(dat[s]))
        }
        else{
          non.dup <- !duplicated(as.matrix(dat[s,]))
        }
		orig.ind <- s[non.dup]
		first.occ <- orig.ind[cumsum(non.dup)]
		#first.occ[non.dup] <- NA
		first.occ[order(s)]
	}
	###############end of loading auxiliary funcs.	

	#####
	#Field data:
	#Select those columns corresponding to the controllable inputs:
	x<- field.data[,model@controllablenames]	
	#Extract the unique part of the field design matrix as required by bayesfit:
	x.unique<- as.data.frame(unique(x))
	names(x.unique) <- model@controllablenames
	my.original<- duplicates(x)
	yf.ordered<- numeric(0)
	nreps<- NULL
	for (i in unique(my.original)){
  		yf.ordered<- c(yf.ordered,model@yf[my.original==i])
  		nreps<- c(nreps, sum(my.original==i))
	}

	#establish the wd:
	wd <- tempdir()
	model@wd <- paste(wd,'/',sep='')	

	#Write to the files
	#file field_data.dat:
	write.table(yf.ordered, file=paste(model@wd,"/field_data.dat", sep=""),
                    col.names=F, row.names=F)
        
	#file field_nreps.dat:
	write.table(nreps, file=paste(model@wd,"/field_nreps.dat", sep=""),
                    col.names=F, row.names=F)
        
	#file field_unique.dat:
	write.table(x.unique, file=paste(model@wd,"/field_unique.dat",sep=""),
                    col.names=F, row.names=F)
	#####
	#####
	#Model data:
	x.m<- model.data[,c(model@controllablenames,model@calibrationnames)]
	#write to the files:
	write.table(x.m, file=paste(model@wd,"/model_inputs.dat",sep=""),
                    col.names=F, row.names=F)
	write.table(as.vector(t(model.data[,model@responsename])),
                    file=paste(model@wd,"/model_data.dat",sep=""),
                    col.names=F, row.names=F)
	#####
	#####
	#Mean response:
    model@xm <- model.matrix(model@meanformula, x.m)
				
	write.table(model@xm,
                    file=paste(model@wd,"/mcmc.field.design.M.matrix.dat",sep=""),
                    col.names=F, row.names=F)
	model@xf <- model.matrix(model@meanformula,x.unique)
   	write.table(model@xf,
                    file=paste(model@wd,"/mcmc.field.design.F.matrix.dat",sep=""),
                    col.names=F, row.names=F)
	#####
	dthetaM <- length(model@calibrationnames) + length(model@controllablenames)
	dthetaM <- dthetaM*2+1
	dthetaL <- ncol(model@xm)
	dthetaF <- 2*length(model@controllablenames)+2
	
	if(length(thtaM)!=dthetaM)
				{model@mle<- NULL; stop("Dimension of MLE for the covariance of the computer model is incorrect\n")}
	if(length(thtaL)!=dthetaL)
				{model@mle<- NULL; stop("Dimension of MLE for the mean of the computer model is incorrect\n")}
	if(length(thtaF)!=dthetaF)
				{model@mle<- NULL; stop("Dimension of MLE for the second stage parameters is incorrect\n")}
				
				

	#####
	#Write the MLE to the appropriate files
	write.table(model@mle$thetaM, file=paste(model@wd,"/thetaM_mle.dat", sep=""),
                    col.names=F, row.names=F)
	write.table(model@mle$thetaL, file=paste(model@wd,"/thetaL_mle.dat", sep=""),
                    col.names=F, row.names=F)
	write.table(model@mle$thetaF, file=paste(model@wd,"/thetaF_mle.dat", sep=""),
                    col.names=F, row.names=F)

#validObject(model, complete=TRUE)
	return(model)	
}

####################	
#load some auxiliary functions
normal<- function(var.name, mean, sd, lower, upper){
			c(var.name, 1, lower, upper, mean, sd^2)
}

uniform<- function(var.name, lower, upper){
			c(var.name, 0, lower, upper, 0, 0)
}

.expand.call <- function(call=sys.call(sys.parent(2)), expand.dots = TRUE)
{
	ans <- as.list(call)
	frmls <- formals(deparse(ans[[1]]))
	if (names(frmls[length(frmls)])=='...') frmls <- frmls[-length(frmls)]
	add <- which(!(names(frmls) %in% names(ans)))
	return(as.call(c(ans, frmls[add])))
}

.deprecate.parameters <- function(call=sys.call(sys.parent(1)), expand.dots = TRUE)
{
	ans <- as.list(.expand.call(call=call))
	frmls <- formals(deparse(ans[[1]]))
	aans <- ans [-1]
	if (names(frmls[length(frmls)])=='...') frmls <- frmls[-length(frmls)]
	#print(names(ans))
	#print(names(frmls))
	deprct <- which(!(names(aans) %in% names(frmls)))
	if (length(deprct) != 0) {
		cat('Warning!: Deprecated parameters:\n')
		print (as.character(names(aans)[deprct]))
		ans <- ans[-(deprct+1)]
	}
	#return (deprct)
	#return(call)
	#print(ans)
	return (ans)
}

