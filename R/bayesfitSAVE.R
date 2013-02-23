##########################################################################
## Bayesfit Method
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################
bayesfit.SAVE <- function (object, prior,mcmcMultmle=1,
		     prob.prop=0.5, method=2, n.iter, nMH=20, n.burnin=0,
         	 n.thin=1){
    	 	
	# To deprecate the unused parameters and to include in the 
	# call() all the default parameters not used in the call
	# to the function
	dprct <- .deprecate.parameters(call=sys.call(sys.parent(1)))
	object@bayesfitcall <- as.call(dprct)
	
	#####
	#Write the file for the prior
	if (!is.null(object@calibrationnames)){
	if ((length(prior)/6)<length(object@calibrationnames))
          {stop("Not enough prior distributions\n")}
	if ((length(prior)/6)>length(object@calibrationnames))
          {stop("Too many prior distributions\n")}
	prior.matrix<- matrix(0, ncol=length(object@calibrationnames), nrow=5)
	colnames(prior.matrix)<- object@calibrationnames
	for (i.name in object@calibrationnames){
		thisplace<- which(prior==i.name)
		prior.matrix[,i.name]<- as.numeric(prior[thisplace+1:5])
	}
	object@prior <- t(prior.matrix)
	write.table(object@prior, file=paste(object@wd,"/bounds.dat",sep=""),
                    col.names=F, row.names=F)
	}

	#####
	#Values of the parameters:
	printOutput<- 0
	#total number of inputs:
	numInputs<- length(c(object@controllablenames,object@calibrationnames))
	#number of calibration inputs:
	numCalibration<- length(object@calibrationnames)
	#dimension of the linear model for the mean of the GP prior
	#        (q in C notation)
	numPModel<- dim(object@xm)[2]
	#dimension of code design (NM in C notation)
	sizeData<- length(object@ym)
#dimension of field design unique inputs (NF in Rui''s notation)
	sizeField<- dim(object@xf)[1]
#Probability of sampling from the prior
	object@method <- method
	object@mcmcMultmle <- mcmcMultmle

	# load the C code
	if(!is.loaded("bayesfit")) {
		lib.file <- file.path(paste("bayesfit",
                                            .Platform$dynlib.ext, sep=""))
		dyn.load(lib.file)
		cat(" -Loaded ", lib.file, "\n")
	}

    #Call to the fucnction:
	if ((object@method !=1) && (object@method != 2)){
		stop ("Wrong type of method introduced as parameter")
	}
	else
	output <- .C('bayesfit',as.integer(printOutput),
                     as.integer(numInputs),as.integer(numCalibration),
                     as.integer(numPModel),as.integer(sizeData),
                     as.integer(sizeField),as.double(prob.prop),
                     as.double(object@mcmcMultmle),as.integer(object@method),
                     as.integer(n.iter),as.integer(nMH),
                     as.character(object@wd))
#cat('The results can be found on ',object@wd,'\n')
#	system(paste('ls ',object@wd,'*.out',sep=''))

	#info:
	cat("Acceptance rate:",scan(file=paste(object@wd,"rate.out",sep="/"), 
			quiet=T)[1],"\n")
     	 	
    #//////
	post.thetaF<- read.table(file=paste(object@wd,"thetaF.out",sep="/"), header=F, 
                col.names=c("lambdaB","lambdaF"))
	
	if (!is.null(object@calibrationnames)){
		post.calparams<- read.table(file=paste(object@wd,"ustar.out",sep="/"), 
			 header=F, col.names=object@calibrationnames)
		auxparams<- cbind(post.calparams,post.thetaF)
	}
	else auxparams<- post.thetaF
    #burnin and thining:
    auxparams<- auxparams[-(1:n.burnin),]
	auxparams<- auxparams[seq(from=1,to=dim(auxparams)[1], by=n.thin),]
	row.names(auxparams)<- 1:(dim(auxparams)[1])
	object@mcmcsample <- as.matrix(auxparams)
	

	object
}

if(!isGeneric("bayesfit")) {
  setGeneric(name = "bayesfit",
             def = function(object, prior,mcmcMultmle=1,
		     prob.prop=0.5, method=2, n.iter, nMH=20, n.burnin=0,
         	 n.thin=1,...) standardGeneric("bayesfit")
             )
}

setMethod("bayesfit", "SAVE", 
          function(object,prior, n.iter,...) {
			 #print(match.call(expand.dots=T))
			 #deprctcall <- .deprecate.parameters()
			 #print(as.call(deprctcall))
			 #object@bayesfitcall <- as.call(deprctcall)
			 bayesfit.SAVE(object=object, prior=prior,
			 n.iter = n.iter, mcmcMultmle=mcmcMultmle,
		     prob.prop=prob.prop, method=method, nMH=nMH, n.burnin=n.burnin,
         	 n.thin=n.thin)
          }
          )
          

