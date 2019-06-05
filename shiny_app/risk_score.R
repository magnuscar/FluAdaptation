########################################
######## function for shiny app ########
########################################
library(deSolve)



########################################
########          model         ########
########################################

model2 <- function(t, cinit, parms){

    with(as.list(c(cinit, parms)), {


		dU = - beta * U * (VE + VK)
	
		dIE = beta * U * VE - delta * IE
		dIK = beta * U * VK - delta * IK
		
		dVE = pE * IE - c * VE
		dVK = pK * IK - c * VK
		

       return(list(c(dU = dU, dIE = dIE, dIK = dIK, dVE = dVE, dVK = dVK)))
  })
}
 

# function to solve ODEs
solve_ode_model2 <- function(cinit, parms, t){

    # solve ODE for a given set of parameters
    out <- ode(y = cinit, times = t, parms = parms, func = model2, method="ode23")
    return(out)
}


predicting_passage_model2 <- function(cinit, parms, t, Npass, Tpass){
	
	# cinit 		vector with initial conditions, e.g. cinit<- c(U = 4e5, IE=0, IK=0, VE=0.8*400, VK=0.2*400 )
	# parms			vector with parameters parms <-  c(beta=2.7e-4,  delta=4, pE = 12/3, pK=12, c=3)
	# t				vector with time steps that need to be calculated t <- seq(0,5,0.1)
	# Npass			number of passages, Npas <- 5
	# Tpass			time (in steps) in simulation that corresponds to a passaging round Tpass<-10
	
	cinit_or <- cinit[ 1:3]
	
	# output matrix:  passage number and %VK
	out<- matrix(NA, nrow = Npass+1, ncol= 2, dimnames=list(rep("", Npass+1), c("passage", "percentPB2.627K")))
	
	
	# generate the first passage
	first <- solve_ode_model2(cinit, parms, t)
	
	# store the values in output matrix
	out[1,1] <- 0
	out[1,2] <- first[1,"VK"]/(first[1,"VK"]+first[1,"VE"]) 
	
	# determination for initial conditions in next round
	cinit_next <- c(cinit_or, first[1+Tpass,"VE"],  first[1+Tpass,"VK"])
	
	# loop over the further passages
	for( i in 1:Npass){
		
		zw <- solve_ode_model2(cinit_next, parms, t)
		
		out[1+i,1] <- i
		out[1+i,2] <- zw[1,"VK"]/(zw[1,"VK"]+zw[1,"VE"])
		
		
		# determination of new initial conditions
		cinit_next <- c(cinit_or, zw[1+Tpass,"VE"],  zw[1+Tpass,"VK"])
	}
	
	out[,2] <- 100*out[,2]
	out
}



########################################
########  calculate risk score  ########
########################################

normalised_risk_score <- function(predictions){
	
	# function that calculates the risk score of evolving a K mutation with range: [-1,1]
	# predictions: 		matrix with columns passage,  percentPB2.627K as generated with function predicting_passage_model2()
	startK <- predictions[which(predictions[,"passage"]==0) , "percentPB2.627K"]
	pred <- predictions[, "percentPB2.627K"] - startK	
	for(i in 1:(length(pred))){
		if(pred[i] < 0){ 
			pred[i] <- 	pred[i]/startK
		}
		else{
			pred[i] <- 	pred[i]/(100-startK)
		}
	}
	sum(pred)/ (dim(predictions)[1] -1)
}

risk_score <- function(fX1, fX3, perc.K = 20, beta = 2.7e-6, delta = 4, c = 3, v.pEi = c(200, 68.50662, 52.34761), v.pKi = c(136.60789, 94.48494, 200) ){
		
		parms <- c(beta = beta, delta = delta, c=c, pE = sum(c(fX1, 1- fX1-fX3, fX3)*v.pEi), pK = sum(c(fX1, 1- fX1-fX3, fX3)*v.pKi))
		normalised_risk_score( predicting_passage_model2(c(U = 4e5, IE=0, IK=0, VE=(1 - perc.K/100)*400, VK=perc.K/100*400 ), parms, t = seq(0,5,0.1), Npass = 5, Tpass = 6) )
	
}


########################################
########      color coding      ########
########################################

col.fred <- function(n, k) c(rev(rgb(0,0,1,alpha=seq(0,1,length=n/2)^k)),"white",rgb(1,0,0,alpha=seq(0,1,length=n/2)^(k)))



########################################
# calculating the riskscore threshold  #
########################################

return.change.riskscore <- function(fX1, data ){
	
	dat <- data[which(data[, "fX1"] == fX1),]
	
	nex <- dim(dat)[1]
	
	xout <- NA
	
	
	# 1st check: are all riskscores negative or positive?
	if( sum(dat[,"riskscore"] < 0) != nex  &  sum(dat[,"riskscore"] > 0) != nex ){
		if(dat[1,"riskscore"] < 0 ){
			zw <- dat[ tail(which(dat[, "riskscore"] < 0),1) + 0:1 , ]
		}
		if(dat[1,"riskscore"] > 0 ){
			zw <- dat[ tail(which(dat[, "riskscore"] > 0),1) + 0:1 , ]
		}
	
	xout <- - zw[1 ,"riskscore"] * (zw[2 ,"fX3"] - zw[1 ,"fX3"]) /(zw[2 ,"riskscore"] - zw[1 ,"riskscore"]) + zw[1 ,"fX3"]
	}
	xout
	
} 


return.change.riskscore.v <- function(FX1,data){ sapply(FX1, return.change.riskscore, data=data) }

   
get.0lines <- function(data){
    	
    	
	dat <- data[complete.cases(data), c("fX1", "fX3", "riskscore") ]
    	
   	FX1 <- unique(dat[,"fX1"])	
    	
    FX3 <- return.change.riskscore.v(FX1,dat)
    
    firstnonNA <- which(!is.na(FX3) )[1]
    
    if(firstnonNA > 1){
    	
    	FX3[1:(firstnonNA -1)]<- -0.2
    }
    
    	
    return(list(fX1= FX1, fX3 = FX3))
}
