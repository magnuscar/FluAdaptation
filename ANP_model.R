########################################################################
################## Virus Dynamics Model for Passaging ##################
########################################################################


###### Required libraries

library(deSolve) 


########################################################################
######                  Definition of functions                   ######
########################################################################

########################################################################
####      Case 1: beta, delta, c are the same for both viruses      ####
########################################################################

### Virus dynamics model

model2 <- function(t, cinit, parms){
  
  # definition of the ODE virus dynamics model 
  
  ### input:
  # t       time
  # cinit   vector of initial conditions, eg. cinit<- c(U = 4e5, IE = 0, IK = 0, VE = 0.8*400, VK = 0.2*400 )
  # parms   vector of parmeter values, eg. parms <- c(beta = 2.7e-6,  delta = 4, pE = 3, pK = 4, c = 3)
  ### output:
  # list with differentials of U, IE, IK, VE, VK
  
  with(as.list(c(cinit, parms)), {
 
    dU = - beta * U * (VE + VK)
    dIE = beta * U * VE - delta * IE
    dIK = beta * U * VK - delta * IK
    dVE = pE * IE - c * VE
    dVK = pK * IK - c * VK
    return(list(c(dU = dU, dIE = dIE, dIK = dIK, dVE = dVE, dVK = dVK)))
  })
}



### Solving the ODE

solve_ode_model2 <- function(cinit, parms, t){
  
  # solve ODE for a given set of parameters
  
  ### input:
  # cinit   vector of initial conditions, eg. cinit<- c(U = 4e5, IE = 0, IK = 0, VE = 0.8*400, VK = 0.2*400 )
  # parms   vector of parmeter values, eg. parms <- c(beta = 2.7e-6,  delta = 4, pE = 3, pK = 4, c = 3)
  # t       vector with time steps, eg. t <- seq(0,5,0.1)
  ### output:
  # matrix with time steps and the concentrations of U, IE, IK, VE, VK
  
  out <- ode(y = cinit, times = t, parms = parms, func = model2, method="ode23")
  return(out)
}



### Predicting virus abundance in passaging experiments

predicting_passage_model2 <- function(cinit, parms, t, Npass, Tpass){
  
  # function to predict percentPB2.627K (mammalian-adapted) virus percentage
  
  ### input: 
  # cinit   vector of initial conditions, eg. cinit<- c(U = 4e5, IE = 0, IK = 0, VE = 0.8*400, VK = 0.2*400 )
  # parms   vector of parmeter values, eg. parms <- c(beta = 2.7e-6,  delta = 4, pE = 3, pK = 4, c = 3)
  # t       vector with time steps, eg. t <- seq(0,5,0.1)
  # Npass		number of passages, eg. Npas <- 5
  # Tpass		time (in steps) in simulation that corresponds to a passaging round, eg. Tpass<-6
  ### output:
  # matrix of dimensions (Npass +1)x2 with passage number and percentPB2.627K
  
  
  
  # extract initial conditions for cells (Uj) and infected cells (Iij)
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



########################################################################
####      Case 2: beta, delta, c are different for each virus       ####
########################################################################

### Virus dynamics model

model2_diffratesforvir <- function(t, cinit, parms){
  
  # definition of the ODE virus dynamics model for different beta, delta, c for the different viral strains
  
  ### input:
  # t       time
  # cinit   vector of initial conditions, eg. cinit<- c(U = 4e5, IE = 0, IK = 0, VE = 0.8*400, VK = 0.2*400 )
  # parms   vector of parmeter values, eg. parms <- c(beta_E = 2.7e-6, beta_K = 2.7e-6, delta_E = 4, delta_K =4 , c_E = 3, c_K = 3, pE = 12, pK = 1)
  ### output:
  # list with differentials of U, IE, IK, VE, VK
  
  with(as.list(c(cinit, parms)), {
    dU = - beta_E * U * VE  - beta_K * U * VK
    dIE = beta_E * U * VE - delta_E * IE
    dIK = beta_K * U * VK - delta_K * IK
    dVE = pE * IE - c_E * VE
    dVK = pK * IK - c_K * VK
    return(list(c(dU = dU, dIE = dIE, dIK = dIK, dVE = dVE, dVK = dVK)))
  })
}



### Solving the ODE

solve_ode_model2_diffratesforvir <- function(cinit, parms, t){
  
  # solve ODE for a given set of parameters for different beta, delta, c for the different viral strains
  
  ### input:
  # cinit   vector of initial conditions, eg. cinit<- c(U = 4e5, IE = 0, IK = 0, VE = 0.8*400, VK = 0.2*400 )
  # parms   vector of parmeter values, eg. parms <- c(beta_E = 2.7e-6, beta_K = 2.7e-6, delta_E = 4, delta_K =4 , c_E = 3, c_K = 3, pE = 12, pK = 1)
  # t       vector with time steps, eg. t <- seq(0,5,0.1)
  ### output:
  # matrix with time steps and the concentrations of U, IE, IK, VE, VK
  
  out <- ode(y = cinit, times = t, parms = parms, func = model2_diffratesforvir, method="ode23")
  return(out)
}



### Predicting virus abundance in passaging experiments

predicting_passage_model2_diffratesforvir <- function(cinit = c(U = 4e5, IE=0, IK=0, VE=0.8*400, VK=0.2*400 ), parms, t, Npass, Tpass){
  
  # function to predict percentPB2.627K (mammalian-adapted) virus percentage for different beta, delta, c for the different viral strains
  
  ### input: 
  # cinit   vector of initial conditions, eg. cinit<- c(U = 4e5, IE = 0, IK = 0, VE = 0.8*400, VK = 0.2*400 )
  # parms   vector of parmeter values, eg. parms <- c(beta_E = 2.7e-6, beta_K = 2.7e-6, delta_E = 4, delta_K =4 , c_E = 3, c_K = 3, pE = 12, pK = 1)
  # t       vector with time steps, eg. t <- seq(0,5,0.1)
  # Npass		number of passages, eg. Npas <- 5
  # Tpass		time (in steps) in simulation that corresponds to a passaging round, eg. Tpass<-6
  ### output:
  # matrix of dimensions (Npass +1)x2 with passage number and percentPB2.627K
  
  # extract initial conditions for cells (Uj) and infected cells (Iij)
  
  cinit_or <- cinit[ 1:3]
  
  # output matrix:  passage number and %VK
  out<- matrix(NA, nrow = Npass+1, ncol= 2, dimnames=list(rep("", Npass+1), c("passage", "percentPB2.627K")))
  
  
  # generate the first passage
  first <- solve_ode_model2_diffratesforvir(cinit, parms, t)
  
  # store the values in output matrix
  out[1,1] <- 0
  out[1,2] <- first[1,"VK"]/(first[1,"VK"]+first[1,"VE"]) 
  
  # determination for initial conditions in next round
  cinit_next <- c(cinit_or, first[1+Tpass,"VE"],  first[1+Tpass,"VK"])
  
  # loop over the further passages
  for( i in 1:Npass){
    
    zw <- solve_ode_model2_diffratesforvir(cinit_next, parms, t)
    
    out[1+i,1] <- i
    out[1+i,2] <- zw[1,"VK"]/(zw[1,"VK"]+zw[1,"VE"])
    
    
    # determination of new initial conditions
    cinit_next <- c(cinit_or, zw[1+Tpass,"VE"],  zw[1+Tpass,"VK"])
  }
  
  out[,2] <- 100*out[,2]
  out
}


########################################################################
##########      Definition of the normalised risk score      ###########
########################################################################


normalised_risk_score <- function(predictions){
  
  # function that calculates the risk score of evolving a K mutation with range: [-1,1]
  
  ### input:
  # predictions: 		matrix with columns passage,  percentPB2.627K as generated with function predicting_passage_model2()
  ### output:
  # normalised risk score
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

