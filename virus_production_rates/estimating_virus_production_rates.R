#######################################################################
##########         Estimating virus production rates         ##########
#######################################################################

###### Required libraries

library()


###### set working directory
# setwd("/path/to/FluAdaptation/") ### <<<--- you have to add the path to the folder on your machine here and uncomment this line!!!

###### source required R scripts
source("ANP_model.R")

########################################################################
######                  Definition of functions                   ######
########################################################################

########################################################################
####      Case 1: beta, delta, c are the same for both viruses      ####
########################################################################

rss_for_estimating_pEpK <- function(pEpK, data, cinit = c(U=4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3 ){
  
  # function to calculate the residual sum of square between an estimated passaging outcome and the real data
  
  ### input:
  # pEpK		  vector with virus production rates of cells infected with E- and K-virus: pEpK <- c(pE, pK)
  # data		  passage experiment data for one ANP ratio, eg.  data <- subset(read.csv("data/input/ANP_passage_experiment.csv"), ANP=="X1")
  # cinit     vector with initial conditions, the ratio between E and K virus is the mean of the data points at t=0
  # t         vector with time steps
  # Npass		  number of passages
  # Tpass		  time (in steps) in simulation that corresponds to a passaging round
  # beta      infection rate
  # delta     infected cell death rate
  # c         death rate of virus
  ### output:
  # one value for the rss
  
  parms <-  c(beta = beta,  delta = delta, pE = pEpK[1], pK=pEpK[2], c=c)
  zwischen <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
  sum((data[,"percentK"] - zwischen[data[,"passage"]+1,"percentPB2.627K"])^2 )
}



estimate_pEpK_ANP <- function(ANP, start.par = c(20,20), fn.data, lower = c(0,0), upper = c(200,200), cinit = c(U=4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3 ){
  
  # function to estimates the virus production rates 
  
  ### input:
  # ANP		    character string of the cell line for which the virus production rates should be estimated, either "X1", "X2" or "X3"
  # start.par vector with starting values for pE and pK; we tested many different starting values which did not have an influence on the final estimates
  # fn.data	  file name of data file, ie. fn.data <- "data/input/ANP_passage_experiment.csv"
  # lower     vector with lower bounds for pE and pK
  # upper     vector with upper bounds for pE and pK
  # remaining modelling parameters see rss_for_estimating_pEpK()
  ### output:
  # list created by optim(), for further information on the objectes contained in this list call help(optim)
  
  anp <- ANP
  data <- subset(read.csv(fn.data), ANP == anp)
  estimate <- optim(start.par, rss_for_estimating_pEpK, data = data, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass, beta = beta, delta = delta, c = c, method="L-BFGS-B", lower=lower,upper=upper)
  estimate
  
}



rss_for_estimating_pEonly <- function(pE, pK=15, data, cinit = c(U=4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  # function to calculate the residual sum of square for a fixed pK value, i.e. the same K-virus production for each cell type expressing one of the 3 ANP splice variants
  # used for sensitivity analysis (see sensitivity_analyses/sensitivityanalyses.R)
  
  ### input:
  # pE, pK		burst sizes of E- and K- viruses, respectively
  # data		  passage experiment data for one ANP ratio, eg.  data <- subset(read.csv("data/input/ANP_passage_experiment.csv"), ANP=="X1")
  # cinit     vector with initial conditions, the ratio between E and K virus is the mean of the data points at t=0
  # t         vector with time steps
  # Npass		  number of passages
  # Tpass		  time (in steps) in simulation that corresponds to a passaging round
  # beta      infection rate
  # delta     infected cell death rate
  # c         death rate of virus
  ### output:
  # one value for the rss
  
  parms <-  c(beta = beta,  delta = delta, pE = pE, pK = pK, c = c)
  zwischen <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
  sum((data[,"percentK"] - zwischen[,"percentPB2.627K"])^2)
}



estimate_pEonly_ANP <- function(ANP, start.pE, fn.data, pK = 15, lower = 0, upper = 200, cinit = c(U=4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  # function to estimate the E-virus production rate based on the experimental passaging data
  # used for sensitivity analysis (see sensitivity_analyses/sensitivityanalyses.R)
  
  ### input:
  # ANP		    character string of the cell line for which the virus production rates should be estimated, either "X1", "X2" or "X3"
  # start.pE  starting value for pE
  # fn.data	  file name of data file, ie. fn.data <- "data/input/ANP_passage_experiment.csv"
  # pK        value for pK
  # lower     vector with lower bound for pE
  # upper     vector with upper bound for pE
  # remaining modelling parameters: see rss_for_estimating_pEonly()
  ### output:
  # list created by optim(), for further information on the objectes contained in this list call help(optim)
  
  bla <- ANP
  data <- subset(read.csv(fn.data), ANP==bla)
  
  estimate <- optim(start.pE, rss_for_estimating_pEonly, pK=pK, data = data, method="L-BFGS-B", lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass, beta = beta, delta = delta, c = c)
  return(estimate)
}




### function for boostrap procedure

estimate.pEi.pKi.bootstrap <- function(fn.data, ANP, fn.intdat, fn.out, start.par = c(20,20), lower= c(0,0), upper = c(200,200),  cinit = c(U = 4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  # this function subsamples the original data and estimates the virus production rates for this new data, which corresponds to one bootstrap step
  
  ### input:
  # fn.data		file name of data fn.data <- "data/input/ANP_passage_experiment.csv"
  # ANP		    character string of the cell line for which the virus production rates should be estimated, either "X1", "X2" or "X3"
  # fn.intdat	file name where the subsampled data is stored
  # fn.out    file name of where the output matrix should be stored
  # remaining modelling parameters see above
  ### output
  # matrix with estimates of the virus production rates in different cell types, the rss of this specific estimation procdure, the convergence indicator as defined for optim(), and the percentage of K virus for the passaging steps
  
  dat <- read.csv(fn.data)
  dat <- dat[which(dat[,"ANP"] == ANP),]
  new_rows <- sample(1: (dim(dat)[1]), dim(dat)[1], replace = T )
  dat <- dat[new_rows, ]
  write.table(dat, fn.intdat, append = F, quote = F, sep =",", row.names=F)
  
  zwischen <- estimate_pEpK_ANP(ANP, start.par, fn.intdat, lower, upper, cinit,  t , Npass , Tpass, beta, delta, c )
  parms <-  c(beta = beta,  delta = delta, pE = zwischen$par[1], pK = zwischen$par[2], c = c)
  v.pred <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
  out <- matrix(c(zwischen$par, zwischen$value, zwischen$convergence, v.pred[,"percentPB2.627K"]),
                nrow = 1, ncol = 10 , dimnames= list("", c(paste0(c("est_pE", "est_pK"), ANP), "rss", "convergence", paste0("percK_pass", 0:Npass) ) ) )
  
  system(paste0("rm ", fn.intdat))
  
  write.table(out, fn.out, append = F, quote = F, sep =",", row.names=F)
}



####################################################
######       Example for function calls       ######
####################################################

# # Note this section is commented out such that this file can be sourced without producing any output.
# 
# # source("virus_production_rates/estimating_virus_production.R")
# # To find the best for estimates for the pre-defined parameters call:
# ANP <- "X1"     # or "X2" or "X3"
# fn.data <- "data/input/ANP_passage_experiment.csv"
# estimate_pEpK_ANP(ANP = ANP, fn.data = fn.data)
# 
# # For the bootstrapping procedure, we called the function estimate.pEi.pKi.bootstrap() 1000 times for each of "X1", "X2", "X3"
# # for easy use, we combined these bootstrap estimates in the files "data/results/estimation_virus_production/CIestimation_X<1,2,3>.csv"
# # here we show one bootstrap call for the pre-defined parameters:
# fn.data <- "data/input/ANP_passage_experiment.csv"
# ANP <- "X1"     # or "X2" or "X3"
# # define the output folder and create this folder:
# fol.out <- "virus_production_rates/test/"
# system(paste0("mkdir ", fol.out))
# fn.intdat <- paste0(fol.out,"intermediatedata.csv")
# fn.out <- paste0(fol.out,"bootstrap_", ANP, "_run1.csv")
# estimate.pEi.pKi.bootstrap(fn.data, ANP, fn.intdat, fn.out)
#
# # From the obtained bootstrap estimates, we determined the 95% CI intervals for the separate parameter estimates as explained in "CIbootstrap.Rmd"
# # A table with the best estimates and the 95% CI is stored in "data/results/estimation_virus_production/bestfits_pEi_pKi_withbootstrapCI.csv"



