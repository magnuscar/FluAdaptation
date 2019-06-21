########################################################################
######                    Generating heatmaps                     ######
########################################################################

###### Required libraries

library(vioplot)
library(graphics)

###### set working directory
# setwd("/path/to/FluAdaptation/") ### <<<--- you have to add the path to the folder on your machine here and uncomment this line!!!

###### source required R scripts
source("ANP_model.R")
source("virus_production_rates/estimating_virus_production_rates.R")


########################################################################
######                  Definition of functions                   ######
########################################################################


list.with.PEK.estimates <- function(v.pEmin=0, v.pKmin=0, v.pEmax, v.pKmax, start.par=c(12,12), fn.data, fn.out, cinit = c(U=4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  #function that generates matrix with different starting values and their estimates for pEi and pKi
  
  ### input:
  # v.pEmin 			vector with minimal values for pE
  # v.pKmin			  vector with minimal values for pK
  # v.pEmax			  vector with maximal values for pE
  # v.pKmax			  vector with maximal values for pK
  # fn.dat				file name of X1, X2, X3 only expressing cells fn.data <- "data/input/ANP_passage_experiment.csv"
  # fn.out				file name of output file
  # remaining modelling parameters: see rss_for_estimating_pEonly()
  ### output
  # matrix with the lower and upper bounds for pEi and pKi, the estimates of pEi and pKi, and the residual sum of squares of between the data and the estimates in cells expressing Xi splice variant, i=1,2,3
  
  out <- matrix(NA, nrow=0, ncol=13, dimnames= list( c(), c("pEmin", "pEmax", "pKmin", "pKmax", paste0("pE", 1:3), paste0("pK", 1:3 ), paste0("valX", 1:3)) ) ) 
  write.table(out,fn.out,quote=F,row.names=F, sep =", ") 
  
  for(pEmin in v.pEmin){
    for(pKmin in v.pKmin){
      for(pEmax in v.pEmax){
        for(pKmax in v.pKmax){
          X1 <- estimate_pEpK_ANP("X1", start.par, fn.data, c(pEmin,pKmin), c(pEmax,pKmax), cinit, t, Npass, Tpass, beta, delta, c)
          X2 <- estimate_pEpK_ANP("X2", start.par, fn.data, c(pEmin,pKmin), c(pEmax,pKmax), cinit, t, Npass, Tpass, beta, delta, c)
          X3 <- estimate_pEpK_ANP("X3", start.par, fn.data, c(pEmin,pKmin), c(pEmax,pKmax), cinit, t, Npass, Tpass, beta, delta, c)
          zw <- matrix(c(pEmin, pEmax,pKmin,pKmax, X1$par[1], X2$par[1], X3$par[1], X1$par[2], X2$par[2], X3$par[2],  X1$value, X2$value, X3$value), nrow=1, ncol=13, dimnames= list( c(), c("pEmin", "pEmax", "pKmin", "pKmax", paste0("pE", 1:3), paste0("pK", 1:3 ), paste0("valX", 1:3)) ) ) 
          
          write.table(zw,fn.out,quote=F, append=T, col.names=F, row.names=F, sep =", ") 
          out <- rbind(out,zw)
        }
      }
    }
  }
  out
}



list.with.PEonly.estimates <- function(start.pE=12, v.pK, lower=0, upper=200, fn.data, fn.out, cinit = c(U=4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  #function that generates matrix with different starting values and their estimates for pEi and pKi, whereas pKi is fixed (and the same for X1, X2, X3)
  
  ### input:
  # start.pE			starting value of pE for fitting procedure
  # v.pK  		    vector with values for pK
  # fn.dat				file name of X1, X2, X3 only expressing cells fn.data <- "data/input/ANP_passage_experiment.csv"
  # fn.out				file name of output file
  # remaining modelling parameters: see rss_for_estimating_pEonly()
  ### output
  # matrix with the lower and upper bounds for pEi and pKi, the estimates of pEi and pKi, and the residual sum of squares of between the data and the estimates in cells expressing Xi splice variant, i=1,2,3
  
  out <- matrix(NA, nrow=0, ncol=13, dimnames= list( c(), c("pEmin", "pEmax", "pKmin", "pKmax", paste0("pE", 1:3), paste0("pK", 1:3 ), paste0("valX", 1:3) ) ) ) 
  write.table(out,fn.out,quote=F,row.names=F, sep =", ") 
  
  for(pK in v.pK){
    X1 <- estimate_pEonly_ANP("X1", start.pE, fn.data, pK, lower, upper, cinit, t, Npass, Tpass, beta, delta, c)
    X2 <- estimate_pEonly_ANP("X2", start.pE, fn.data, pK, lower, upper, cinit, t, Npass, Tpass, beta, delta, c)
    X3 <- estimate_pEonly_ANP("X3", start.pE, fn.data, pK, lower, upper, cinit, t, Npass, Tpass, beta, delta, c)
    zw <- matrix(c(lower, upper, pK, pK, X1$par, X2$par, X3$par, pK, pK, pK, X1$value, X2$value, X3$value), nrow=1, ncol=13, dimnames= list( c(), c("pEmin", "pEmax", "pKmin", "pKmax", paste0("pE", 1:3), paste0("pK", 1:3 ), paste0("valX", 1:3)) ) ) 
    
    write.table(zw,fn.out,quote=F, append=T, col.names=F, row.names=F, sep =", ") 
    out <- rbind(out,zw)
  }
  out
}



pEpK_passagingconc_riskscore <- function(fn.ANPratios, fac=2, fn.pEpK, v.startpercK, fn.out, U=4e5, IE=0, IK=0, VEK = 400, t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  # function to calculate riskscores for passage predictions in different animals

  ### input
  # fn.ANPratios		  file name of ANP ratios, fn.ANPratios <- "data/input/ANPratios.csv"
  # fac				        factor of how much the rss is allowed to deviate from the minimnal rss value
  # fn.pEpK			      file name of estimates of pE1-3 and pK1-3, fn.pEpK <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/pEpKestimates_varying_starting_conditions.csv"
  # v.startpercK			vector with starting values for K strain in risk analysis and passaging prediction plots
  # fn.out			      file name of where the outputted table should be stored
  # remaining parameters are the usual model parameters
  ### output
  # csv file with matrix, 20 columns with ("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" ), in the rows the different estimates for pE123, pK123
  
  ANPratios <- read.csv(fn.ANPratios)
  estpek <- read.csv(fn.pEpK)
  
  # reduction of all outcomes to those that have small deviance from the minimal RSS-values
  estpek <- estpek[ which(estpek[,"valX1"]< fac*min(estpek[,"valX1"]) & estpek[,"valX2"]< fac*min(estpek[,"valX2"]) & estpek[,"valX3"]<fac*min(estpek[,"valX3"])) , ]
 
  # table with animal passaging predictions and riskscore
  pprs <- matrix(NA, nrow=0, ncol= 20, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" )) )
  
  # definition of initial parameters
  cinit = c(U=U, IE=IE, IK=IK, VE= VEK, VK=VEK )
  
  for(animal in dimnames(ANPratios)[[1]]){
    for(startpercK in v.startpercK ){
      for(i in 1:(dim(estpek)[1])){
        zw.pprs <- matrix(c(animal,startpercK, ANPratios[animal,], estpek[i, c(paste0("pE", 1:3), paste0("pK", 1:3)) ],rep(NA,9)), nrow=1, ncol= 20, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" )) )
        
        zw.pE <- sum(ANPratios[animal,]/100*estpek[i, c(paste0("pE", 1:3))])
        zw.pK <- sum(ANPratios[animal,]/100*estpek[i, c(paste0("pK", 1:3))])
        zw.pprs[1, "pE"] <- zw.pE
        zw.pprs[1, "pK"] <- zw.pK
        zw.pprs[1, "ratiopEpK"] <- zw.pE/zw.pK
        
        cinit["VE"] <- (1 - startpercK/100)*VEK
        cinit["VK"] <- (startpercK/100)*VEK				
        parms <-  c(beta=beta,  delta=delta, pE = zw.pE, pK= zw.pK, c=c)
        zwpredictions <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
        zw.pprs[1, paste0("percKpass", 1:Npass)] <- zwpredictions[1 + 1:Npass,"percentPB2.627K"]
        
        zw.pprs[1,"riskscore"] <- normalised_risk_score(zwpredictions) 
        
        pprs <- rbind(pprs, zw.pprs)
      }
    }
  }	
  write.table(pprs, fn.out, quote=F, append=F, col.names=T, row.names=F, sep =", ")
}



reduction_passaging_estimates_blackbirdselE <- function(fn.pred, fn.out){
  
  # function to reduce the passaging predictions to only those parameters estimates for which blackbird selects for E
  
  ### input
  # fn.pred		file name of predictions for all parameter estimates, e.g. fn.pred <- "sensitivity_analyses/test/sensana1_riskscores.csv"
  # fn.out		file name of predicitons subset fn.out <- "sensitivity_analyses/test/sensana1_riskscores_blackbirdselE.csv"
  ### output
  # csv file with matrix with matrix, 20 columns with ("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" ), in the rows the different estimates for pE123, pK123
  #     which select for E in blackbirds
   
  pprs <- read.csv(fn.pred)
  #subset for blackbird
  pprsbb <- pprs[which(pprs[,"animal"]=="Blackbird"),]
  pprsbb <- pprsbb[which(pprsbb[,"riskscore"] <0 & pprsbb[,"startpercK"]==unique(pprsbb[,"startpercK"])[1] ), ] 
  
  pprsnew<-matrix(NA,nrow=0,ncol= (dim(pprsbb)[2]), dimnames=list(c(), dimnames(pprsbb)[[2]] ))
  
  for(i in 1:(dim(pprsbb)[1])){
    pprsnew <- rbind(pprsnew, pprs[ which(pprs[,"pE1"] == pprsbb[i,"pE1"] &
                                            pprs[,"pE2"] == pprsbb[i,"pE2"] &
                                            pprs[,"pE3"] == pprsbb[i,"pE3"] &
                                            pprs[,"pK1"] == pprsbb[i,"pK1"] &
                                            pprs[,"pK2"] == pprsbb[i,"pK2"] &
                                            pprs[,"pK3"] == pprsbb[i,"pK3"] 
    ), ] )		
  }
  
  pprsnew2 <- pprsnew[ order(pprsnew[,"animal"], pprsnew[,"startpercK"]), ]
 
  write.table(pprsnew2, fn.out, quote = F, sep = ", ", row.names=F)
}



violin_plot_riskscore <- function(fn.data, startpercK, fn.pdf){
  
  # function to generate violin plots of risk scores
  
  ### input
  # fn.data 		  file name of passaging data for all animals fn.data <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/estimation_of_passaging_diff_animals.csv"
  # startpercK		starting value for K strain in risk analysis 
  # fn.pdf		    file name of pdf file
  ### output 
  # pdf file with violin plot
  
  dat <- read.csv(fn.data)
  dat <- dat[which(dat[,"startpercK"]==startpercK), ]
  
  pdf(fn.pdf,width=12, height=12)
    vioplot(  as.numeric(dat[ which( dat[,"animal"] == "Mammals" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Swallow" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Magpie" ), "riskscore" ]),	
            as.numeric(dat[ which( dat[,"animal"] == "Blackbird" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Goose" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Swan" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Turkey" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Quail" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Chicken" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Duck" ), "riskscore" ]),
            as.numeric(dat[ which( dat[,"animal"] == "Gull" ), "riskscore" ]),
            col="lightgray",
            drawRect =F,
            ylim=c(-1,1),
            names = c("Mammals", "Swallow", "Magpie", "Blackbird", "Goose", "Swan", "Turkey", "Quail", "Chicken",  "Duck","Gull"),
            cex.lab = 0.75
    )
  
    xmax <- 12				
    xred <- c(0,xmax,xmax,0)
    yred <- c(0,0,1,1)
    xblue <-c(0,xmax,xmax,0)
    yblue <- c(0,0,-1,-1)
    
    polygon(xred, yred , col= rgb(1,0,0,0.1), border=NA)
    polygon(xblue, yblue, col= rgb(0,0,1,0.1), border=NA)
  dev.off()
}



estimate_globalpEpK_allANP_diffrates <- function(v.beta, v.delta, v.c, fn.data, fn.out, start.pEpK = rep(100,6), lower=rep(0,2), upper=rep(200,2), cinit = c(U = 4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6 ){
  
  # function to estimate all pEi and pKi values in case beta, delta and c are the same for the two viruses
  
  ### input
  # v.beta		  vector with values of beta that shall be tested
  # v.delta		  vector with death rates
  # v.c			    vector with clearance rates
  # fn.data     file name of data fn.data <- "data/input/ANP_passage_experiment.csv"
  # fn.out      file name of output file
  # start.pEpK  6 dimensional vector with starting values for pEi and pKi (i=1,2,3)
  # lower       vector with lower bounds for pEi and pKi
  # upper       vector with upper bounds for pEi and pKi 
  # remaining parameters are the parameters of the model
  ### output
  # csv file with matrix with 11 columns and one row per parameter combination with the estimation results
  
  out <- matrix(NA, nrow=0, ncol= 11, dimnames = list(c(),c("beta", "delta", "c", paste0("est_pE",1:3), paste0("est_pK",1:3), "sum_rss", "convergence")))
  write.table(out, fn.out, quote=F, sep =",", row.names=F, col.names=T)
  
  for(beta in v.beta){
    for(delta in v.delta){
      for(cvirion in v.c){
        
        zwischen1 <- estimate_pEpK_ANP("X1", c(start.pEpK[1], start.pEpK[4]), fn.data, lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass, beta = beta, delta = delta, c = cvirion)
        zwischen2 <- estimate_pEpK_ANP("X2", c(start.pEpK[2], start.pEpK[5]), fn.data, lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass, beta = beta, delta = delta, c = cvirion)
        zwischen3 <- estimate_pEpK_ANP("X3", c(start.pEpK[3], start.pEpK[6]), fn.data, lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass, beta = beta, delta = delta, c = cvirion)
        
        newrow <-  matrix(c(beta, delta, cvirion, zwischen1$par[1], zwischen2$par[1], zwischen3$par[1], zwischen1$par[2], zwischen2$par[2], zwischen3$par[2], zwischen1$value + zwischen2$value +zwischen3$value, zwischen1$convergence + zwischen2$convergence +zwischen3$convergence), nrow=1, ncol= 11, dimnames = list(c(),c("beta", "delta", "c", paste0("est_pE",1:3), paste0("est_pK",1:3), "sum_rss", "convergence")))
        
        write.table(newrow, fn.out, append =T, quote=F, sep =",", row.names=F, col.names=F)
        
        out <- rbind(out, newrow)
      }
    }
  }
}



estimate_allpEpK_diffratesforvir_diffrates <- function(v.betaE, v.deltaE, v.cE, v.betaK, v.deltaK, v.cK, fn.data, fn.out, start.pEpK = rep(100,6), lower=rep(0,2), upper=rep(200,2), cinit = c(U = 4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6 ){
  
  # function to estimate all pEi and pKi values in case beta, delta and c are differ for the two viruses
  
  ### input
  # v.betaE		  vector with infection rates with E-virus
  # v.deltaE	  vector with death rates of cells infected with E-virus 
  # v.cE			  vector with E-virus clearance rates
  # v.betaK		  vector with infection rates with K-virus
  # v.deltaK  	vector with death rates of cells infected with K-virus
  # v.cK		  	vector with K-virus clearance rates
  # fn.data     file name of data fn.data <- "data/input/ANP_passage_experiment.csv"
  # fn.out      file name of output file
  # start.pEpK  6 dimensional vector with starting values for pEi and pKi (i=1,2,3)
  # lower       vector with lower bounds for pEi and pKi
  # upper       vector with upper bounds for pEi and pKi 
  # remaining parameters are the parameters of the model
  ### output
  # csv file with matrix with 14 columns and one row per parameter combination with the estimation results
  
  out <- matrix(NA, nrow=0, ncol= 14, dimnames = list(c(),c("beta_E","beta_K", "delta_E", "delta_K", "c_E", "c_K", paste0("est_pE",1:3), paste0("est_pK",1:3), "sum_rss", "convergence")))
  write.table(out, fn.out, quote=F, sep =",", row.names=F, col.names=T)
  
  for(beta_E in v.betaE){
    for(beta_K in v.betaK){
      for(delta_E in v.deltaE){
        for(delta_K in v.deltaK){
          for(c_E in v.cE){
            for(c_K in v.cK){
              
              zwischen1 <- estimate_pEpK_ANP_diffratesforvir("X1", start.pEpK, fn.data, lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass,  beta_E = beta_E, beta_K = beta_K, delta_E = delta_E, delta_K = delta_K , c_E = c_E, c_K = c_K)
              zwischen2 <- estimate_pEpK_ANP_diffratesforvir("X2", start.pEpK, fn.data, lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass,  beta_E = beta_E, beta_K = beta_K, delta_E = delta_E, delta_K = delta_K , c_E = c_E, c_K = c_K)
              zwischen3 <- estimate_pEpK_ANP_diffratesforvir("X3", start.pEpK, fn.data, lower=lower, upper=upper, cinit = cinit, t = t, Npass = Npass, Tpass = Tpass,  beta_E = beta_E, beta_K = beta_K, delta_E = delta_E, delta_K = delta_K , c_E = c_E, c_K = c_K)
              
              newrow <-  matrix(c(beta_E, beta_K, delta_E, delta_K, c_E, c_K, zwischen1$par[1], zwischen2$par[1], zwischen3$par[1], zwischen1$par[2], zwischen2$par[2], zwischen3$par[2], zwischen1$value + zwischen2$value +zwischen3$value, zwischen1$convergence + zwischen2$convergence +zwischen3$convergence), nrow=1, ncol= 14, dimnames = list("", c("beta_E","beta_K", "delta_E", "delta_K", "c_E", "c_K", paste0("est_pE",1:3), paste0("est_pK",1:3), "sum_rss", "convergence")))
              
              write.table(newrow, fn.out, append =T, quote=F, sep =",", row.names=F, col.names=F)
              
              out <- rbind(out, newrow)
              
            }
          }
        }
      }
    }
  }	
}



pEpK_passagingconc_riskscore_diffrates <- function(fn.ANPratios, fac = 1.1, fn.pEpK, v.startpercK, fn.out, U = 4e5, IE=0, IK=0, VEK = 400, t = seq(0,5,0.1), Npass = 5, Tpass = 6){
  
  # function to calculate riskscores for passage predictions in different animals
  
  ### input
  # fn.ANPratios		  file name of ANP ratios, fn.ANPratios <- "data/input/ANPratios.csv"
  # fac				        factor of how much the rss is allowed to deviate from the minimnal rss value
  # fn.pEpK			      file name of estimates of pE1-3 and pK1-3, fn.pEpK <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/pEpKestimates_varying_starting_conditions.csv"
  # v.startpercK			vector with starting values for K strain in risk analysis and passaging prediction plots
  # fn.out			      file name of where the outputted table should be stored
  # remaining parameters are the usual model parameters
  ### output
  # csv file with matrix, 18 + Npass columns with ("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" ), in the rows the different estimates for pE123, pK123

  ANPratios <- read.csv(fn.ANPratios)
 
  # reduction of the estimates for which the rss is not bigger than fac * minimal rss 
  rest <- read.csv(fn.pEpK)
  rest <- rest[which( (rest[,"convergence"] == 0) & (rest[ , "sum_rss"] <= fac* min(rest[ , "sum_rss"], na.rm=T)) ) , ]
  
  # table with animal passaging predictions and riskscore
  pprs <- matrix(NA, nrow=0, ncol= 18+Npass, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), "beta", "delta", "c", paste0("est_pE", 1:3), paste0("est_pK", 1:3 ), "est_pE", "est_pK", "ratiopEpK", paste0("est_pErcKpass", 1:Npass), "riskscore" )) )
  
  # definition of initial parameters
  cinit = c(U=U, IE=IE, IK=IK, VE= VEK, VK=VEK )
  
  for(animal in dimnames(ANPratios)[[1]]){
    for(startpercK in v.startpercK ){
      for(i in 1:(dim(rest)[1])){

        v.est_pK <- rest[i, paste0("est_pK",1:3)]

        zw.pprs <- matrix(c(animal,startpercK, ANPratios[animal,],rep(NA,3), rest[i, paste0("est_pE", 1:3)], v.est_pK, rep(NA,4+Npass)), nrow=1, ncol= 18 + Npass, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), "beta", "delta", "c", paste0("est_pE", 1:3), paste0("est_pK", 1:3 ), "est_pE", "est_pK", "ratiopEpK", paste0("est_pErcKpass", 1:Npass), "riskscore" )) )
        
        zw.pE <- sum(ANPratios[animal,]/100*rest[i, c(paste0("est_pE", 1:3))])
        
        zw.pK <- ANPratios[animal,"X1"]/100 * rest[i,"est_pK1"] + ANPratios[animal,"X2"]/100 * rest[i,"est_pK2"] + ANPratios[animal,"X3"]/100 * rest[i,"est_pK3"]

        zw.pprs[1, "est_pE"] <- zw.pE
        zw.pprs[1, "est_pK"] <- zw.pK
        zw.pprs[1, "ratiopEpK"] <- zw.pE/zw.pK
        
        cinit["VE"] <- (1 - startpercK/100)*VEK
        cinit["VK"] <- (startpercK/100)*VEK				
        parms <-  c(beta = as.numeric(rest[i,"beta"]),  delta = as.numeric(rest[i,"delta"]), pE = as.numeric(zw.pE), pK = as.numeric(zw.pK), c = as.numeric(rest[i,"c"]))					
        
        zw.pprs[1,"beta"] <- as.numeric(rest[i,"beta"])
        zw.pprs[1, "delta"] <- as.numeric(rest[i,"delta"])
        zw.pprs[1, "c"] <- as.numeric(rest[i,"c"])
        
        zwpredictions <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
        
        zw.pprs[1, paste0("est_pErcKpass", 1:Npass)] <- zwpredictions[1 + 1:Npass,"percentPB2.627K"]
        
        zw.pprs[1,"riskscore"] <- normalised_risk_score(zwpredictions) 
        
        pprs <- rbind(pprs, zw.pprs)
      }
    }
  }	
  write.table(pprs, fn.out, quote=F, append=F, col.names=T, row.names=F, sep =", ")
}



pEpK_passagingconc_riskscore_diffrates_diffbetadeltac <- function(fn.ANPratios, fac, fn.pEpK, v.startpercK, fn.out, U = 4e5, IE=0, IK=0, VEK = 400, t = seq(0,5,0.1), Npass = 5, Tpass = 6){
  
  # function to calculate riskscores for passage predictions in different animals with varying beta, delta, c for the two viral strains
  
  ### input
  # fn.ANPratios		  file name of ANP ratios, fn.ANPratios <- "data/input/ANPratios.csv"
  # fac				        factor of how much the rss is allowed to deviate from the minimnal rss value
  # fn.pEpK			      file name of estimates of pE1-3 and pK1-3, fn.pEpK <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/pEpKestimates_varying_starting_conditions.csv"
  # v.startpercK			vector with starting values for K strain in risk analysis and passaging prediction plots
  # fn.out			      file name of where the outputted table should be stored
  # remaining parameters are the usual model parameters
  ### output
  # csv file with matrix, 21 + Npass columns with ("animal","startpercK", paste0("X", 1:3), "beta_E", "beta_K", "delta_E", "delta_K", "c_E", "c_K", paste0("est_pE", 1:3), paste0("est_pK", 1:3 ), "est_pE", "est_pK", "ratiopEpK", paste0("est_pErcKpass", 1:Npass), "riskscore" )), in the rows the different estimates for pE123, pK123
  
  ANPratios <- read.csv(fn.ANPratios)
  
  # reduction of the estimates for which the rss is not bigger than fac * minimal rss 
  rest <- read.csv(fn.pEpK)
  rest <- rest[which( (rest[,"convergence"] == 0) & (rest[ , "sum_rss"] <= fac* min(rest[ , "sum_rss"], na.rm=T)) ) , ]
  
  # definition of initial parameters
  cinit = c(U=U, IE=IE, IK=IK, VE= VEK, VK=VEK )
  
  # table with animal passaging predictions and riskscore
  pprs <- matrix(NA, nrow=0, ncol= 21 + Npass, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), "beta_E", "beta_K", "delta_E", "delta_K", "c_E", "c_K", paste0("est_pE", 1:3), paste0("est_pK", 1:3 ), "est_pE", "est_pK", "ratiopEpK", paste0("est_pErcKpass", 1:Npass), "riskscore" )) )
  write.table(pprs, fn.out, quote=F, append=F, col.names=T, row.names=F, sep =", ")
  
  
  for(animal in dimnames(ANPratios)[[1]]){
    for(startpercK in v.startpercK ){
      for(i in 1:(dim(rest)[1])){
        
        v.est_pK <- rest[i, paste0("est_pK",1:3)]
        
        
        zw.pprs <- matrix(c(animal,startpercK, ANPratios[animal,],rep(NA,6), rest[i, paste0("est_pE", 1:3)], v.est_pK, rep(NA,4+Npass)), nrow=1, ncol= 21 + Npass, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), "beta_E", "beta_K", "delta_E", "delta_K", "c_E", "c_K", paste0("est_pE", 1:3), paste0("est_pK", 1:3 ), "est_pE", "est_pK", "ratiopEpK", paste0("est_pErcKpass", 1:Npass), "riskscore" )) )
        
        zw.pE <- sum(ANPratios[animal,]/100*rest[i, c(paste0("est_pE", 1:3))])	
        zw.pK <- ANPratios[animal,"X1"]/100 * rest[i,"est_pK1"] + ANPratios[animal,"X2"]/100 * rest[i,"est_pK2"] + ANPratios[animal,"X3"]/100 * rest[i,"est_pK3"]
        
        
        zw.pprs[1, "est_pE"] <- zw.pE
        zw.pprs[1, "est_pK"] <- zw.pK
        zw.pprs[1, "ratiopEpK"] <- zw.pE/zw.pK
        
        cinit["VE"] <- (1 - startpercK/100)*400
        cinit["VK"] <- (startpercK/100)*400				
        
        #parms <-  c(beta = as.numeric(rest[i,"beta"]),  delta = as.numeric(rest[i,"delta"]), pE = as.numeric(zw.pE), pK = as.numeric(zw.pK), c = as.numeric(rest[i,"c"]))		
        parms <- c(beta_E = as.numeric(rest[i,"beta_E"]), beta_K = as.numeric(rest[i,"beta_K"]), delta_E = as.numeric(rest[i,"delta_E"]), delta_K = as.numeric(rest[i,"delta_K"]), c_E = as.numeric(rest[i,"c_E"]), c_K = as.numeric(rest[i,"c_K"]),  pE = as.numeric(zw.pE), pK = as.numeric(zw.pK))
        
        zwpredictions <- predicting_passage_model2_diffratesforvir(cinit, parms, t, Npass, Tpass)
        
        zw.pprs[1,"beta_E"] <- as.numeric(rest[i,"beta_E"])
        zw.pprs[1, "delta_E"] <- as.numeric(rest[i,"delta_E"])
        zw.pprs[1, "c_E"] <- as.numeric(rest[i,"c_E"])
        zw.pprs[1,"beta_K"] <- as.numeric(rest[i,"beta_K"])
        zw.pprs[1, "delta_K"] <- as.numeric(rest[i,"delta_K"])
        zw.pprs[1, "c_K"] <- as.numeric(rest[i,"c_K"])
        
        zw.pprs[1, paste0("est_pErcKpass", 1:Npass)] <- zwpredictions[1 + 1:Npass,"percentPB2.627K"]
        
        zw.pprs[1,"riskscore"] <- normalised_risk_score(zwpredictions) 
        
        write.table(zw.pprs, fn.out, quote=F, append=T, col.names=F, row.names=F, sep =", ")
        pprs <- rbind(pprs, zw.pprs)
      }
    }
  }	
}



###################################################
######       Example for function call       ######
###################################################

## Note: this section is commented out such that this file can be sourced without producing any output.

# # source("sensitivity_analyses/sensitivityanalyses.R")
#
# fol.out <- "sensitivity_analyses/test/"
# system(paste0("mkdir ", fol.out))
#
# ### Sensitivity analysis 1: testing different assumptions on the range of the virus production rates pEi and pKi (SuppFigure 3a)
# # (i) generate a list of pEi estimates for fixed pKi (which ranges between 1 and 200) WARNING: this takes about 10min to run
# start.pE <-12
# v.pK <- 1:200
# lower <- 0
# upper <- 200
# fn.data <- "data/input/ANP_passage_experiment.csv"
# fn.out <- paste0(fol.out, "sensana1_estimatespEpK_pKconst.csv")
# # for the remaining parameters we used the pre-defined ones
# list.with.PEonly.estimates(start.pE, v.pK, lower, upper, fn.data, fn.out)
# # (ii)  generate a list of pEi and pKi estimates WARNING: this takes about 1h to run
# v.pEmin <- 0
# v.pKmin <- 0
# v.pEmax <- seq(5, 200, by=5)
# v.pKmax <- seq(5, 200, by=5)
# start.par <- c(4,4)
# fn.data <- "data/input/ANP_passage_experiment.csv"
# fn.out <- paste0(fol.out, "sensana1_estimatespEpK_pKvar.csv")
# list.with.PEK.estimates(v.pEmin, v.pKmin, v.pEmax, v.pKmax, start.par, fn.data, fn.out)
# # (iii) combine the two outputted tables from (i) and (ii), the combined table we generated can be found under "data/results/sensana1_pEpKestimates.csv"
# # (iv)  calculate the riskscores for different animals
# fn.ANPratios <- "data/input/ANPratios.csv"
# fac <- 2
# fn.pEpK <- "data/results/sensitivity_analyses/sensana1_pEpKestimates.csv"
# v.startpercK <- 1
# fn.out <-  paste0(fol.out, "sensana1_riskscores.csv")
# pEpK_passagingconc_riskscore(fn.ANPratios, fac, fn.pEpK, v.startpercK, fn.out)
# # (v)   reduction of pE pK estimates to those that select for E in blackbirds
# fn.pred <-  paste0(fol.out, "sensana1_riskscores.csv")
# fn.out <-  paste0(fol.out, "sensana1_riskscores_blackbirdselE.csv")
# reduction_passaging_estimates_blackbirdselE(fn.pred, fn.out)
# # (vi)  generate violinplots 
# fn.data <-  paste0(fol.out, "sensana1_riskscores_blackbirdselE.csv")
# startpercK <-1
# fn.pdf <- paste0(fol.out, "sensana1_violinplot_riskscores.pdf")
# violin_plot_riskscore(fn.data, startpercK, fn.pdf)
# #
# # ### Sensitivity analysis 2: testing the influence of the parameters beta, delta, c when they are varying but are the same for the two variants (SupplFig 3b) 
# # 
# # (i) create the list with pEi, pKi estimates !!! This function is very slow. The table we created with this function call can be found under "data/results/sensana2_pEpKestimates.csv"
# v.beta <- seq(1e-6, 9e-5, length.out = 10)
# v.delta <- 1:10
# v.c <- 1:10
# fn.data <- "data/input/ANP_passage_experiment.csv"
# fn.out <- paste0(fol.out, "sensana2_pEpKestimates.csv")
# estimate_globalpEpK_allANP_diffrates(v.beta, v.delta, v.c, fn.data, fn.out)
# # (ii) calculate the riskscores for different animals
# fn.ANPratios <- "data/input/ANPratios.csv"
# fac <- 1.1
# fn.pEpK <- "data/results/sensitivity_analyses/sensana2_pEpKestimates.csv"
# v.startpercK <- 1
# fn.out <-  paste0(fol.out, "sensana2_riskscores.csv")
# pEpK_passagingconc_riskscore_diffrates(fn.ANPratios, fac, fn.pEpK, v.startpercK, fn.out)
# # (iii) violin plots
# fn.data <- paste0(fol.out, "sensana2_riskscores.csv")
# startpercK <- 1
# fn.pdf <- paste0(fol.out, "sensana2_violinplot_riskscores.pdf")
# violin_plot_riskscore(fn.data, startpercK, fn.pdf)
# #
# ### Sensitivity analysis 3: testing the influence of beta, delta, c, when they can even differ for the two variants (SupplFig 3c)
# #
# # (i) create the list with pEi, pKi estimates !!! This function is very slow. The table we created with this function call can be found under "data/results/sensana3_pEpKestimates.csv"
# v.betaE <- c(2.7e-6, seq(8.8e-7,1e-5, length.out = 3))
# v.betaK <- v.betaE
# v.deltaE <- c(4, seq(2.6, 6.1, length.out = 3))
# v.deltaK <- v.deltaE
# v.cE <- c(2.4, 2.8, 3, 3.6)
# v.cK <- v.cE
# fn.data <- "data/input/ANP_passage_experiment.csv"
# fn.out <- paste0(fol.out, "sensana3_pEpKestimates.csv")
# estimate_allpEpK_diffratesforvir_diffrates(v.betaE, v.deltaE, v.cE, v.betaK, v.deltaK, v.cK, fn.data, fn.out)
# # (ii) calculate the riskscores for different animals
# fn.ANPratios <- "data/input/ANPratios.csv"
# fac <- 1.1
# fn.pEpK <- "data/results/sensitivity_analyses/sensana3_pEpKestimates.csv"
# v.startpercK <- 1
# fn.out <-  paste0(fol.out, "sensana3_riskscores.csv")
# pEpK_passagingconc_riskscore_diffrates_diffbetadeltac(fn.ANPratios, fac, fn.pEpK, v.startpercK, fn.out)
# # (iii) 
# fn.data <- paste0(fol.out, "sensana3_riskscores.csv")
# startpercK <- 1
# fn.pdf <- paste0(fol.out, "sensana3_violinplot_riskscores.pdf")
# violin_plot_riskscore(fn.data, startpercK, fn.pdf)

