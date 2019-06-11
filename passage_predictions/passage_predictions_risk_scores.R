#######################################################################
##########       Predicting passage outcomes in species      ##########
#######################################################################

###### Required libraries

library(graphics)


###### set working directory
# setwd("/path/to/FluAdaptation/") ### <<<--- you have to add the path to the folder on your machine here and uncomment this line!!!

###### source required R scripts
source("ANP_model.R")



########################################################################
######                  Definition of functions                   ######
########################################################################

host_species_passaging_predictions_CIwithbootstrap <- function(fn.estCI, fn.ANP, fn.out, v.startpercK, Npass = 5, U = 4e5, IE=0, IK=0, sumVEVK=400, t = seq(0,5,0.1), Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  # function to plot passage predictions
  
  ### input
  # fn.estCI	    file name with best estimates from passage experiments and CI fn.estCI <- "data/results/estimation_virus_production/bestfits_pEi_pKi_withbootstrapCI.csv"
  # fn.ANP		    file name of ANP ratios in species, fn.ANP <- "data/input/ANPratios.csv"
  # fn.out		    start of the file name for the outputed .pdf; the animal name and .pdf are automatically attached
  # v.startpercK	vector with starting percentages of K-virus
  # remaining variables: model parameters
  ### output
  # one .pdf with passage predictions for each species in fn.ANP and for all starting percentages of K-virus specified in v.startpercK
  
  est <- read.csv(fn.estCI)
  anp <- read.csv(fn.ANP)
  
  for(animal in dimnames(anp)[[1]]){
    
    fn.pdf <- paste0(fn.out, animal, ".pdf")
    
    ratioXi <- as.numeric(anp[animal,])
    
    v.pE_lower <- est[which(est[,"kind"] == "lower"), "pEi" ]
    v.pK_lower <- est[which(est[,"kind"] == "lower"), "pKi" ]
    pE_lower <- sum( v.pE_lower * ratioXi/sum(ratioXi) )
    pK_lower <- sum( v.pK_lower * ratioXi/sum(ratioXi) )
    
    v.pE_estimate <- est[which(est[,"kind"] == "estimate"), "pEi" ]
    v.pK_estimate <- est[which(est[,"kind"] == "estimate"), "pKi" ]
    pE_estimate <- sum( v.pE_estimate * ratioXi/sum(ratioXi) )
    pK_estimate <- sum( v.pK_estimate * ratioXi/sum(ratioXi) )	
    
    v.pE_upper <- est[which(est[,"kind"] == "upper"), "pEi" ]
    v.pK_upper <- est[which(est[,"kind"] == "upper"), "pKi" ]
    pE_upper <- sum( v.pE_upper * ratioXi/sum(ratioXi) )
    pK_upper <- sum( v.pK_upper * ratioXi/sum(ratioXi) )
    
    parms_lower <-  c(beta = beta,  delta = delta, pE = pE_lower, pK = pK_lower, c = c)
    parms_est <-  c(beta = beta,  delta = delta, pE = pE_estimate, pK = pK_estimate, c = c)
    parms_upp <-  c(beta = beta,  delta = delta, pE = pE_upper, pK = pK_upper, c = c)	
    
    
    colline <- c(rgb(0.8,0.5,0), rgb(0.5,0.5,0.5), rgb(0,0.8,0.8), rgb(0.8,0,0.8))
    colCI <- c(rgb(0.8,0.5,0,0.2), rgb(0.5,0.5,0.5,0.2), rgb(0,0.8,0.8,0.2), rgb(0.8,0,0.8,0.2))
    
    
    pdf(fn.pdf, width=4, height=4)
    
      plot(c(0,Npass), c(0,100), xlab="passage", ylab="% PB2-627K", axes=F, type="n", main = animal )
      axis(1,lwd=2)
      axis(2,lwd=2)
    
      i<-0			
      for(startVK in v.startpercK){
      
        i <- i+1
      
        cinit <- c(U = U, IE=IE, IK=IK, VE=(1-startVK/100)*400, VK=startVK/100*400 )
      
        est_lower <- predicting_passage_model2(cinit, parms_lower, t, Npass, Tpass)
        est_est <- predicting_passage_model2(cinit, parms_est, t, Npass, Tpass)
        est_upp <- predicting_passage_model2(cinit, parms_upp, t, Npass, Tpass)
        
        polygon(c(0:Npass, Npass:0), c(est_lower[, "percentPB2.627K"], rev(est_upp[, "percentPB2.627K"]) ), col= colCI[i], border=NA )
        lines(0:Npass, est_est[,"percentPB2.627K"], lwd=2, col = colline[i] )
      
    }
    dev.off()
  }
}



display_risk_scores_CIwithbootstrap <- function(fn.estCI, fn.ANP, v.animals, startpercK, fn.out, Npass = 5, U = 4e5, IE=0, IK=0, sumVEVK=400, t = seq(0,5,0.1), Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  # function to plot risk scores
  
  ### input
  # fn.estCI	    file name with best estimates from passage experiments and CI fn.estCI <- "data/results/estimation_virus_production/bestfits_pEi_pKi_withbootstrapCI.csv"
  # fn.ANP		    file name of ANP ratios in species, fn.ANP <- "data/input/ANPratios.csv"
  # v.animals 	  vector with animal names from ANPratios.csv in the order they should appear on the x-axis
  # startpercK	  starting percentage of K-virus 
  # fn.out		    file name of output .pdf file
  # remaining variables: model parameters
  ### output
  # one .pdf with the risk scores with 95%CI

  est <- read.csv(fn.estCI)
  anp <- read.csv(fn.ANP)
  
  nanim <- length(v.animals)
  
  pdf(fn.out,width=8, height=8,useDingbats=FALSE)
  
    plot(c(0.5,nanim+0.5), c(-1,1), xlab="", ylab="risk for human adaptation",axes=F, type="n")
    
    polygon(c(0.5,0.5,nanim+0.5,nanim+0.5), c(0,1,1,0) , col= rgb(1,0,0,0.1), border=NA)
    polygon(c(0.5,0.5,nanim+0.5,nanim+0.5), c(0,-1, -1,0), col= rgb(0,0,1,0.1), border=NA)
    
    axis(2,lwd=2)
    axis(1,lwd=2, at=1:nanim, labels=v.animals , las=2)
    
  
    i <- 0
    for(animal in v.animals){
    
      i <- i+1
    
      ratioXi <- as.numeric(anp[animal,])
      
      v.pE_lower <- est[which(est[,"kind"] == "lower"), "pEi" ]
      v.pK_lower <- est[which(est[,"kind"] == "lower"), "pKi" ]
      pE_lower <- sum( v.pE_lower * ratioXi/sum(ratioXi) )
      pK_lower <- sum( v.pK_lower * ratioXi/sum(ratioXi) )
      
      v.pE_estimate <- est[which(est[,"kind"] == "estimate"), "pEi" ]
      v.pK_estimate <- est[which(est[,"kind"] == "estimate"), "pKi" ]
      pE_estimate <- sum( v.pE_estimate * ratioXi/sum(ratioXi) )
      pK_estimate <- sum( v.pK_estimate * ratioXi/sum(ratioXi) )			
        
      v.pE_upper <- est[which(est[,"kind"] == "upper"), "pEi" ]
      v.pK_upper <- est[which(est[,"kind"] == "upper"), "pKi" ]
      pE_upper <- sum( v.pE_upper * ratioXi/sum(ratioXi) )
      pK_upper <- sum( v.pK_upper * ratioXi/sum(ratioXi) )
    
      parms_lower <-  c(beta = beta,  delta = delta, pE = pE_lower, pK = pK_lower, c = c)
      parms_est <-  c(beta = beta,  delta = delta, pE = pE_estimate, pK = pK_estimate, c = c)
      parms_upp <-  c(beta = beta,  delta = delta, pE = pE_upper, pK = pK_upper, c = c)	
      
      cinit <- c(U = U, IE=IE, IK=IK, VE=(1-startpercK/100)*400, VK=startpercK/100*400 )
      
      est_lower <- predicting_passage_model2(cinit, parms_lower, t, Npass, Tpass)
      est_est <- predicting_passage_model2(cinit, parms_est, t, Npass, Tpass)
      est_upp <- predicting_passage_model2(cinit, parms_upp, t, Npass, Tpass)
      
      #calculate risk scores
      risk_lower <- normalised_risk_score(est_lower)
      risk_est <- normalised_risk_score(est_est)
      risk_upp <- normalised_risk_score(est_upp)
      
      lines(rep(i,2),c(risk_lower,risk_upp), lwd=1)
      lines(i + c(-.1,.1), rep(risk_lower,2), lwd=1)
      lines(i + c(-.1,.1), rep(risk_upp,2), lwd=1)
      points(i, risk_est, pch=16, cex = 1.75)
  }
  dev.off()
}

####################################################
######       Example for function calls       ######
####################################################

# # Note this section is commented out such that this file can be sourced without producing any output.
# 
# # source("passage_predictions/passage_predictions_risk_scores.R")
#
# fol.out <- "passage_predictions/test/"
# system(paste0("mkdir ", fol.out))
# 
# ### Plot passaging predictions; these are Figure 5 and SupplFig 2
# fn.estCI <- "data/results/estimation_virus_production/bestfits_pEi_pKi_withbootstrapCI.csv"
# fn.ANP <- "data/input/ANPratios.csv"
# fn.out <- paste0(fol.out, "prediction_")
# v.startpercK <- c(1,80,20,95)
# host_species_passaging_predictions_CIwithbootstrap(fn.estCI, fn.ANP, fn.out, v.startpercK)
#
# ### Plot risk scores with confidence intervals; the graph for startpercK <- 1 is Figure 6a
# fn.estCI <- "data/results/estimation_virus_production/bestfits_pEi_pKi_withbootstrapCI.csv"
# fn.ANP <- "data/input/ANPratios.csv"
# v.animals <- c("Mammals", "Swallow", "Magpie", "Blackbird", "Goose", "Swan", "Turkey", "Quail", "Chicken",  "Duck","Gull")
# for( startpercK in c(1,80,20,95)){
#   fn.out<- paste0(fol.out, "risk_scores_withbootstrapCI_startpercK", startpercK,".pdf")
#   display_risk_scores_CIwithbootstrap(fn.estCI, fn.ANP, v.animals, startpercK, fn.out)
# }


