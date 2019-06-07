########################################################################
######    Displaying best estimates of virus production rates     ######
########################################################################


###### Required libraries
library(graphics)



###### set working directory
# setwd("/path/to/FluAdaptation/") ### <<<--- you have to add the path to the folder on your machine here and uncomment this line!!!



###### source required R scripts
source("ANP_model.R")



########################################################################
######                  Definition of functions                   ######
########################################################################



# function to plot passaging data alongside best fits

plot.passigingdata.withCIbootstrap <- function(fn.pass, fn.estCI, fn.pdf, cinit = c(U = 4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 ), t = seq(0,5,0.1), Npass = 5, Tpass = 6, beta = 2.7e-6, delta = 4, c = 3, average = T){
  
  # function to plot the passaging data in ANP_passage_experiment.csv alongside best estimates (see estimating_virus_production_rates.R) and CI (see CIbootstrap.Rmd)
  
  ### input:
  # fn.pass 			file name of passaging data fn.pass <- "data/2018-08-16_ANP_passage_experiment.csv"
  # fn.estCI			file name with best estimates and CI fn.estCI <- "results/revisions/CIwithbootstrap/bestfits_pEi_pKi_withbootstrapCI.csv"
  # fn.pdf			  file name of plot
  # ...           parameters as described in estimating_virus_production_rates.R
  # average			  if TRUE, the average of the passaging data is drawn as a dotted line
  ### output:
  # pdf file with data and best estimates plotted
  
  ordata <- read.csv(fn.pass)
  est <- read.csv(fn.estCI)
 
  pdf(fn.pdf, width=4, height=4, useDingbats=FALSE)
  
    plot(c(0,Npass), c(0,100), xlab="passage", ylab="% PB2-627K", axes=F, type="n" )
    axis(1,lwd=2)
    axis(2,lwd=2)
  
    for(i in 1:3){
    
      xi <- paste0("X",i)

      parms_lower <-  c(beta = beta,  delta = delta, pE = est[ which(est[,"Xi"] == i & est[,"kind"] == "lower" ), "pEi" ], pK =  est[ which(est[,"Xi"] == i & est[,"kind"] == "lower" ), "pKi" ], c = c)
      parms_est <-  c(beta = beta,  delta = delta, pE = est[ which(est[,"Xi"] == i & est[,"kind"] == "estimate" ), "pEi" ], pK =  est[ which(est[,"Xi"] == i & est[,"kind"] == "estimate" ), "pKi" ], c = c)
      parms_upp <-  c(beta = beta,  delta = delta, pE = est[ which(est[,"Xi"] == i & est[,"kind"] == "upper" ), "pEi" ], pK =  est[ which(est[,"Xi"] == i & est[,"kind"] == "upper" ), "pKi" ], c = c)						
    
      est_lower <- predicting_passage_model2(cinit, parms_lower, t, Npass, Tpass)
      est_est <- predicting_passage_model2(cinit, parms_est, t, Npass, Tpass)
      est_upp <- predicting_passage_model2(cinit, parms_upp, t, Npass, Tpass)
    
    
      col = list(rgb(1,0,0,0.1), rgb(0,1,0,0.1), rgb(0,0,1,0.1))[[i]]
      polygon(c(0:Npass, Npass:0), c(est_lower[, "percentPB2.627K"], rev(est_upp[, "percentPB2.627K"]) ), col= col, border=NA )
      lines(0:Npass, est_est[, "percentPB2.627K"] , col= c("red", "green", "blue")[i])
    
      zw <- ordata[which(ordata[,"ANP"]==xi),]
      points(zw[,"passage"], zw[,"percentK"], pch= c(15,19,17)[i], col = c("red", "green", "blue")[i])	
    
      if(average){	
        x <- c()
        y <- c()	
        for(pa in unique(zw[,"passage"])){
          x <- c(x,pa)
          y <- c(y, mean(zw[which(zw[,"passage"]==pa) ,"percentK"]))
        }		
        lines(x,y, pch= c(15,19,17)[i], col = c("red", "green", "blue")[i], lty=2)
      }
    }
  dev.off()
}	



# function to plot passaging data in mixed ANP32A expressing cells alongside the model prediction

plot.data.with.prediction.withCIbootstrap  <- function(fn.dat, fn.estCI, ratioXi, line, fn.pdf, addaverage = F, Npass = 5, U = 4e5, IE=0, IK=0, sumVEVK=400, t = seq(0,5,0.1), Tpass = 6, beta = 2.7e-6, delta = 4, c = 3){
  
  ### input:
  # fn.dat		  file name of data, eg. fn.dat <- "data/arificial_cell_lines/Line_7_9.csv"
  # fn.estCI		file name with best estimates and CI fn.estCI <- "results/revisions/CIwithbootstrap/bestfits_pEi_pKi_withbootstrapCI.csv"
  # ratioXi		  ratios of X1, X2, X3
  # addaverage	True if a line between the average points should be drawn	
  ### output:
  # a pdf file 
  
  dat <- read.csv(fn.dat)
  subdat <- dat[which(dat[,"line"] == line), ]
  
  est <- read.csv(fn.estCI)
  
  pdf(fn.pdf, width=4, height=4, useDingbats=FALSE)
  
    plot(c(0,Npass), c(0,100), xlab="passage", ylab="% PB2-627K", axes=F, type="n")
    axis(1,lwd=2)
    axis(2,lwd=2)
  
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
    
    col <- c(rgb(0,0,0,0.1),rgb(1,0,0,0.1))
    
    startVK <- mean(subdat[which(subdat[,"passage"]==0), "percentPB2.627K"])/100
    cinit <- c(U = U, IE=IE, IK=IK, VE=(1-startVK)*400, VK=startVK*400 )
    
    est_lower <- predicting_passage_model2(cinit, parms_lower, t, Npass, Tpass)
    est_est <- predicting_passage_model2(cinit, parms_est, t, Npass, Tpass)
    est_upp <- predicting_passage_model2(cinit, parms_upp, t, Npass, Tpass)
    
    polygon(c(0:Npass, Npass:0), c(est_lower[, "percentPB2.627K"], rev(est_upp[, "percentPB2.627K"]) ), col= col[1], border=NA )
    lines(0:Npass, est_est[,"percentPB2.627K"], lwd=2, col = 1 )
    points(subdat[,"passage"], subdat[,"percentPB2.627K"], pch=16, col= 1)
  
    if(addaverage){
      
      perc.av <- c()
      
      for(pass in sort(unique(subdat[,"passage"]))){
        perc.av <- c(perc.av, mean( subdat[ which(subdat[,"passage"] ==pass) ,"percentPB2.627K"] ))
      }
    
      lines(0:pass, perc.av, lwd=2, col = 1, lty = 3 )
    }
    
  dev.off()
}


####################################################
######       Example for function calls       ######
####################################################

# # Note this section is commented out such that this file can be sourced without producing any output.
# 

# # This is the script for reproducing Figure 4b of the paper:
# fol.out <- "virus_production_rates/test/"
# system(paste0("mkdir ", fol.out))
# 
# fn.pass <- "data/input/ANP_passage_experiment.csv"
# fn.estCI <- "data/results/estimation_virus_production/bestfits_pEi_pKi_withbootstrapCI.csv"
# fn.pdf <- paste0(fol.out, "bestfits_pEi_pKi_withbootstrapCI_averageline.pdf")
# plot.passigingdata.withCIbootstrap(fn.pass = fn.pass, fn.estCI = fn.estCI, fn.pdf = fn.pdf, average=T)



