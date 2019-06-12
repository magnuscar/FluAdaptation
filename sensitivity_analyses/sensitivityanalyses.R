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



pEpK_passagingconc_riskscore <- function(fn.data, fn.ANPratios, fac=2, fn.pEpK, plotXifits=T, v.startpercK, fol.out){
  
  # fn.data			file name of passaging data with Xi only cells fn.data <- "~/work/projects/Hale_InfluenzaAdaptation/data/2018-08-16_ANP_passage_experiment_X2_no_rep2.csv"
  # fn.ANPratios		file name of ANP ratios, fn.ANPratios <- "~/work/projects/Hale_InfluenzaAdaptation/data/ANPratios.csv"
  # fac				factor of how much the estimation of Xi only data should deviate form the minimal value
  # fn.pEpK			file name of estimates of pE1-3 and pK1-3, fn.pEpK <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/pEpKestimates_varying_starting_conditions.csv"
  # v.startpercK			vector with starting values for K strain in risk analysis and passaging prediction plots
  # fol.out			folder name where the different analyses parts should be stored, eg fol.out <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/"
  
  ###########################
  # definition of initial conditions 
  # initial conditions: 0.1639 is the average fraction of K at the start of the passaging experiment
  cinit<- c(U = 4e5, IE=0, IK=0, VE=(1-0.1639)*400, VK=0.1639*400 )
  beta <- 2.7e-6
  delta <- 4  
  c <- 3
  t <- seq(0,5,0.1)
  Npass <- 5
  Tpass <- 6
  
  ###########################
  
  
  
  
  ordata <- read.csv(fn.data)
  
  ANPratios <- read.csv(fn.ANPratios)
  estpek <- read.csv(fn.pEpK)
  
  # reduction of all outcomes to those that have small deviance from the minimal RSS-value
  # minimal RSS values: X1 9.97, X2 314.09, X3 22.46
  estpek <- estpek[ which(estpek[,"valX1"]< fac*min(estpek[,"valX1"]) & estpek[,"valX2"]< fac*min(estpek[,"valX2"]) & estpek[,"valX3"]<fac*min(estpek[,"valX3"])) , ]
  
  
  # print ranges of combined pE/pK for the different animals
  for(animal in dimnames(ANPratios)[[1]]){
    
    pEnew <- ANPratios[animal,"X1"]/100 * estpek[,"pE1"] + ANPratios[animal,"X2"]/100 * estpek[,"pE2"] + ANPratios[animal,"X3"]/100 * estpek[,"pE3"]
    pKnew <- ANPratios[animal,"X1"]/100 * estpek[,"pK1"] + ANPratios[animal,"X2"]/100 * estpek[,"pK2"] + ANPratios[animal,"X3"]/100 * estpek[,"pK3"]
    
    print(animal)
    print(range(pEnew/pKnew))
    
  }
  
  
  # plot the fits of the chosen parameter estimates to the original passaging data 
  fn.passest <- paste0(fol.out,"estimates_passaging.pdf")
  if(plotXifits){
    pdf(fn.passest, width=12, height=4)
    
    
    par(mfrow=c(1,3))
    
    for(i in 1:3){
      
      xi <- paste0("X",i)
      pEi <- paste0("pE",i)
      pKi <- paste0("pK",i)
      plot(c(0,Npass), c(0,100), xlab="passage", ylab="% PB2-627K", axes=F, main=xi, type="n" )
      axis(1,lwd=2)
      axis(2,lwd=2)
      for(j in 1:(dim(estpek)[1])){
        parms <-  c(beta=beta,  delta=delta, pE = estpek[j,pEi], pK= estpek[j,pKi], c=c)				
        zwpredictions <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
        
        lines(zwpredictions[,"passage"], zwpredictions[,"percentPB2.627K"],lwd=0.5,col="lightgray" )
        
      }
      zw <- ordata[which(ordata[,"ANP"]==xi),]
      points(zw[,"passage"], zw[,"percentK"], pch=16)	
      
      
    }
    dev.off()
  }
  
  
  # table with animal passaging predictions and riskscore
  pprs <- matrix(NA, nrow=0, ncol= 20, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" )) )
  
  for(animal in dimnames(ANPratios)[[1]]){
    for(startpercK in v.startpercK ){
      for(i in 1:(dim(estpek)[1])){
        zw.pprs <- matrix(c(animal,startpercK, ANPratios[animal,], estpek[i, c(paste0("pE", 1:3), paste0("pK", 1:3)) ],rep(NA,9)), nrow=1, ncol= 20, dimnames = list(c(), c("animal","startpercK", paste0("X", 1:3), paste0("pE", 1:3), paste0("pK", 1:3 ), "pE", "pK", "ratiopEpK", paste0("percKpass", 1:Npass), "riskscore" )) )
        
        zw.pE <- sum(ANPratios[animal,]/100*estpek[i, c(paste0("pE", 1:3))])
        zw.pK <- sum(ANPratios[animal,]/100*estpek[i, c(paste0("pK", 1:3))])
        zw.pprs[1, "pE"] <- zw.pE
        zw.pprs[1, "pK"] <- zw.pK
        zw.pprs[1, "ratiopEpK"] <- zw.pE/zw.pK
        
        cinit["VE"] <- (1 - startpercK/100)*400
        cinit["VK"] <- (startpercK/100)*400				
        parms <-  c(beta=beta,  delta=delta, pE = zw.pE, pK= zw.pK, c=c)
        zwpredictions <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
        zw.pprs[1, paste0("percKpass", 1:Npass)] <- zwpredictions[1 + 1:Npass,"percentPB2.627K"]
        
        zw.pprs[1,"riskscore"] <- normalised_risk_score(zwpredictions) 
        
        pprs <- rbind(pprs, zw.pprs)
      }
      
    }
    
  }	
  
  write.table(pprs, paste0(fol.out,"estimation_of_passaging_diff_animals.csv"),quote=F, append=F, col.names=T, row.names=F, sep =", ")
  
  
  # passaging prediction graphs
  fn.passpred <- paste0(fol.out,"estimation_of_passaging_diff_animals.pdf")
  
  # color scheme for passaging predictions
  v.col <- c("lightgray","darkgray","lightblue","darkblue")
  
  pdf(fn.passpred, width=2*4, height=4*ceiling(dim(ANPratios)[1]/2))
  
  par(mfrow = c(ceiling(dim(ANPratios)[1] / 2),2) )
  for(animal in dimnames(ANPratios)[[1]]){
    
    zw <- pprs[which(pprs[,"animal"]==animal),]
    
    plot(c(0,Npass), c(0,100), xlab="passage", ylab="% PB2-627K", axes=F, main=animal, type="n" )
    axis(1,lwd=2)
    axis(2,lwd=2)
    
    icol <- 0
    for(startpercK in v.startpercK){
      icol <- icol+1
      zwpercK <- zw[which(zw[, "startpercK"] == startpercK ),]
      for(j in 1:(dim(zwpercK)[1])){
        lines(0:Npass, zwpercK[j, c("startpercK", paste0("percKpass",1:Npass))], lwd=0.5, col=v.col[icol])
        
      }				
    }
    
  }
  
  dev.off()
  
  
}



reduction_passaging_estimates_blackbirdselE <- function(fn.pred, fn.out){
  
  # function to reduce the passaging predictions to only those parameters estimates for which blackbird selects for E
  
  ### input
  # fn.pred		file name of predictions for all parameter estimates "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/estimation_of_passaging_diff_animals.csv"
  # fn.out		file name of predicitons subset fn.out<-"~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/estimation_of_passaging_diff_animals_only_blackbirdselectionE.csv"
  ### output
  
  
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




violin_plot_riskscore <- function(fn.data, startpercK, horizontal = F, fn.pdf){
  
  # fn.data 		file name of passaging data for all animals fn.data <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/estimation_of_passaging_diff_animals.csv"
  # startpercK		starting value for K strain in risk analysis 
  # horizontal		whether violin plot should be drawn horizontal (if yes =T)
  # fn.pdf		file name of pdf file
  
  dat <- read.csv(fn.data)
  dat <- dat[which(dat[,"startpercK"]==startpercK), ]
  
  pdf(fn.pdf,width=10, height=10)
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
            horizontal=horizontal,
            col="lightgray",
            drawRect =F,
            ylim=c(-1,1),
            names = c("Mammals", "Swallow", "Magpie", "Blackbird", "Goose", "Swan", "Turkey", "Quail", "Chicken",  "Duck","Gull")#,
            #las = 2
    )
  
    xmax <- 12				
    xred <- c(0,xmax,xmax,0)
    yred <- c(0,0,1,1)
    xblue <-c(0,xmax,xmax,0)
    yblue <- c(0,0,-1,-1)
    
  
    if(horizontal){
      ymax <-13
      yred <- c(0,ymax,ymax,0)
      xred <- c(0,0,1,1)
      yblue <-c(0,ymax,ymax,0)
      xblue <- c(0,0,-1,-1)
    }
    polygon(xred, yred , col= rgb(1,0,0,0.1), border=NA)
    polygon(xblue, yblue, col= rgb(0,0,1,0.1), border=NA)
  dev.off()
}




###################################################
######       Example for function call       ######
###################################################

## Note: this section is commented out such that this file can be sourced without producing any output.

# # source("sensitivity_analyses/sensitivityanalyses.R")
#
# fol.out <- "sensitivity_analyses/test/"
# system(paste0("mkdir ", fol.out))
# # Sensitivity analysis 1: testing different assumptions on the range of the virus production rates pEi and pKi (SuppFigure 3a)
# 
# # (i) generate a list of pEi estimates for fixed pKi (which ranges between 1 and 200) WARNING: this takes about 10min to run
# 
start.pE <-12
v.pK <- 1:1
lower <- 0
upper <- 200
fn.data <- "data/input/ANP_passage_experiment.csv"
fn.out <- paste0(fol.out, "sensana1_estimatespEpK_pKconst.csv")
# for the remaining parameters we used the pre-defined ones
list.with.PEonly.estimates(start.pE, v.pK, lower, upper, fn.data, fn.out)
# # (ii) generate a list of pEi and pKi estimates WARNING: this takes about 1h to run
v.pEmin <- 0
v.pKmin <- 0
v.pEmax <- seq(5, 200, by=5)
v.pKmax <- seq(5, 200, by=5)
start.par <- c(4,4)
fn.data <- "data/input/ANP_passage_experiment.csv"
fn.out <- paste0(fol.out, "sensana1_estimatespEpK_pKvar.csv")
list.with.PEK.estimates(v.pEmin, v.pKmin, v.pEmax, v.pKmax, start.par, fn.data, fn.out)


# 
# 
# 
# 
# 
# 
# 
# 
# #### 11.10.2018
# #### creating the plots for the "send-out version" of the manuscript
# 
# # riskscore for startpercK =95% and 99%
# # outsourced on cluster
# '~/work/projects/Hale_InfluenzaAdaptation/scripts/clusterjob_passpred.R'
# system("bsub -W 24:00 R CMD BATCH /cluster/home/magnusc/haleflu/clusterjob_passpred.R")
# 
# 
# 
# #### this plot is independent of the starting concentration, as we only look at pE/pK ratio
# source('~/work/projects/Hale_InfluenzaAdaptation/scripts/flu_virusdyn_ANP32A.R')
# fn.data <- "~/work/projects/Hale_InfluenzaAdaptation/data/2018-08-16_ANP_passage_experiment.csv"
# fn.pred <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/estimation_of_passaging_diff_animals.csv"
# fol.out<-"~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/"
# # passaging experiments
# pEpK_estimationplots_blackbirdselE(fn.data, fn.pred,  fol.out)
# 
# 
# 
# # reducing predicitons to thoses parameter estimates for which blackbird selects for E
# source('~/work/projects/Hale_InfluenzaAdaptation/scripts/flu_virusdyn_ANP32A.R')
# fn.pred <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/estimation_of_passaging_diff_animals.csv"
# fn.out <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/estimation_of_passaging_diff_animals_only_blackbirdselE.csv"
# reduction_passaging_estimates_blackbirdselE(fn.pred, fn.out)
# 
# 
# 
# ### for the following we use the passaging predictions with 1,20,80,95 start conc
# source('~/work/projects/Hale_InfluenzaAdaptation/scripts/flu_virusdyn_ANP32A.R')
# fn.pred <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/estimation_of_passaging_diff_animals_only_blackbirdselE.csv"
# fol.out <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/passest_1208095/"
# # prediction of passaging experiments for different host species
# v.startpercK <- c(1,80,20,95)
# host_species_passaging_predictions_blackbirdselE(fn.pred, v.startpercK, fol.out)
# 
# fol.out <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/passest_12080/"
# # prediction of passaging experiments for different host species
# v.startpercK <- c(1,80,20)
# host_species_passaging_predictions_blackbirdselE(fn.pred, v.startpercK, fol.out)
# 
# 
# # violin plots
# fn.data <- "~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/estimation_of_passaging_diff_animals_only_blackbirdselE.csv"
# 
# 
# for(startpercK in c(1,80,20,95, 99)){
#   fn.pdf <- paste0("~/work/projects/Hale_InfluenzaAdaptation/results/uncertainty_analysis/allXidata/only_blackbirdselectionE/violinplots_riskscore_startpercK",startpercK,".pdf") 
#   violin_plot_riskscore(fn.data, startpercK, horizontal = F, fn.pdf)
# }
# 
# 
# 
# # Sensitivity analysis 2: testing different assumptions on the model parameters beta, delta, c; here, they stay the same for the different virus variants (SuppFigure 3b)
# # 1) when estimating all parameters
# ### sensitivity analysis for different beta, delta and c values for the two different virus strains
# # see script ~/Dropbox/projects_atIMV/Hale_InfluenzaAdaptation/scripts/cluster_2019/sens_ana_cluster.R
# 
# # 1) post processing of diff beta, delta c : sent to cluster scripts/cluster_2019/job_postprocess_diffbetadeltac.R
# source("scripts/flu_virusdyn_ANP32A.R")
# #setwd("/Users/magnus.carsten/Dropbox/projects_atIMV/Hale_InfluenzaAdaptation")
# fn.est <- "results/revisions/sensitivityana_diffpK/estimates_smallrangeforbetadeltac/estimatates_varyingbetadeltac_globalpEpK_allANP_4.csv"
# facrss <- 1.1
# pKone <- F
# fn.ANPexp <- "data/2018-08-16_ANP_passage_experiment.csv"
# fn.ANPratios <- "data/ANPratios.csv"
# fol.out <-"results/revisions/sensitivityana_diffpK/estimates_smallrangeforbetadeltac/"
# 
# postprocessing_sensitivityAna(fn.est, facrss, pKone, fn.ANPexp, fn.ANPratios, fol.out)
# 
# # Sensitivity analysis 3: testing different assumptions on the model parameters beta, delta, c; here, they vary for the different virus variants (SuppFigure 3c)
# # 3) post-processing: the estimates for different beta, delta, c values for the two different viral strains: on cluster. See job_postprocess_diffbetadeltac.R
# ### 04.04.2019
# # As there is no violinplot - package installed on the cluster, this part could not have been performed. Thus, this analysis will run here now:
# setwd("/Users/magnus.carsten/Dropbox/projects_atIMV/Hale_InfluenzaAdaptation")
# source("scripts/flu_virusdyn_ANP32A.R")
# fn.data <- "results/revisions/sensitivityana_diffbetadeltac/estimates_smallrange1.1/estimation_of_passaging_diff_animals.csv"
# fol.out <- "results/revisions/sensitivityana_diffbetadeltac/estimates_smallrange1.1/"
# postprocessing_sensitivityAna_diffbetadeltac_violinplotsonly( fn.data, fol.out)
# 
