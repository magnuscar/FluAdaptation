#######################################################################
##########     Testing host adaptation in sequence data      ##########
#######################################################################

###### Required libraries

library(DescTools)
library(graphics)

###### set working directory
# setwd("/path/to/FluAdaptation/") ### <<<--- you have to add the path to the folder on your machine here and uncomment this line!!!

###### source required R scripts
source("ANP_model.R")




########################################################################
######                  Definition of functions                   ######
########################################################################


prop.test.table <- function(fn.data, test_animal = "Gull", fn.out){
  
  # function to perform a proportion test on the surveillance data
  
  ### input
  # fn.data			  file name of input data fn.data <- "results/revisions/hostadaptation/hostadaptation.csv"
  # test_animal		the animal with which the others should be compared
  ### output
  # matrix with pvalues of the proportion test with abundance data of test_animal as baseline
  
  dat <- read.csv(fn.data)
  v.animals <- dimnames(dat)[[1]]
  
  # extract the mutations
  v.mut <- unlist(strsplit( dimnames(dat)[[2]] , "_")) 
  v.mut <- unique( v.mut[!(v.mut %in% c("numb", "total"))] )

  out <- matrix(NA, nrow = dim(dat)[1], ncol = length(v.mut), dimnames = list(v.animals, paste0("pvalue_", v.mut)))
  
  for(animal in v.animals[v.animals != test_animal ]){
     for(mut in v.mut){
      
      zw <- prop.test(c(dat[test_animal, paste0("numb_", mut)], dat[animal, paste0("numb_", mut)]), c(dat[test_animal, paste0("total_", mut)], dat[animal, paste0("total_", mut)]), alternative = "less")
      out[animal, paste0("pvalue_",mut)] <- signif(zw$p.value,4)
      
    }
  }
  
  write.table(out, fn.out, quote=F, sep = ",", row.names=T, col.names=T)
}





display_pMammal_CI_vers3 <- function(fn.data, v.mut = c("K", "N", "Q", "D", "I") , ref_animal = "Gull", test_animal="Gull",fn.pvalues, fol.out){
  
  # function to display the estimates of a binomial estimates of the population wide abundance of the observed mutations with CI and significance in comparison of a reference animal
  
  ### input:
  # fn.data			  file name of input data fn.data <- "data/input/hostadaptation/hostadaptation.csv"
  # v.mut				  mutations to look at
  # ref_animal		reference animal for blue/red coloring and proportion tests of significance
  # test_animal		animal for testing and determination of the pvalue from the following table
  # fn.pvalues	 	file name of where pvalues from proportion test are stored
  # fn.out			  file name of output
  ### output:
  # graphical outputs, i.e. one .pdf file with an overview of the two statistical analyses for each mutation in v.mut
  
  data <- read.csv(fn.data)
  
  v.animals <- dimnames(data)[[1]]
  
  iref <- which(v.animals == test_animal)
  
  # load pvalues
  pv <- read.csv(fn.pvalues)
  pv[is.na(pv)] <- 10
  
  for(mut in v.mut){
    fn.out <- paste0(fol.out,mut,".pdf")
    pdf(fn.out, height= 4, width = 4, useDingbats=FALSE)
    
    num_mut <- paste0("numb_", mut)
    total_mut <- paste0("total_", mut)		
    #lower.p <- BinomCI(data[ref_animal, num_mut], data[ref_animal, total_mut], method = "wilson", sides = "t")[1,"lwr.ci"]
    upp.p <- BinomCI(data[ref_animal, num_mut], data[ref_animal, total_mut], method = "wilson", sides = "t")[1,"upr.ci"]
    
    plot(  c(0,1), c(1, length(v.animals)), xlab="abundance of mammal adapted viruses", ylab="", axes=F, type="n", main= paste0("Mutation: ",mut), xlim = c(0,0.8))
    
      axis(1, lwd=2)
      axis(2, lwd = 2, labels = v.animals, at=1:length(v.animals), las=2, cex=0.7)
      
      polygon(c(upp.p,0.8, 0.8,upp.p), c(0.5,0.5,length(v.animals)+0.5,length(v.animals) +0.5 ), col= rgb(1,0,0,0.1), border = NA )
      polygon(c(0,upp.p,upp.p,0), c(0.5,0.5,length(v.animals)+0.5,length(v.animals) +0.5 ), col= rgb(0,0,1,0.1), border = NA )
      i <- 0
      k <- 0
      for(animal in v.animals){
        
        i<- 1+i
        
        zw <-  BinomCI(data[animal, num_mut], data[animal, total_mut], method = "wilson", sides = "t")
        
        points( zw[1,"est"], i, pch = 16) 
        
        lines(zw[1, c("lwr.ci","upr.ci")], rep(i,2), lwd=1)
        lines(rep(zw[1, c("lwr.ci")],2), i + c(-.1,.1), lwd=1)
        lines(rep(zw[1, c("upr.ci")],2), i + c(-.1,.1), lwd=1)
        
        
        ####
        # add stars for p-value
        if(pv[animal,paste0("pvalue_",mut)] < 0.05){
          k <- k+1
          
          x <- 0.62 + 0.05*k
          lines(rep(x ,2) , c(i, iref ), lwd=0.5, col = "black")
          lines(rep(x ,2) + c(0,-0.02), c(i, i), lwd=0.5, col = "black")
          lines(rep(x ,2) + c(0,-0.02), c(iref, iref ), lwd=0.5, col = "black")
          
          # one star
          if( pv[animal,paste0("pvalue_",mut)] > 0.01 & pv[animal,paste0("pvalue_",mut)] <0.05){
            points(x + 0.02, (iref- i)/2 + i, pch=8 , col="black", cex=0.5)
          }
          # two stars 
          if( pv[animal,paste0("pvalue_",mut)] > 0.001 & pv[animal,paste0("pvalue_",mut)] <0.01){
            k <- k+0.25
            points(x + c(0.02,0.04), rep((iref- i)/2 + i, 2), pch=8 , col="black", cex=0.5)
          }
          
          # three stars
          if(pv[animal,paste0("pvalue_",mut)] <0.001){
            k <- k+0.5
            points(x + c(0.02,0.04,0.06), rep((iref- i)/2 + i, 3), pch=8 , col="black", cex=0.5)
          }
        }
    }
    dev.off()
  }
}



####################################################
######       Example for function calls       ######
####################################################

# # Note this section is commented out such that this file can be sourced without producing any output.
# 
# we generate a test folder:
fol.out <- "surveillance_analysis/test/"
system(paste0("mkdir ", fol.out))
# Figure 7
# 1. Generate table with p-values from proportion test, with "Gull" as reference
# source('surveillance_analysis/host_adaptation.R')
fn.data <- "data/input/hostadaptation.csv"
fn.out <- paste0(fol.out, "pvalues_reference_Gull_hostadaptation.csv")
prop.test.table(fn.data, test_animal = "Gull", fn.out)
# 2. Graphical overview of the two statistical analyses
fn.data <- "data/input/hostadaptation.csv"
fol.out_new <- paste0(fol.out, "hostadaptation_final_mut_")
fn.pvalues <- paste0(fol.out, "pvalues_reference_Gull_hostadaptation.csv")
display_pMammal_CI_vers3(fn.data=fn.data,  fol.out=fol.out_new, fn.pvalues=fn.pvalues)