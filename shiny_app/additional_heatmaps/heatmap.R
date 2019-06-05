##########################################################
##########         GENERATING HEATMAPS          ##########
##########################################################

###### Required libraries

library(grid)
library(graphics)


######  Functions

### Calculate risk scores

risk_score_heatmap <- function(v.pE, v.pK, squarenumber = 21, betaK=2.7e-6, deltaK=4, cK=3, betaE=2.7e-6, deltaE=4, cE=3, cinit = c(U = 4e5, IE=0, IK=0, VE=(1-0.2)*400, VK=0.2*400 ), t= seq(0,5,0.1), Npass=5, Tpass=6, fn.out){
  
  # This function calculates a matrix of riskscores and stores the matrix
  
  ### input:
  # v.pK		    	vector with values for pKi
  # v.pE			    vector with values for pEi
  # squarenumber	number of points on x (fX1) and y (fX2) axis
  # betaK         infection rate with K-virus
  # betaE         infection rate with E-virus
  # deltaK        infected cell death rate of cells infected with K-virus
  # deltaE        infected cell death rate of cells infected with E-virus
  # cK            death rate of K-virus
  # cE            death rate of E-virus
  # cinit         vector of initial conditions
  # t             vector with time steps
  # Npass		      number of passages
  # Tpass		      time (in steps) in simulation that corresponds to a passaging round
  # fn.out	      file name of output file .csv
  ### output:
  # matrix with fX1, fX3 frequencies and riskscore; this matrix is also stored under fn.out
  
  risk.matrix <- matrix(NA, nrow = 0, ncol = 3, dimnames=list( c(), c("fX1", "fX3","riskscore")  ))
  write.table(risk.matrix, fn.out, quote = F, sep = ",", row.names=F)

  for(fX1 in seq(0,1,length.out = squarenumber)){
    for(fX3 in seq(0,1,length.out = squarenumber)){
      
      newrow <- matrix(c(fX1, fX3, NA), nrow = 1, ncol = 3, dimnames=list( "", c("fX1", "fX3","riskscore")  ))
      
      if(fX1 + fX3 <= 1){
        
        pE <- sum( c(fX1, 1 - fX1 - fX3, fX3) * v.pE )
        pK <- sum( c(fX1, 1 - fX1 - fX3, fX3) * v.pK )
        
        if( betaK==betaE & deltaK==deltaE & cK==cE){
          parms <- c(beta = betaK, delta = deltaK, pE = pE, pK = pK, c = cK)
          zwischen <- predicting_passage_model2(cinit, parms, t, Npass, Tpass)
        }
        else{
          parms <- c(beta_E = betaE, beta_K = betaK, delta_E = deltaE, delta_K = deltaK , c_E = cE, c_K = cK, pE = pE, pK = pK)
          zwischen <- predicting_passage_model2_diffratesforvir(cinit, parms, t, Npass, Tpass)
        } 
        newrow[1,"riskscore"] <- normalised_risk_score( zwischen ) 
      }
      
      write.table(newrow, fn.out, append =T, quote = F, sep = ",", row.names=F, col.names=F)
      risk.matrix <- rbind(risk.matrix, newrow)
    }
  }
  return(risk.matrix)
}



###### Help-functions to display the heatmap

### Calculate for which fX3 value the sign of the riskscore changes

return.change.riskscore <- function(fX1, data ){
  
  # this function finds the value of fX3 where the sign of the riskscore changes for a given value of fX1
  
  ### input:
  # fX1     fraction of X1
  # data    a matrix of riskscores as created with function risk_score_heatmap
  ### output:
  # NA if the sign of the risk score does not change, otherwise the value of fX3 where this change happens
  
  dat <- data[which(data[, "fX1"] == fX1),]
  nex <- dim(dat)[1]
  xout <- NA
  
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



### Vectorization of this function in fX1

return.change.riskscore.v <- function(FX1,data){ sapply(FX1, return.change.riskscore, data=data) }



### Output the values of fX1 and fX3 as a list

get.0lines <- function(data){
  
  ### input:
  # data        a matrix of riskscores as created with function risk_score_heatmap
  ### output:
  # list with fX1 and fX3 values 
  
  dat <- data[complete.cases(data), c("fX1", "fX3", "riskscore") ]
  FX1 <- unique(dat[,"fX1"])	
  FX3 <- return.change.riskscore.v(FX1,dat)
  firstnonNA <- which(!is.na(FX3) )[1]
  
  if(firstnonNA > 1){
    FX3[1:(firstnonNA -1)]<- -0.2
  }
  return(list(fX1= FX1, fX3 = FX3))
}



######### threshold line that marks the jumps from negative to positive in steps

threshold_line_stepfunction <- function(riskmatrix){
  
  # function that returns the x and values for a threshold line
  
  ### input: 
  # riskmatrix    a matrix of riskscores as created with function risk_score_heatmap
  ### output: 
  # list with fX1 and fX3 values
  
  
  rm <- riskmatrix[!is.na(riskmatrix[,"riskscore"]), ]
  
  # initialize v_x and v_y
  v_x <- c()
  v_y <- c()
  
  FX1 <- unique(rm[,"fX1"])
  
  for(i_x in 1:length(FX1)){
    
    FX3 <- unique(rm[which(rm[,"fX1"] == FX1[i_x]), "fX3"] )
    
    risk <- -1
    i_y <- 0
    while(risk < 0){
      i_y <- i_y + 1
      risk <- rm[which(rm[,"fX1"] == FX1[i_x] & rm[,"fX3"] == FX3[i_y]  ), "riskscore"]			
      
      if(i_y == length(FX3)){
        risk <- 1
      }
    }
    
    if(i_x != length(FX1)){
      v_x <- c(v_x, FX1[i_x + 0:1])
      v_y <- c(v_y, rep(FX3[i_y],2))
    }
    else{
      v_x <- c(v_x, FX1[i_x])
      v_y <- c(v_y, FX3[i_y])
    }		
  }
  return(list(fX1 = v_x - 0.025,  fX3 = v_y - 0.025))
}



### costumized color scheme

col.fred <- function(n, k) c(rev(rgb(0,0,1,alpha=seq(0,1,length=n/2)^k)),"white",rgb(1,0,0,alpha=seq(0,1,length=n/2)^(k)))



### Create a heatmap

heatmap_withseparatingline <- function(fol.dat,  betaK, deltaK, cK, betaE, deltaE, cE, startpercK, fn.ANP, fn.out , thresholdline = "average"){
  
  # this function reads in a riskscore matrix from the folder fol.dat and produces a pdf of a heatmap
  
  ### input:
  # fol.dat       folder where the riskscore matrices are stored, e.g. fol.dat <- "data/heatmaprawdata/"
  # startpercK		1, 20, 40, 80
  # betaK, betaE		8.8e-7, 2.7e-6, 1e-5 
  # deltaK, deltaE	2.6, 4, 6.1
  # cK, cE			2.4, 3, 3.6
  # fn.ANP		file name of ANP ratios in species, fn.ANP <- "data/ANPratios.csv"
  # thresholdline		if "average", it averages between positive and negative jumps in risks
  #					if "step", it marks the jumps
  ### output:
  # pdf with the heatmap 
  
  
  anp <- read.csv(fn.ANP)
  
  risk.matrix <- read.csv(paste0(fol.dat, "riskmatrix_betaK",betaK,"_deltaK", deltaK,"_cK", cK,"_betaE", betaE,"_deltaE",deltaE,"_cE", cE, "_percK", startpercK,".csv"))
  
  if(thresholdline == "average"){
    threshold <- get.0lines(risk.matrix)
  }
  if(thresholdline == "step"){
    threshold <- threshold_line_stepfunction(risk.matrix) 
  }
  
  plot.int <- levelplot(riskscore ~ fX1*fX3, risk.matrix, at = seq(-1,1,0.02), col.regions = col.fred(100,1.1),
                        panel = function(...){
                          panel.levelplot(...)
                          grid.points(anp$X1 / 100, anp$X3 / 100, pch=16)
                          panel.lines(threshold$fX1, threshold$fX3, lwd=1, col= "gray")
                          grid.text(dimnames(anp)[[1]], x=anp[,"X1"]/100, y=anp[,"X3"]/100, just = c("right", "centre"))
                        }
  )
  pdf(fn.out, width=8, height=8, useDingbats=F)
  plot.int
  dev.off()
  
}





