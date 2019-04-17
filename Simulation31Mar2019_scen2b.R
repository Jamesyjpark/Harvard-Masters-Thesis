################################################################################
# Initial Setup
################################################################################

# Set a seed for the simulations for reproducibility
set.seed(1)

################################################################################
# Load libraries
################################################################################
library(nlme) # for lme

################################################################################
# Set parameters
################################################################################

nSim <- 2000 # Number of simulations
n <- 46 # Number of subjects
nDay <- 21 # Number of days
nObs <- 20 # Number of observations
nRep <- nDay*nObs #number of measurements per subject

meanVACs <- 20.04  # Mean VACS score from the data
varVACs <- 61.42029 # VACS score variance from the data

baseStep <- 37.91353  #37.91353 Baseline mean for steps from the data
# Variance and Covariance values obtained from random intercept model in R
varSteps <- 1323.1 # 1323.1 ~ (Residual StdDev)**2 + (Intercept StdDev)**2

# Set effect sizes for VACs with steps
betaVI <- c(seq(0,0.1,length.out=4),0.2,0.5)
alpha <- 0.05

#cut off to exclude steps
stepCut<-54

################################################################################
#Store results for simulations
################################################################################

#store results for 4 methods and number of observations and subjects in reduced sets
matR <- matrix(0,ncol=6,nrow=length(betaVI))
colnames(matR) <- c("LMall","LMMall", "LMcut","LMMcut","aveRed","nRed")

################################################################################
# Cycle through effect of VACs on steps betaV
################################################################################

for(kk in 1:length(betaVI)){
  betaV<-betaVI[kk]
  
  ################################################################################
  # Run simulation
  ################################################################################
  
  # Simulation loop outside so that different dataset generated for different simulation
  for(ss in 1:nSim){
    
    cutPrint<-50
    if(floor(ss/cutPrint)==ceiling(ss/cutPrint)){print(paste(kk,"of",length(betaVI),"Betas and",ss,"in",nSim,"sims"))}
    
    ################################################################################
    # Generate data
    ################################################################################
    #generate the vacs score from a normal distribution
    VACs <- rnorm(n,mean=meanVACs,sd=sqrt(varVACs))
    VACs[VACs < 0] <- 0
    
    # Mean value of steps as a function of VACs
    meanStep <- baseStep + VACs*betaV
    
    # Initialize the empty matrix
    # Row indicates each person (1-46) 
    # columns indicate the step counts (same number of measurements for all)
    steps <- matrix(0,nrow=n,ncol=nRep)
    
    for(jj in 1:n){
      # Simulate steps data from multivariate normal
      steps[jj,] <- rnorm(nRep, mean=meanStep[jj], sd=sqrt(varSteps))
    }
    
    # check data generation
    # c(min( steps[1,]),max( steps[1,]),mean( steps[1,]),meanStep[1]) 
    
    # Each simulation iteration generates different steps per person
    colnames(steps) <- paste("m",c(1:nRep),sep="")
    steps[steps < 0] <- 0
    
      ################################################################################
      # Add VACS-unrelated steps 
      ################################################################################
    
      # Try different values for the noisecut (0, 5, 10, 15...)
      noisecut <- 0
      
      e <- list() # Empty list to hold random values
      AddObs <- 1 # Number of observations added
      varnoise <- 1 # Change variance in different scenarios
      
      for (i in 1:AddObs){
        e[i] <- list(rnorm(n, noisecut, sqrt(varnoise)))
      }

      MeAdd <- 
       function(dd){
          # Function to add different VACS-unrelated measurement at the end of each day's last observation
          # change the number of added e[[dd]] by the changing AddObs number
          
          # 1) AddObs = 1
          cbind(steps[,(1:nObs)+nObs*(dd-1)], e[[1]])
         
          # 2) AddObs = 5
          # cbind(steps[,(1:nObs)+nObs*(dd-1)], e[[1]], e[[2]], e[[3]], e[[4]], e[[5]])
         
          # 3) AddObs = 20
          #cbind(steps[,(1:nObs)+nObs*(dd-1)], e[[1]], e[[2]], e[[3]], e[[4]], e[[5]], e[[6]],e[[7]],e[[8]],e[[9]]
          #      ,e[[10]],e[[11]],e[[12]],e[[13]],e[[14]],e[[15]],e[[16]],e[[17]],e[[18]],e[[19]],e[[20]])
       }
  
      # Bind the new matrix for each day to have a total matrix
      Sim2bStep <-
        cbind(MeAdd(1),MeAdd(2),MeAdd(3),MeAdd(4),MeAdd(5),MeAdd(6),MeAdd(7),MeAdd(8),MeAdd(9),MeAdd(10),
              MeAdd(11),MeAdd(12),MeAdd(13),MeAdd(14),MeAdd(15),MeAdd(16),MeAdd(17),MeAdd(18),MeAdd(19),MeAdd(20),MeAdd(21))
    
      Sim2bStep[Sim2bStep < 0] <- 0 # Make steps nonnegative
      # ncol(Sim2bStep) - more columns added 
      
    ################################################################################
    # linear regression: ALL steps
    ################################################################################
    
    # Empty matrix with each row indicating new person & each column indicating day (1-21)
    stepD<-matrix(0,nrow=n,ncol=nDay)
    
    for(dd in 1:nDay){
        # this function now needs to account for the added column(s)
        stepD[,(dd)]<-rowSums(Sim2bStep[, (1+(AddObs*(dd-1))):(nObs+dd*AddObs)+nObs*(dd-1)])
    }
    
    # take the average number of steps per day
    # Slightly lower rowMeans compared to different simulation scenarios
    stepDA<-rowMeans(stepD)
    
    # linear regression for average total step per day
    #index reject matrix if significant association
    if(summary(lm(stepDA~VACs))$coef[2,4]<alpha){
      matR[kk,"LMall"]<-	matR[kk,"LMall"]+1
    }
    
    ################################################################################
    # linear regression: Excluding steps less than 54
    ################################################################################
    
    # EXCLUDE steps less than cut off
    stepsM2 <- Sim2bStep
    stepsM2[Sim2bStep < stepCut] <- NA
  
    # Empty matrix
    stepDM2<-matrix(0,nrow=n,ncol=nDay)
  
    for(dd2 in 1:nDay){
      stepDM2[,dd2]<-rowSums(stepsM2[,(1+(AddObs*(dd-1))):(nObs+dd*AddObs)+nObs*(dd-1)],na.rm=T)
    }
  
    # take the average number of steps per day
    # Lower than stepDA 
    stepDAM2<-rowMeans(stepDM2)
    
    #note all NAs give a row sum of 0, exclude those subjects
    stepDAM2na<-stepDAM2[stepDAM2>0]
    VACsNA<-VACs[stepDAM2>0]
    
    #check sample size for reduced set
    # c(length(stepDAM2na),length(stepDAM2na))
    matR[kk,"nRed"]<-matR[kk,"nRed"]+length(stepDAM2na)
    
    # linear regression for average total step per day
    pVacsLMcut<-summary(lm(stepDAM2na~VACsNA))$coef[2,4]
    
    #index reject matrix if significant association
    if(pVacsLMcut<alpha){
      matR[kk,"LMcut"]<-	matR[kk,"LMcut"]+1
    }
    
    ################################################################################
    # LMM: Random Intercept: ALL steps
    ################################################################################
    #reformat the matrix for LMM
    
    # dimension must change to incorporate added measurements
    matLMM<-matrix(NA,nrow=n*ncol(Sim2bStep),ncol=3)
    colnames(matLMM)<-c("subject","Vacs","Steps")

    # subject IDs
    matLMM[,"subject"]<-rep(c(1:n),each=ncol(Sim2bStep))
    
    #vacs per subject
    matLMM[,"Vacs"]<-rep(VACs,each=ncol(Sim2bStep))
    
    # New dimension check - returns 'TRUE'
    # ncol(Sim2bStep)*nrow(Sim2bStep) == nrow(matLMM)
  
    #steps per subject
    matLMM[,"Steps"] <- as.vector(t(Sim2bStep))
    
    pVACsLMM<-summary(lme(matLMM[,"Steps"]~matLMM[,"Vacs"], random = ~1 | matLMM[,"subject"],method = "REML",control = lmeControl(opt = "optim")))$tTable[2,"p-value"]
    
    #index reject matrix if significant association
    if(pVACsLMM<alpha){
      matR[kk,"LMMall"]<-	matR[kk,"LMMall"]+1
    }
    
    ################################################################################
    # LMM: Random Intercept: EXCLUDE steps below cutoff
    ################################################################################
    #reformat the matrix for LMM
    matLMMc<-matrix(NA,nrow=n*ncol(Sim2bStep),ncol=3)
    colnames(matLMMc)<-c("subject","Vacs","Steps")

    # subject IDs
    matLMMc[,"subject"]<-rep(c(1:n),each=ncol(Sim2bStep))
    
    #vacs per subject
    matLMMc[,"Vacs"]<-rep(VACs,each=ncol(Sim2bStep))
    
    #steps per subject excluded steps
    matLMMc[,"Steps"]<-as.vector(t(stepsM2))
    
    #matrix without NA - Smaller than unreduced matrix
    matLMMcNA<-na.omit(matLMMc) 
  
    # check dimension of reduced and not reduced matrices
    matR[kk,"aveRed"]<-matR[kk,"aveRed"]+dim(matLMMcNA)[1]
    
    pVACsLMMcut2<-summary(lme(matLMMcNA[,"Steps"]~matLMMcNA[,"Vacs"], random = ~1 | matLMMcNA[,"subject"],method = "REML",control = lmeControl(opt = "optim")))$tTable[2,"p-value"]
    
    #index reject matrix if significant association
    if( pVACsLMMcut2<alpha){
      matR[kk,"LMMcut"]<-	matR[kk,"LMMcut"]+1
    }    
    
    ################################################################################
    # End simulations loop
    ################################################################################
  } #end of sim loop
} #end of beta loop

################################################################################
# Compile the results
################################################################################

# Reject matrix
matR <- matR/nSim

# CHECK results
matR

# Save results for scenario 2b 
write.table(matR,file=paste("SCENnObs",nObs,"baseStep",floor(baseStep),"MeanNoiseCut",noisecut,"VarNoise",varnoise,"NumNoisePerDay",AddObs,".txt",sep=""),quote=F,row.names=F)

################################################################################
# Create and Save plot
################################################################################ 

pdf(file=paste("PLOTnObs",nObs,"baseStep",floor(baseStep),"MeanNoiseCut",noisecut,"VarNoise",varnoise,"NumNoisePerDay",AddObs,".pdf",sep=""))

# Iterate through each method to draw a line for each regrsssion model
nn <- 4

plot(-1, 1, xlim = c(min(betaVI), max(betaVI)), ylim=c(0,1), xlab = "BetaV", ylab = "", main = " ")

for (pp in 1:nn){
  lines(betaVI, matR[,pp], pch = pp, col = pp, type = "b")
}

legend("topleft", c("LM: all steps","LMM: all steps", "LM: cutoff for steps","LMM: cutoff for steps"), col = c(1:nn), pch = c(1:nn), lwd = 1)

dev.off()  

################################################################################ 

matR
