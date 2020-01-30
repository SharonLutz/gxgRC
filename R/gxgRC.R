gxgRC <-
function(n=1000,nSim=1000,MAF1=0.2,gamma0=0,gammaX1=0.2,
beta0=0,betaX1=0.1,betaX2=0.1,betaI=seq(from=0.1,to=0.5,by=0.1),varY=1,
alpha_level=0.05,plot.pdf=T,plot.name="gxgRC.pdf",SEED=1){
  
##############################################################
  set.seed(SEED)
 
##############################################################
# Error checks
##############################################################
  if(MAF1<0|MAF1>1){stop("Error: MAF1 must be between 0 and 1")}
  if(!(n>0)){stop("Error: n must be positive")}
  if(floor(n)!=ceiling(n)){stop("Error: n must be an integer")}
  if(!(nSim>0)){stop("Error: nSim must be positive")}
  if(floor(nSim)!=ceiling(nSim)){stop("Error: nSim must be an integer")}
  if(alpha_level<0|alpha_level>1){stop("Error: alpha_level must be between 0 and 1")}
  if(!(varY>0)){stop("Error: varY must be greater than 0")}
   
##############################################################
# Create a total results matrix
##############################################################
 mat_results <- matrix(0,nrow=length(betaI),ncol=5)
  colnames(mat_results) <- c("model0:X1","model1:X1","model2:X1&XI","model2:X1","model2:XI")
     
 
##############################################################
#cycle through simulations
##############################################################
  for(ii in 1:nSim){
if(floor(ii/100)==ceiling(ii/100)){print(paste("simulation",ii,"of",nSim))}
    
##############################################################
#cycle through betaI
##############################################################
    for(bi in 1:length(betaI)){
  
##############################################################
#Generate the data, X1 SNP of interest, X2 SNP in LD, outcome Y
##############################################################     
      # Generate  SNP of interest    
      X1 <- rbinom(n,1,MAF1) 
      
      #generate SNP 2 from SNP 1
      logitX2<-gamma0+gammaX1*X1
      pX2<-exp(logitX2)/(1+exp(logitX2))
      X2 <- rbinom(n,1,pX2)
       
      # Generate Y based on X1 and X2
      muY<-beta0+betaX1*X1+betaX2*X2+betaI[bi]*X1*X2
      Y <- rnorm(n,muY,sqrt(varY))
      
##############################################################
# Fit Model 1: No interaction
##############################################################
      model1 <- lm(Y~X1+X2) #summary(fit1F)$coefficients[2,4]
      if(summary(model1)$coefficients[2,4] < alpha_level){mat_results[bi,"model1:X1"] <- mat_results[bi,"model1:X1"] +1 }
      
##############################################################
# Fit model 2: With interaction
##############################################################
      fitR <- lm(Y~X2) #reduced model
      model2 <- lm(Y~X1+X2+X1*X2)
      if(anova(fitR,model2)$P[2] < alpha_level){mat_results[bi,"model2:X1&XI"] <- mat_results[bi,"model2:X1&XI"] +1 }
      if(summary(model2)$coefficients[2,4] < alpha_level){mat_results[bi,"model2:X1"] <- mat_results[bi,"model2:X1"] +1 }
      if(summary(model2)$coefficients[4,4] < alpha_level){mat_results[bi,"model2:XI"] <- mat_results[bi,"model2:XI"] +1 }
##############################################################
# Fit model 0: Just X1
##############################################################
model0 <- lm(Y~X1) 
      if(summary(model0)$coefficients[2,4] < alpha_level){mat_results[bi,"model0:X1"] <- mat_results[bi,"model0:X1"] +1 }
      
##############################################################
# end sims and save results
##############################################################
    } # End of betaI loop
    
  } # End of nSim
  
  # Divide by the number of simulations 
 mat_results <-mat_results/nSim
  
##############################################################
#make plot and output results
##############################################################
  if(plot.pdf){
    # Put plot code here
    pdf(plot.name)
    plot(-1,-1, xlim=c(min(betaI),max(betaI)), ylim=c(0,1),xlab="betaI values",ylab="")
    points(betaI,mat_results[,"model0:X1"],type="b",lty=2,col=1,pch=1)
    points(betaI,mat_results[,"model1:X1"],type="b",lty=3,col=2,pch=2)
    points(betaI,mat_results[,"model2:X1&XI"],type="b",lty=4,col=3,pch=3)
    points(betaI,mat_results[,"model2:X1"],type="b",lty=5,col=4,pch=4)
    points(betaI,mat_results[,"model2:XI"],type="b",lty=6,col=5,pch=5)
    legend("topleft",lty=c(2:6),col=c(1:5),pch=c(1:5),legend=c("model0:X1","model1:X1","model2:X1&XI","model2:X1","model2:XI"))
    dev.off()
  }
  
  # Write Table
  list(mat_results)
  
}
