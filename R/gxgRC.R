gxgRC <-
function(n=1000,nSim=5000,MAF1=0.5,gamma0=0,gammaX1=0.3,
beta0=0,betaX1=0.3,betaX2=0.3,betaI=seq(from=0.3,to=1,by=0.05),varY=1,
alpha_level=0.00000005,plot.pdf=T,plot.name="gxgRC.pdf",SEED=1){
  
##############################################################
#Set the seed
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
  colnames(mat_results) <- c("scenario1","scenario2","scenario3","scenario4","scenario5")
  
  #M0: Y~X1 Ho: no association between X1 and Y
  #M1: Y~X1+X2 Ho: no association between X1 and Y
  #M2: Y~X1+X2+X1*X2 Ho: no association between X1*X2 and Y
  
  #scenario 1: reject M0, M1, M2
  #scenario 2: reject M0, M1
  #scenario 3: reject M0, M2
  #scenario 4: reject M0
  #scenario 5: fail to reject M0
   
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
      # Generate  binary SNP of interest    
      X1 <- rbinom(n,1,MAF1) 
      
      #generate SNP 2 from SNP 1
      logitX2<-gamma0+gammaX1*X1
      pX2<-exp(logitX2)/(1+exp(logitX2))
      X2 <- rbinom(n,1,pX2)
       
      # Generate Y based on X1 and X2
      muY<-beta0+betaX1*X1+betaX2*X2+betaI[bi]*X1*X2
      Y <- rnorm(n,muY,sqrt(varY))

##############################################################
# Fit 3 models
##############################################################
	#M0: Y~X1 Ho: no association between X1 and Y
	m0 <- summary(lm(Y~X1))$coef[2,4] 
    
    #M1: Y~X1+X2 Ho: no association between X1 and Y
    m1 <- summary(lm(Y~X1+X2))$coef[2,4] 
    
    #M2: Y~X1+X2+X1*X2 Ho: no association between X1*X2 and Y
    m2 <- summary(lm(Y~X1+X2+X1*X2))$coef[4,4]  
      
##############################################################
# Determine scenarios
##############################################################
  #scenario 1: reject M0, M1, M2
 if((m0<alpha_level)&(m1<alpha_level)&(m2<alpha_level)){
 	mat_results[bi,"scenario1"]<-mat_results[bi,"scenario1"]+1}
   
  #scenario 2: reject M0, M1
 if((m0<alpha_level)&(m1<alpha_level)&(m2>alpha_level)){
 	mat_results[bi,"scenario2"]<-mat_results[bi,"scenario2"]+1}
   
  #scenario 3: reject M0, M2
 if((m0<alpha_level)&(m1>alpha_level)&(m2<alpha_level)){
 	mat_results[bi,"scenario3"]<-mat_results[bi,"scenario3"]+1}
   
  #scenario 4: reject M0
if((m0<alpha_level)&(m1>alpha_level)&(m2>alpha_level)){
	mat_results[bi,"scenario4"]<-mat_results[bi,"scenario4"]+1}
   
  #scenario 5: fail to reject M0
if((m0>alpha_level)){
	mat_results[bi,"scenario5"]<-mat_results[bi,"scenario5"]+1}
 
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
    plot(-1,-1, xlim=c(min(betaI),max(betaI)), ylim=c(0,1.05),
    xlab=expression(beta[I]),ylab="")
    points(betaI,mat_results[,"scenario1"],type="b",lty=2,col=1,pch=1)
    points(betaI,mat_results[,"scenario2"],type="b",lty=3,col=2,pch=2)
    points(betaI,mat_results[,"scenario3"],type="b",lty=4,col=3,pch=3)
    points(betaI,mat_results[,"scenario4"],type="b",lty=5,col=4,pch=4)
    points(betaI,mat_results[,"scenario5"],type="b",lty=6,col=6,pch=5)

legend("topleft",lty=c(2:6),col=c(1:4,6),pch=c(1:5),horiz=T,
legend=c("S1","S2","S3","S4","S5"))
    dev.off()
  }
  
  # Write Table
  list(mat_results) 
}
##############################################################
