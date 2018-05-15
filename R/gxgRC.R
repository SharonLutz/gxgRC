gxgRC <-
function(n=1000,betaB=0.1,beta0=0,beta1=0.2,beta2=0.2,betaI=c(0.1,0.2),MAF1=0.2,MAF2=0.2,varY=1,
                  alpha_level=0.05,plot.pdf=T,plot.name="gxgRC.pdf",nSim=1000){
  
  # Error checks
  if(MAF1<0|MAF1>1){stop("Error: MAF1 must be between 0 and 1")}
  if(MAF2<0|MAF2>1){stop("Error: MAF2 must be between 0 and 1")}
  if(!(n>0)){stop("Error: n must be positive")}
  if(floor(n)!=ceiling(n)){stop("Error: n must be an integer")}
  if(alpha_level<0|alpha_level>1){stop("Error: alpha_level must be between 0 and 1")}
  if(!(varY>0)){stop("Error: varY must be greater than 0")}
  
  # Create a total results matrix
  mat_total <- matrix(0,nrow=length(betaI),ncol=4)
  colnames(mat_total) <- c("naiveP","intP0","intP1","intP2")
  
  for(i in 1:nSim){
    
    # Create a results matrix for each simulation
    mat_results <- matrix(0,nrow=length(betaI),ncol=4)
    colnames(mat_results) <- c("naiveP","intP0","intP1","intP2")
    
    for(bI.ind in 1:length(betaI)){
      
      # Generate the data: X1, X2 and Y
      X2 <- rbinom(n,2,MAF2) # Generate replicated SNP
      X2[X2==2] <- 1
      
      # Function=0 to find the root to get beta0 based on P(X1=1)=MAF1^2+2*MAF1*(1-MAF1)
      # P(X1=1)=P(X1=1|X2=1)P(X2=1)+P(X1=1|X2=0)P(X2=0)
      betaB0F<-function(betaB0){
        # PX11GX2=P(X1=1|X2)=exp(beta0+beta1*X2)/(1+exp(beta0+beta1*X2))
        PX11GX20<-exp(betaB0+betaB*0)/(1+exp(betaB0+betaB*0))
        PX11GX21<-exp(betaB0+betaB*1)/(1+exp(betaB0+betaB*1))
        PX11GX21*(MAF2^2+(2*MAF2*(1-MAF2)))+PX11GX20*(1-MAF2)^2-(MAF1^2+2*MAF1*(1-MAF1)) #P(X1=1)-P(X1=1)=0
      }
      
      betaB0 <<- uniroot(betaB0F,c(-1000,10))$root
      PX11GX20 <- exp(betaB0+betaB*0)/(1+exp(betaB0+betaB*0)) #P(X1=1|X2=0)
      PX11GX21 <- exp(betaB0+betaB*1)/(1+exp(betaB0+betaB*1)) #P(X1=1|X2=1)
      
      X1 <- rep(1,n)
      ind2 <- runif(n,0,1) # Generate a random uniform to generate X1
      X1[ind2>PX11GX20&X2==0] <- 0 #length(X1[X2==0&X1==1])/length(X1[X2==0])
      X1[ind2>PX11GX21&X2==1] <- 0 #length(X1[X2==1&X1==1])/length(X1[X2==1])
      
      # Generate Y based on X1 and X2
      Y <- rnorm(n,beta0+beta1*X1+beta2*X2+betaI*X1*X2,sqrt(varY))
      
      # fit1: No interaction
      fit1 <- lm(Y~X1+X2)
      if(summary(fit1)$coefficients[2,4] < alpha_level){mat_results[bI.ind,"naiveP"] <- mat_results[bI.ind,"naiveP"] +1 }
      
      # fit2: With interaction
      fit2 <- lm(Y~X1+X2+X1*X2)
      # For X2=0, beta1
      if(summary(fit2)$coefficients[2,4] < alpha_level){mat_results[bI.ind,"intP0"] <- mat_results[bI.ind,"intP0"] +1 }
      
      # For X2=1, beta1+betaI
      # var(ax+by)=a^2var(x) + b^2var(y) + 2abcov(x,y)
      # var(beta1+betaI)=var(beta1) + var(betaI) + 2cov(beta1,betaI)
      est1 <- summary(fit2)$coefficients[2,1] + summary(fit2)$coefficients[4,1]
      var1 <- vcov(fit2)[2,2] + vcov(fit2)[4,4] + 2*vcov(fit2)[2,4]
      t1 <- est1/sqrt(var1)
      pval1 <- 2*pt(abs(t1), summary(fit2)$df[2], lower=FALSE)
      if(pval1 < alpha_level){ mat_results[bI.ind,"intP1"] <- mat_results[bI.ind,"intP1"] +1 }
      
      # For X2=2, beta1+2*betaI
      # var(ax+by)=a^2var(x) + b^2var(y) + 2abcov(x,y)
      # var(beta1+2betaI)=var(beta1) + 4var(betaI) + 4cov(beta1,betaI)
      est2 <- summary(fit2)$coefficients[2,1] + 2*summary(fit2)$coefficients[4,1]
      var2 <- vcov(fit2)[2,2] + 4*vcov(fit2)[4,4] + 4*vcov(fit2)[2,4]
      t2 <- est2/sqrt(var2)
      pval2 <- 2*pt(abs(t2), summary(fit2)$df[2], lower=FALSE)
      if(pval2 < alpha_level){ mat_results[bI.ind,"intP2"] <- mat_results[bI.ind,"intP2"] +1 }
      
    } # End of betaI loop
    
    # Add the current matrix to the total matrix
    mat_total <- mat_total + mat_results
    
  } # End of nSim
  
  # Divide by the number of simulations to find power
  mat_total <- mat_total/nSim
  
  if(plot.pdf){
    # Put plot code here
    pdf(plot.name)
    plot(-1,-1, xlim=c(min(betaI),max(betaI)), ylim=c(0,1),xlab="betaI values",ylab="")
    points(betaI,mat_total[,"naiveP"],type="b",lty=2,col=1,pch=1)
    points(betaI,mat_total[,"intP0"],type="b",lty=3,col=2,pch=2)
    points(betaI,mat_total[,"intP1"],type="b",lty=4,col=3,pch=3)
    points(betaI,mat_total[,"intP2"],type="b",lty=5,col=4,pch=4)
    legend("topleft",lty=c(2:5),col=c(1:4),pch=c(1:4),legend=c("naiveP","intP0","intP1","intP2"),cex = 0.6)
    dev.off()
  }
  
  # Write Table
  list(mat_total)
  
}
