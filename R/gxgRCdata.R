gxgRCdata <-
function(X1,X2,Y){
  
  # Error checks
  if(length(X1)!=length(X2)){stop("Error: X1 and X2 must be the same length")}
  if(length(X1)!=length(Y)){stop("Error: X1 and X2 must be the same length as Y")}
  
  # Create a results matrix for each simulation
  mat_results <- matrix(0,nrow=1,ncol=8)
  colnames(mat_results)<-c("naiveE","naiveP","intE0","intP0","intE1","intP1","intE2","intP2")
  
  # fit1: No interaction
  fit1 <- lm(Y~X1+X2)
  mat_results[1,"naiveE"] <- summary(fit1)$coefficients[2,1]
  mat_results[1,"naiveP"] <- summary(fit1)$coefficients[2,4]
  
  # fit2: With interaction
  fit2 <- lm(Y~X1+X2+X1*X2)
  # For X2=0, beta1
  mat_results[1,"intE0"] <- summary(fit2)$coefficients[2,1]
  mat_results[1,"intP0"] <- summary(fit2)$coefficients[2,4]
  # For X2=1, beta1+betaI
  mat_results[1,"intE1"] <- summary(fit2)$coefficients[2,1]+summary(fit2)$coefficients[4,1]
  est1 <- summary(fit2)$coefficients[2,1] + summary(fit2)$coefficients[4,1]
  var1 <- vcov(fit2)[2,2] + vcov(fit2)[4,4] + 2*vcov(fit2)[2,4]
  t1 <- est1/sqrt(var1)
  pval1 <- 2*pt(abs(t1), summary(fit2)$df[2], lower=FALSE)
  mat_results[1,"intP1"] <- pval1
  # For X2=2, beta1+2*betaI
  mat_results[1,"intE2"] <- summary(fit2)$coefficients[2,1]+2*summary(fit2)$coefficients[4,1]
  est2 <- summary(fit2)$coefficients[2,1] + 2*summary(fit2)$coefficients[4,1]
  var2 <- vcov(fit2)[2,2] + 4*vcov(fit2)[4,4] + 4*vcov(fit2)[2,4]
  t2 <- est2/sqrt(var2)
  pval2 <- 2*pt(abs(t2), summary(fit2)$df[2], lower=FALSE)
  mat_results[1,"intP2"] <- pval2
  
  # Write Table
  list(mat_results)
  
}
