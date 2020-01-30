\name{gxgRC}
\alias{gxgRC}
\title{gxgRC}
\description{Tests for the effect of a SNP on the outcome in the presence of epistasis.}
\usage{
gxgRC(n=1000,nSim=1000,MAF1=0.2,gamma0=0,gammaX1=0.2,
beta0=0,betaX1=0.1,betaX2=0.1,betaI=seq(from=0.1,to=0.5,by=0.1),varY=1,
alpha_level=0.05,plot.pdf=T,plot.name="gxgRC.pdf",SEED=1)
}
\arguments{
  \item{n}{
is the sample size
}
    \item{nSim}{
is the number of simulations
}
    \item{MAF1}{
is the MAF for X1
}
  \item{gamma0}{
is the intercept for the association between the SNP X1 and the SNP X2 where logit(X2)=gamma0+gammaX1*X1
}
    \item{gammaX1}{
is the effect of SNP X1 on SNP X2 where logit(X2)=gamma0+gammaX1*X1
}
  \item{beta0}{
is the intercept for E[Y]=beta0+betaX1*X1+betaX2*X2+betaI*X1*X2
}
  \item{betaX1}{
is the association between X1 and Y for E[Y]=beta0+betaX1*X1+betaX2*X2+betaI*X1*X2
}
  \item{betaX2}{
is the association between X2 and Y for E[Y]=beta0+betaX1*X1+betaX2*X2+betaI*X1*X2
}
  \item{betaI}{
is the interaction between X1 and X2 on Y for E[Y]=beta0+betaX1*X1+betaX2*X2+betaI*X1*X2
}
  \item{varY}{
is the variance for Y
}
  \item{alpha_level}{
is the significance level, default 0.05
}
  \item{plot.pdf}{
equals T to output a plot, equals F to not output a plot
}
  \item{plot.name}{
is the name of the plot
}
  \item{SEED}{
sets the seed
}
}
\details{
Tests for the effect of a SNP on the outcome in the presence of epistasis.
}
\value{
a matrix of estimates and p-values
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sharon Lutz
}


