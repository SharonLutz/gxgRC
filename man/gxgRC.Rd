\name{gxgRC}
\alias{gxgRC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gxgRC
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gxgRC(n = 1000, betaB = 0.1, beta0 = 0, beta1 = 0.2, beta2 = 0.2, betaI = c(0.1, 0.2), MAF1 = 0.2, MAF2 = 0.2, varY = 1, alpha_level = 0.05, plot.pdf = T, plot.name = "gxgRC.pdf", nSim = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
is the sample size
}
  \item{betaB}{
is the association between X1 and X2
}
  \item{beta0}{
is the intercept for Y
}
  \item{beta1}{
is the association between X1 and Y
}
  \item{beta2}{
is the association between X2 and Y
}
  \item{betaI}{
is the interaction between X1 and X2 on Y
}
  \item{MAF1}{
is the MAF for X1
}
  \item{MAF2}{
is the MAF for X2
}
  \item{varY}{
is the variance for Y
}
  \item{alpha_level}{
is the significance level
}
  \item{plot.pdf}{
is T to output a plot, is F to not output a plot
}
  \item{plot.name}{
is the name of the plot
}
  \item{nSim}{
is the number of simulations
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a matrix of estimates and p-values
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
