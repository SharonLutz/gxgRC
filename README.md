## gxgRC
The gxgRC R package computes the power to determine the association of a SNP of interest with the outcome given that another SNP in the region is in weak LD with the SNP of interest.

## Installation
```
install.packages("devtools") # devtools must be installed first and only once
```

```
devtools::install_github("SharonLutz/gxgRC") #install for new updates
```

## Input
For n subjects, SNP X1 is generated from a binomial distribution with a MAF specified by the user (input: MAF1). The second SNP X2 is generated to have an association with X1 specified by the user (input: betaB) and to have an MAF specified by the user (input: MAF2). The outcome Y is generated from a normal distribution with mean as follows:

E\[Y\] = &beta;<sub>0</sub> + &beta;<sub>1</sub> X<sub>1</sub> + &beta;<sub>2</sub> X<sub>2</sub> + &beta;<sub>I</sub> X<sub>1</sub> X<sub>2</sub>   

See the manpage for more detail regarding the input of the gxgRC function.

```
library(gxgRC)
?gxgRC # For details on this function
```

## Simulation Scenario
For 1,000 subjects, we generate X1 to have a MAF of 0.20 and X2 to have an MAF of 0.05. The interaction between X1 and X2 on Y varies from 0, 0.1, to 0.2.

```
gxgRC(n=1000,nSim=1000,MAF1=0.2,gamma0=0,gammaX1=0.2,
beta0=0,betaX1=0,betaX2=0.2,betaI=seq(from=0.1,to=0.5,by=0.1),varY=1,
alpha_level=0.05,plot.pdf=T,plot.name="gxgRC.pdf",SEED=1)
```

## Simulation Scenario Output
For this example, we get the following matrix and corresponding plot:

```
     model0:X1 model1:X1 model2:X1&XI model2:X1 model2:XI
[1,]     0.135     0.112        0.123     0.047     0.097
[2,]     0.330     0.271        0.381     0.055     0.244
[3,]     0.560     0.500        0.699     0.058     0.461
[4,]     0.807     0.764        0.912     0.051     0.693
[5,]     0.933     0.919        0.990     0.055     0.878
```
<img src="https://github.com/SharonLutz/gxgRC/blob/master/gxgRCplot.png" width="500">

