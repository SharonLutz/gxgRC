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

## Example
For 1,000 subjects, we generate X1 to have a MAF of 0.20 and X2 to have an MAF of 0.05. The interaction between X1 and X2 on Y varies from 0, 0.1, to 0.2.

```
gxgRC(n = 1000, betaB = 0.1, beta0 = 0, beta1 = 0.1, beta2 = 0.1, betaI = c(0,0.1, 0.2), MAF1 = 0.2, 
MAF2 = 0.05, varY = 1, alpha_level = 0.05, plot.pdf = T, plot.name = "gxgRC.pdf",nSim = 1000, SEED=1)
```

## Output
For this example, we get the following matrix and corresponding plot:

```
     NOintP intP0 intP1 intP2
[1,]   0.30  0.12  0.07  0.04
[2,]   0.53  0.15  0.26  0.08
[3,]   0.79  0.19  0.74  0.46
```
<img src="https://github.com/SharonLutz/gxgRC/blob/master/gxgRC.png" width="600">
