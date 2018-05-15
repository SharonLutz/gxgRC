## gxgRC
The gxgRC R package computes power analyses for gene by gene interactions of both rare and common variants.

## Installation
```
install.packages("devtools") # devtools must be installed first

devtools::install_github("SharonLutz/gxgRC")
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
For 1,000 subjects, we generate X1 to have a MAF of 0.2, we generate X2 to have an MAF of 0.2, and to have an association of 0.15 with X1. The interaction between X1 and X2 on Y is generated at 0.1 and 0.2.

```
library(gxgRC)
gxgRC(n=1000,betaB=0.15,betaI=c(0.1,0.2),MAF1=0.2,MAF2=0.2,alpha_level=0.05,plot.pdf=T,
plot.name="gxgRC.pdf",nSim=100)
```

## Output
For this example, we get the following matrix and corresponding plot:

```
     naiveP intP0 intP1 intP2
[1,]  0.970 0.646 0.907 0.599
[2,]  0.981 0.681 0.901 0.577
```
<img src="https://github.com/SharonLutz/gxgRC/blob/master/gxgRC.png" width="600">
