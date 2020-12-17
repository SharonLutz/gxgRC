## gxgRC


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
gxgRC(n=1000,nSim=1000,MAF1=0.5,gamma0=0,gammaX1=0.3,
beta0=0,betaX1=0.3,betaX2=0.3,betaI=seq(from=0.3,to=1,by=0.05),varY=1,
alpha_level=0.00000005,plot.pdf=T,plot.name="gxgRCexample.pdf",SEED=1)
```

## Simulation Scenario Output
For this example, we get the following matrix and corresponding plot:

```
      scenario1 scenario2 scenario3 scenario4 scenario5
 [1,]     0.000     0.956         0     0.027     0.017
 [2,]     0.000     0.991         0     0.008     0.001
 [3,]     0.010     0.983         0     0.002     0.005
 [4,]     0.027     0.973         0     0.000     0.000
 [5,]     0.059     0.941         0     0.000     0.000
 [6,]     0.133     0.867         0     0.000     0.000
 [7,]     0.214     0.786         0     0.000     0.000
 [8,]     0.330     0.670         0     0.000     0.000
 [9,]     0.505     0.495         0     0.000     0.000
[10,]     0.650     0.350         0     0.000     0.000
[11,]     0.803     0.197         0     0.000     0.000
[12,]     0.878     0.122         0     0.000     0.000
[13,]     0.933     0.067         0     0.000     0.000
[14,]     0.976     0.024         0     0.000     0.000
[15,]     0.990     0.010         0     0.000     0.000
```
<img src="https://github.com/SharonLutz/gxgRC/blob/master/gxgRCexample.png" width="500">

