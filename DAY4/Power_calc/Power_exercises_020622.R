## Power and sample size calculations, exercises 


## Exercise 1
## Fetal effects, asymptotic power and sample size

## a) 300 case triads are available for a single-SNP analysis. The minor allele frequency is 0.1. 
## 	  What is the power to detect a relative risk of 1.5 at the 5% significance level?

## In the following, you would like 90% power to detect a relative risk of 2 at the 5% significance level in a diallelic situation. 
## Assume a minor allele frequency of 0.05.

## b) How many case and control children are needed if you would like to sample twice as many controls as cases?

## c) Case and control parents are also genotyped. You would still like to sample twice as many controls as cases. 
##    What is the required number of case and control triads?

## d) How many individuals need to be genotyped in b) and c)? What is the most efficient design in this situation?


## Exercise 2
## Parent-of-origin (PoO) effects

## a) Execute the command below. What is the power to detect the given PoO effect applying 300 case-parent triads?

.power_PoO <- hapPowerAsymp(cases = c(mfc=300), haplo.freq = c(0.1,0.9), RRcm = c(1.5,1), RRcf = c(1,1))


## b) What is the power to detect the PoO effect if one choose to genotype 450 case-mother dyads instead?

## c) Which design, case-parent triads or case-mother dyads, is the most efficient in this situation?
