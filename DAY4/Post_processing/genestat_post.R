## NORBIS GENESTAT COURSE

#############################
## POST-PROCESSING OF DATA ##
#############################

## Open project Genstat2022.Rproj
library(here)
library(Haplin)
library(ggplot2)
library(tidyverse)
## Bioconductor is not a "normal" package
## Installing qvalue
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.14") ## 3.14 for R version 4.1. R version 4.2 requires 3.15

setwd(here("DAY4","Post_processing"))
source(here("DAY4","Reading_data_in_Haplin","helper_functions.R"))

## Load data from before lunch
qc_all <- genDataLoad(filename = "qc_all_preproc", dir.in = "../Reading_data_in_Haplin")

## Look at data
qc_all

## Exercise 1
## a) What is the total number of individuals/families?
##
## b) What is the number of SNPs? 
##

## Run haplin
haplinRuns <- haplinSlide(data=qc_all,winlength=1,response="mult",reference="ref.cat")

## Class of result
class(haplinRuns)

## Length of result
length(haplinRuns)

## Works only when "winlength=1"
result <- toDataFrame(haplinRuns,reduce=T)

dim(result)
head(result)

## Get RR and p-value
pRR.result <- result[,c("marker","RR.est.","RR.p.value")]

head(pRR.result)

## Look at p-values
## x-axis: SNP names
## y-axis: p-value
## Cutoff at p=0.05
plotPValues( haplinRuns, plot.signif.only = TRUE, signif.thresh = 0.05 )

## Plotting RRs with confidence intervals and allele frequencies for SNPs with lowest p-values
plot(haplinRuns, plot.signif.only = TRUE, signif.thresh = 0.01)

#######################################################################
## q-values
## Storey & Tibshirani. Statistical significance for genomewide studies
## PNAS 2003
#######################################################################

## Loading qvalue
library(qvalue)
?qvalue

## Create q-values from the p-values
q <- qvalue(p=pRR.result$RR.p.value)

## Look at results
summary(q)

## Extract information
## Names of entries in q
names(q)

## Extract q-values and add to matrix of p-values and RRs
pRR.result <- cbind(pRR.result,q=q$qvalues)
head(pRR.result,20)

## Plot q-value according to p-value
## Low p-values correspond to low q-values
plot(x = pRR.result$RR.p.value, y = pRR.result$q, xlab = "p", ylab = "q",
     pch = 20, ylim = 0:1, xlim = 0:1)

###############################################################
## Note that several p-values correspond to the same q-value ##
## This is normal because of how the q-value is calculated   ##
## See Appendix B in Storey & Tibshirani for details         ##
###############################################################

## Common to consider SNPs with q<0.1
## SNPs with lowest q-values
pRR.sort <- pRR.result[order(pRR.result$RR.p.value),]
head(pRR.sort,20)


## Exercise
## a) Create new vector, ps, with the p-values from pRR.results
##
## b) Add new SNP with a p-value of 0.00001 to ps. Name the new SNP "rs1"
##
## c) Create new vector, qs, with the q-values from ps
##
## d) Write summary(qs) and compare with summary(q)
##
## e) Make a vector, q.sort, where qs is sorted according to ps, and 
##     compare with pRR.sort
##
## f) Add red points ("new" q-values vs. p-values) to the "old" plot


#######################################################################
## QQ-plot
#######################################################################
## X-axis: Expected p-value on -log10 scale. 
## Y-axis: Observed p-value on -log10 scale. 
## -log10 scale means that we count the number of zeros in the p-value
##  E.g., -log10(1)=0, -log(0.1)=1, ..., -log10(0.0000001)=7
##  Hence, large numbers mean small p-values
## Red line: What we get if the observed p-values are exactly as 
##            we would expect if H0 was true for ALL markers
## Points above line: Observed p-values are smaller than expected if 
##                     H0 was true for all markers
## Points below line: Observed p-values are larger than expected if 
##                     H0 was true for all markers
## We focus on the smallest p-values

## "Old" p-values
pQQ(ps[-length(ps)], lim = c(0,5.2))

## "New" p-values
pQQ(ps, lim = c(0,5.2))

## Note that the expected values of rs361 and the other SNPs jumped to
##  the red line when we added rs1. Why?

#######################################################################
## Volcano plot
#######################################################################
## X-axis: log2(RR)
## Y-axis: -log10(p-value)
## We focus on values with large RRs and low p-values
## These are found in the top corners
x <- log2(pRR.result$RR.est.)
y <- -log10(pRR.result$RR.p.value)
red <- y>1.5&abs(x)>.5
plot(x = x, y = y, xlab = "log2(RR)", ylab = "-log10(p)",col=1+(red),
     xlim = c(-.7,.7))
text(x = x[red], y = y[red],labels = pRR.result$marker[red],pos=4,col="red")

## Exercise
## a) Create new vector, RRs, with the RRs from pRR.results
##
## b) Add new SNP with RR=3 RRs. Name the new SNP "rs1", and make sure it 
##    is in the position of the corresponding p-value
##
## c) Create a volcano plot using RRs and ps. Do not add xlim or ylim. 



