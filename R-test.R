## This script lets you install the packages necessary to participate in GENESTAT
## Unfortunately, we will not have time to help you "live" during the lectures. 
## Therefore, it is very important that you ensure that this script runs
## successfully BEFORE the first lecture on May 30

## Installs necessary packages
install.packages(pkgs = c("Haplin",
	"dplyr",
	"tidyr",
	"readr",
	"feather",
	"glmnet",
	"doParallel"))
if (!requireNamespace("BiocManager", quietly = TRUE)){
	install.packages("BiocManager")
}
BiocManager::install("snpStats", update = FALSE)

## Loading packages and checking that they are loaded correctly
## Loading
library(Haplin)
library(dplyr)
library(tidyr)
library(readr)
library(glmnet)
## Checking
## packageVersion should return version the version of the package
## If packageVersion returns error or warning, please consult your nearest R expert
packageVersion("Haplin")
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("readr")
packageVersion("glmnet")

## If any of the below lines give warnings or errors, please consult your nearest R expert
## Plotting some data based on data stored in mtcars
plot(mpg~hp, data = mtcars, col = cyl, lwd = 2)

## Checking that the RTools-dependency works
exmpls.dir <- system.file("extdata", package = "Haplin")
exmpl.file1 <- file.path(exmpls.dir, "HAPLIN.trialdata.txt")

my.haplin.data <- genDataRead(file.in = exmpl.file1,
  file.out = "my_haplin_data", dir.out = ".", format = "haplin", n.vars = 0)

rm(my.haplin.data)
unlink(list.files(pattern = "my_haplin_data*"))
