# Get libraries for multiprocessors to speed things up
require(parallel)

# Set seed so that we get the same resuts
set.seed(1234)

# Create a fake DNA methylation matrix: 10000 CpGs (columns) and 100 samples (rows)
dnam<-matrix(runif(10000*100, min=0, max=1), ncol=10000, nrow=100)

# Create a fake cell type matrix as given by the old "Houseman procedure"
celltype<-matrix(runif(6*100, min=0, max=1), ncol=6, nrow=100)

# Make fake CpG labels as given by the Illumina platform
cpglab<-paste("cg0",sep="",formatC(1:10000, digits=6, flag="0"))

# Fake cell type labels
celllab<-c("NK", "CD4", "CD8", "Gran", "Mono", "Bcell")


# Add CpG and cell type labels as column names to their respective matrices
colnames(dnam)<-cpglab
colnames(celltype)<-celllab

############ Make fake phenotypes and covariates

# Linear slope with Gaussian random noise
phenotype<-1:100+rnorm(100, mean=0, sd = 6)

# Sex covariate selected from a Bernoulli distribution with 50% chance of male/female
cov_sex<-rbinom(100, 1,0.5)

# Smoking variable also selected from a Bernoulli distribution with 10% chance of being a smoker
cov_smk<-rbinom(100, 1,0.1)

# Age covariate chosen from a Poisson distribution with mean age 35
cov_age<-rpois(100, 35)


# Function that performs regression to be called in parallel from each core/processor
# y - column from each DNA methylation matrix
# x - phenotype
# c1 - cov_sex covariate
# c2 - cov_smk covariate
# c3 - cov_age covariate
# celltype - celltype matrix
#
# The function returns a row consiting of slope parameter (est) standard error (se) and p value (pval)

lmFun=function(y,x, c1, c2, c3, cellmatr){
	mod=lm(y~x+c1+c2+c3+cellmatr)
	coefs=summary(mod)$coeff
	est=coefs[2,1]
	se=coefs[2,2]
	pval=coefs[2,4]
	out=c(est, se, pval)
	names(out)=c("est", "se", "pval")
	return(out)
}

# Start cluster process, assuming 8 threads (can be more but may be slower)
cl <- makeCluster(getOption("cl.cores", 8))

# Call parallel apply function (parApply) extracting columns from dnam (2)
# parApply calls lmFun above with covariates cov_sex, cov_smk, cov_age and celltype matrix
# parApply is transposed (t) so that the data will fit into a CpG X regression parameter matrix

reg_res=t(parApply(cl, dnam ,2 ,lmFun, phenotype, cov_sex, cov_smk, cov_age, celltype))

# Should parallel fail:
# reg_res=t(apply( dnam ,2 ,lmFun, phenotype, cov_sex, cov_smk, cov_age, celltype))

# When finished remember to stop cluster process otherwise threads could be running and slow things down
stopCluster(cl)

# Convert matrix returned from parApply call to data frame for easier handling
reg_res<-data.frame( reg_res )

# Add CpG labels to results matrix so we know which CpG is significant
row.names( reg_res )<-cpglab

# Sort with respect to p values (lowest on top)
reg_res_sort<-reg_res[order(reg_res$pval),]

# Only include p < 0.05
reg_res_sort[reg_res_sort$pval<0.05,]

# Mixed-effects model also useful, but not straightforward

lmeFun=function(y,x, c1, c2, c3, cellmatr){
	mod=lme(y~x+c1+c2+c3, random = ~ 1|c3)
	coefs=summary(mod)$tTable
	est=coefs[2,1]
	se=coefs[2,2]
	pval=coefs[2,5]
	out=c(est, se, pval)
	names(out)=c("est", "se", "pval")
	return(out)
}

stopCluster(cl)
cl <- makePSOCKcluster(getOption("cl.cores", 8))
# Export library to all cores
clusterEvalQ(cl, library(nlme))
# As well as all objects used...
clusterExport(cl, c("dnam", "phenotype","cov_sex", "cov_smk", "cov_age", "celltype"))
# Run as previously
reg_res=t(parApply(cl, dnam, 2 ,lmeFun, phenotype, cov_sex, cov_smk, cov_age, celltype))
stopCluster(cl)
# Convert matrix returned from parApply call to data frame for easier handling
reg_res<-data.frame( reg_res )

# Add CpG labels to results matrix so we know which CpG is significant
row.names( reg_res )<-cpglab

# Sort with respect to p values (lowest on top)
reg_res_sort<-reg_res[order(reg_res$pval),]

# Only include p < 0.05
reg_res_sort[reg_res_sort$pval<0.05,]
