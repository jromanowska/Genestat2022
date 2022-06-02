rm(list=ls());gc(TRUE)

### Install R packages.
#install.packages("feather")
#install.packages("glmnet")
#install.packages("doParallel")

### Set your working directory.
setwd("C://Users/YULE/OneDrive - Folkehelseinstituttet/Desktop/GENSTAT/")

###Load phenotypic data.
info=as.data.frame(feather::read_feather("GSE87571_info.feather"))

### Redefine variables.
info$characteristics_ch1.1.age[info$characteristics_ch1.1.age=="NA"]=NA
info$Age=as.numeric(info$characteristics_ch1.1.age)
info$Sex=info$characteristics_ch1.0.gender
info$Sample_ID=substr(info$supplementary_file,68,95)

### Plot the distribution of chronological age.
windows();hist(info$Age,xlab="Chronological age",main="")

### Remove samples with missing chronological age.
info=subset(info,is.na(info$Age)==F)

### Select portion of data (for shorter processing time).
set.seed(1)
info=info[sample(1:nrow(info),300),]

### Load DNAm data (already background corrected using the noob method).
expr=as.data.frame(feather::read_feather("GSE87571_expr_noob.feather"))
print(head(expr[,1:4]))

### Restructure the DNAm data.
rownames(expr)=expr$ID;expr$ID=c();expr=as.matrix(expr);expr=t(expr)
expr=expr[match(info$Sample_ID,rownames(expr)),]
table(info$Sample_ID==rownames(expr))

### Perform x-fold cross-validation.
doParallel::registerDoParallel(cores = 3)
set.seed(1)
cv_result=glmnet::cv.glmnet(x=expr,y=info$Age,
                            nfolds = 3, alpha=0.5, trace.it = T, standardize=T,parallel = T)
windows();plot(cv_result)

### Fit glmnet one more time.
set.seed(1)
result=glmnet::glmnet(x=expr,y=info$Age,
                            alpha=0.5, trace.it = T, standardize=T,lambda = cv_result$lambda.min)

beta_hat=data.frame(as.matrix(coef(result)))
beta_hat$cpgs=rownames(beta_hat);colnames(beta_hat)[1]="coef";beta_hat=beta_hat[,c("cpgs","coef")]
beta_hat=subset(beta_hat,beta_hat$coef!=0)
head(beta_hat);dim(beta_hat)

### Calculate epigenetic age in the training set.
info$p_Age=cbind(1,expr[,setdiff(beta_hat$cpgs,"(Intercept)")]) %*% beta_hat$coef
windows();plot(info$p_Age, info$Age, xlab="Epigenetic age", ylab="Chronological age")
abline(0,1,col="red")

feather::write_feather(beta_hat,"clock.feather")

