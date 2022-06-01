rm(list=ls());gc(TRUE)

### Set your working directory.
setwd("C://Users/YULE/OneDrive - Folkehelseinstituttet/Desktop/GENSTAT/")

### Load epigenetic clock.
beta_hat=as.data.frame(feather::read_feather("clock.feather"))

### Load phenotypic data (test set).
info_test=as.data.frame(feather::read_feather("GSE51032_info.feather"))

### Redefine chronological age.
info_test$Age=as.numeric(info_test$characteristics_ch1.1.age)
info_test$Sample_ID=substr(info_test$supplementary_file,68,95)

### Load DNAm data (test set).
expr_test=as.data.frame(feather::read_feather("GSE51032_expr_noob.feather"))
rownames(expr_test)=expr_test$ID;expr_test$ID=c();expr_test=as.matrix(expr_test);expr_test=t(expr_test)
expr_test=expr_test[match(info_test$Sample_ID,rownames(expr_test)),]
table(info_test$Sample_ID==rownames(expr_test))

### Validate the epigenetic clock.
info_test$p_Age=cbind(1,expr_test[,setdiff(beta_hat$cpgs,"(Intercept)")])%*% beta_hat$coef
windows(height=8,width=16);par(mfrow=c(1,2))
plot(info_test$Age, info_test$p_Age, xlab="Epigenetic age", ylab="Chronological age")
abline(0,1,col="red")
plot(info_test$p_Age, info_test$Age, ylab="Epigenetic age", xlab="Chronological age")
abline(0,1,col="red")

















