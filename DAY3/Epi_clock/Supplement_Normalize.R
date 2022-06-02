rm(list=ls());gc(TRUE)
options(stringsAsFactors = F)
setwd("/home/yunsung/third_project/data")
### Start renormalization
gse=commandArgs(TRUE)
cat("\n\n\nRenormalize : ",gse,"\n\n\n")

### Unzip the downloaded tar file
cat("Untar unzip the raw data:",gse,"\n")
if(paste(gse,"_idat",sep="") %in% list.files()==FALSE){
    system(paste("mkdir ",gse,"_idat",sep=""))
    system(paste("tar -xvf ",gse,"_RAW.tar -C ",gse,"_idat/",sep=""))
}

### Load info data.
require(feather)
targets=data.frame(read_feather(paste(gse,"_info.feather",sep="")))
cat("Display part of Basenames\n")
print(head(targets$Basename))


### Redefine an working directory
setwd(paste("/home/yunsung/third_project/data/",gse,"_idat/",sep=""))

### Load minfi and feather
suppressMessages(require(minfi));require(feather)

### Extract Sentrix ID from the extracted idat files
#files=list.files()
#idat_files=files[substr(files,nchar(files)-6,nchar(files)) == "idat.gz"]
#idat_files=substr(idat_files,1,nchar(idat_files)-12)
#freq=data.frame(table(idat_files))
#if(sum(freq$Freq!=2) == 0){cat("Good! All samples have both green and red.\n")}else{cat("Bad! Some samples are missing either green or #red.\n")}
#idat_files=unique(idat_files)
#cat("Sample size: ",length(idat_files),"\n")

### Create a temporal phenotypic file
#targets=data.frame(Basename=idat_files)
cat("Display part of targets\n")
print(head(targets[,c("Age","Sex","Basename")]))
#targets=data.frame(read_feather(path=paste("/home/yunsung/third_project/data/refined/",gse,"_info.feather",sep="")))

### Start normalization
cat("Create a RGset\n")
RGset   <- read.metharray.exp(base=NULL, targets = targets,recursive=TRUE, force=TRUE,verbose=TRUE)
#t1=estimateCellCounts(RGset, compositeCellType="CordBloodNorway")
cat("Create a Mset\n")
Mset <- preprocessNoob(RGset, verbose=TRUE) 
cat("Create a gset\n")
gset <- mapToGenome(Mset) 
cat("Create a grset\n")
grset=ratioConvert(gset,what="both")
cat("Extract beta values\n")
beta=getBeta(grset)
cat("Create expr\n")
expr=data.frame(ID=row.names(beta),beta)
cat("write expr\n")
write_feather(expr,path=paste(gse,"_expr_noob.feather",sep=""))


