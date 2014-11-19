## Run this script to give a rough prediction of purity based on a dilution experiment

## REMEMBER TO USE YOUR OWN SOURCE-BUILT Rscript at /home/chris_w/apps/R-3.1.1/bin/Rscript

## Load libraries
library(VariantAnnotation)

## Gather commandline arguments
args = commandArgs(trailingOnly = TRUE)
inputFile=args[1] # path to input VCF
tumourname=args[2] # name of tumour sample

## Read in data
vcf = readVcf(file=inputFile,"b37")

## Find location of maximum peak in VAF plot
returnvaf = function(vcf){
    tumour=which(colnames(geno(vcf)$AD)==tumourname)
    depth=geno(vcf)$DP[,tumour]
    alt=sapply(geno(vcf)$AD[,tumour],"[",2)
    vaf=alt/depth
    return(vaf)
}
vcfdensity=density(returnvaf(vcf))
maxpeak=vcfdensity$x[which(vcfdensity$y==max(vcfdensity$y))]



## Construct linear model and make prediction
peaks=c(0.06379740,0.05019292,0.07195961,0.13958002,0.23181078,0.40310612)
purities=c(5,10,20,40,70,100)
model=lm(purities~peaks)
#purity=min(predict(model,newdata=data.frame(peaks=maxpeak)),100)
purity=predict(model,newdata=data.frame(peaks=maxpeak),interval="confidence")
filename=(paste0(tumourname,".txt"))
write.table(purity,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",file=filename)
