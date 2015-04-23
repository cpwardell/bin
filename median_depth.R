## Run this script after running using samtools to produce a pileup over an exome
## samtools mpileup -l exome.bed dedup.bam | awk '{print $4}' > pileup.txt 

## Gather commandline argumetns
args = commandArgs(trailingOnly = TRUE)
print(args)
inputfile=args[1]
outputname=args[2]
exome=args[3]

## Read in data
depth = read.table(inputfile,header=FALSE)[,1]
exome = read.table(exome,header=FALSE)

## Pileup only returns sites with coverage; how many sites have a 
## depth of zero?
exomelength=sum(exome[,3]-exome[,2])
shortfall=exomelength-length(depth)

## We can now calculate median depth
sorteddepthzeroes=c(rep(0,shortfall),sort(depth))
mediandepth=median(sorteddepthzeroes)

## Produce output file
write.table(mediandepth,file=paste0(outputname,".depth.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)


