## Run this script after running bedtools on your data like so:
## coverageBed -hist -abam bamfile -b exome.bed | grep ^all > all.only.txt

## Gather commandline argumetns
args = commandArgs(trailingOnly = TRUE)
print(args)
inputFile=args[1]
outputName=args[2]

## Read in data
depth = read.table(file = inputFile)

## THIS IS MEDIAN DEPTH
medianDepth = depth[which(cumsum(depth[,3]) >= sum(depth[,3])/2)[1],2]
write.table(medianDepth,file=paste0(outputName,".depth.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)

png(paste0(outputName,".depth.png"), h=1000, w=1000, pointsize=20)

plot(depth[, 2], depth[,5], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth",
     ylim=c(0,1.0),xlim=c(0,200), main=paste0("Target Region Depth, median = ",medianDepth))
abline(v = c(20,50,80,100), col = "gray60")
abline(h = c(0.50,0.90), col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.50,0.90), labels=c(0.50,0.90)) 
lines(1-cumsum(depth[,5]),lwd=3,col=rgb(7,142,83,max=255))

dev.off()
