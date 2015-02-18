## Load the PSCBS package
#install.packages("PSCBS")
library(PSCBS)

## Gather commandline arguments
args = commandArgs(trailingOnly = TRUE)
normal=args[1] # May be gzipped
tumour=args[2] # May be gzipped
output=args[3]

## Function to convert raw Illumina data to PSCBS data
pscbsmaker = function(normal,tumour){
      ## Read in data
        n=read.delim(normal,stringsAsFactors=FALSE)
	  t=read.delim(tumour,stringsAsFactors=FALSE)
	    
	      ## Clean data ; cast ratios to numeric, remove NA rows in BOTH and also unwanted chromosomes
	        n$Log.R.Ratio=as.numeric(n$Log.R.Ratio)
		  t$Log.R.Ratio=as.numeric(t$Log.R.Ratio)
		    n$B.Allele.Frequency=as.numeric(n$B.Allele.Frequency)
		      t$B.Allele.Frequency=as.numeric(t$B.Allele.Frequency)
		        badrows=unique(c(which(is.na(n$Log.R.Ratio)),which(is.na(t$Log.R.Ratio))))
			  badrows=unique(c(badrows,which(!n$Chr%in%c(1:22,"X","Y")),which(!t$Chr%in%c(1:22,"X","Y"))))
			    n=n[-badrows,]
			      t=t[-badrows,]
			          
				    ## Create output object; note that absolute copy number is calculated
				      ## using the log R ratio
				        output=as.data.frame(matrix(nrow=nrow(n),ncol=6))
					  rownames(output)=n$Name
					    colnames(output)=c("chromosome","x","CT","betaT","CN","betaN")
					      output$chromosome=n$Chr
					        output$x=n$MapInfo
						  output$CT=2*2^t$Log.R.Ratio
						    output$betaT=t$B.Allele.Frequency
						      output$CN=2*2^n$Log.R.Ratio
						        output$betaN=n$B.Allele.Frequency
							     return(output)
}
pnt = pscbsmaker(normal,tumour)

## Run the analysis
data = pnt
data = dropSegmentationOutliers(data) # OPTIONAL
gaps = findLargeGaps(data, minLength = 1e+06) # find gaps
knownSegments = gapsToSegments(gaps) # convert to segments
fit = segmentByPairedPSCBS(data, knownSegments = knownSegments,preserveScale = FALSE, seed = 48879,
verbose = -10)

#date() 
#fit <- callROH(fit, verbose = -10) # DO NOT CALL ROH; it interferes with AB calls
#date() # 6 minutes
deltaAB <- estimateDeltaAB(fit, scale=1)
#date() # 1 second
fit <- callAB(fit, delta=deltaAB)
#date() # 10 minutes
fit = callLOH(fit, delta=estimateDeltaLOH(fit, midpoint=0.7)) ## midpoint=0.5 is default value
#date() # instant
kappa <- estimateKappa(fit)
#date() # instant
deltaCN <- estimateDeltaCN(fit, scale=1, kappa=kappa)
#date() # instant
fit <- callNTCN(fit, delta=deltaCN, verbose=-10)
#date() # 1 minute

## Plot and write to file - doesn't work on HPC, R compiled without Cairo support
#png(file=paste0(output,".pscbs.plot.png"),width = 1000, height = 1000, units = "px")
#plotTracks(fit)
#dev.off()

# Get segments and remove NA rows
seggies = getSegments(fit, simplify = TRUE)
seggies=seggies[-which(is.na(seggies$chromosome)),]
seggies=seggies[-which(is.na(seggies$abCall)),]

## Write this to a file
write.table(seggies,file=paste0(output,".pscbs.segments.txt"),sep="\t",col.names=TRUE,
row.names=FALSE,quote=FALSE)

## Finally, save the session
save.image(paste0(output,".pscbs.RData"))

# "rohCall".  TRUE = run of homozygosity.  A germline LOH in the NORMAL sample
# "abCall" # allelic balance ; FALSE means a copy number event
# "lohCall" # TRUE means a loss of heterozygosity.  Could be uniparental disomy (UPD)
# "ntcnCall" # neutral copy number.  Only important for giving us ranges to compare to

## Copy number events are where abCall is true.
## Gains have tcnMean closer to the higher value of the ntcn range
## Losses have tcnMean closer to the lower value of the ntcn range

#plot(x=1:nrow(seggies),seggies$tcnMean,col=colorc(seggies$tcnMean),main="xxx")
#plot(x=1:150,seggies[seggies$chromosome=="1","tcnMean"],col=colorc(seggies$tcnMean),main="xxx")

### If we need to enforce a minimum window, we can use rle to find outliers
#distancetohigher=abs(seggies$tcnMean-fit$params$ntcnRange[2])
#distancetolower=abs(seggies$tcnMean-fit$params$ntcnRange[1])
#events=rep("neutral",nrow(seggies))
#events[distancetohigher<distancetolower & !seggies$abCall]="gain"
#events[distancetohigher>distancetolower & !seggies$abCall]="loss"
#seggies$events=events






