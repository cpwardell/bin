#!/usr/bin/env python2.7

## Purpose; a re-write of our idtest script originally in R.  Should be
## much faster in Python using pysam

import argparse
import sys
import logging
import csv
import pysam
import os

## Define some globals
#listofsnps = "/home/chris_w/resources/idtestsnps/6668.id.snps.b37.txt"
listofsnps = "/home/chris_w/resources/idtestsnps/6668.id.snps.b37.tiny.txt"

## Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("-n", type=str, help="full path to normal BAM",required=False)
parser.add_argument("-t", type=str, help="full path to tumour BAM",required=False)
parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")

args.n="/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_L/dedup/dedup.bam"
args.t="/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_c/dedup/dedup.bam"

def main():
    try:
	
	## Overwrite existing output files with empty files
	#clearoutput(args.n)
	#clearoutput(args.t)
	## Iterate through the list of SNPs to check..
	snpiterate()
	


    except:
	logging.exception("Error in main")
	sys.exit()

def clearoutput(path):
    file=os.path.dirname(os.path.abspath(path))[0:-5]+"idtest.txt"
    logging.debug("Deleting previous file: "+file)
    if os.path.exists(file):
	os.remove(file)

def basecounter(bam,CHROM,POS):
    try:
	## Dict to hold counts of each base
	basecount={"A":0,"C":0,"G":0,"T":0,"N":0}
	## Open the bam file
	samfile=pysam.Samfile(bam,"rb") # rb = "read bam"
	for alignedread in samfile.fetch(CHROM,POS,POS+1):
	    if(alignedread.is_proper_pair):
		## Which base in the read is at the position we want?
		## Use the "aligned_pairs" list of tuples to determine this
		offset = [item for item in alignedread.aligned_pairs if item[1] == POS][0][0]
		## offset == None when there is an indel at the site of the SNV
		if(offset!=None):
		    basecount[alignedread.seq[offset]]+=1
		else:
		    basecount["N"]+=1
	return(basecount)
    except:
	logging.exception("Error in basecounter")
	sys.exit()

def snpiterate():
    try:
	## A dictionary to store all the results in
	judgements={"MATCH":0,"MISMATCH":0,"FAIL":0}
	for row in csv.reader(open(listofsnps),delimiter="\t"):
	    logging.debug(row)
	    CHROM=row[0]
	    POS=int(row[1])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1
	    allele1=row[3]
	    allele2=row[4]
	    snpid=row[5]

	    ## Fetch counts of each base at SNV location
	    normalcounts=basecounter(args.n,CHROM,POS)
	    tumourcounts=basecounter(args.t,CHROM,POS)

	    ## Perform some calculations
	    normaldepth=sum(normalcounts.values())
	    tumourdepth=sum(tumourcounts.values())
	    ## Do not bother with sites with depth < 10 in either tumour or normal
            ## or sites where the suggested alleles don't account for at least 90% of the reads
	    if(normaldepth<10 or tumourdepth<10 or \
		normalcounts[allele1]+normalcounts[allele2]/normaldepth < 10 or \
		tumourcounts[allele1]+tumourcounts[allele2]/tumourdepth < 10):
		pass
	    else:
		nallele=allelecall(normalcounts,allele1,allele2)
		tallele=allelecall(tumourcounts,allele1,allele2)
		judgement=alleletest(nallele,tallele)
		judgements[judgement]+=1
		## Write raw data to normal and tumour dirs
		#writeallele(args.n,snpid,nallele,allele1,allele2)
		#writeallele(args.t,snpid,tallele,allele1,allele2)
	## Write summary to tumour dir
	#writejudgement(judgements)
	print judgements

    except:
	logging.exception("Error in snpiterate")
	sys.exit()

def writejudgement(judgements):
    total=sum(judgements.values())
    fpc = float(judgements["FAIL"])/total
    mispc = float(judgements["MISMATCH"])/total
    matpc = float(judgements["MATCH"])/total
    outfile=os.path.dirname(os.path.abspath(args.t))[0:-5]+"idteststats.txt"
    out = open(outfile,"w")
    out.write("MATCH:\t"+str(judgements["MATCH"])+"\t"+str(matpc)+"\n"+\
    "MISMATCH:\t"+str(judgements["MISMATCH"])+"\t"+str(mispc)+"\n"+\
    "FAIL:\t"+str(judgements["FAIL"])+"\t"+str(fpc)+"\n")
    

def writeallele(path,snpid,allele,allele1,allele2):
    outfile=os.path.dirname(os.path.abspath(path))[0:-5]+"idtest.txt"
    
    if allele=="HOM1":
	alleles=allele2+"/"+allele2
    elif allele=="HOM2":
	alleles=allele1+"/"+allele1
    else:
	alleles=allele1+"/"+allele2
    results = open(outfile,"a")
    results.write(snpid+"\t"+allele+"\t"+alleles+"\n")
    results.close()

def allelecall(counts,allele1,allele2):
    ratio=float(counts[allele1])/(counts[allele1]+counts[allele2])
    if(ratio<0.25):
	allele="HOM1"
    elif(ratio>0.75):
	allele="HOM2"
    else:
	allele="HET"
    return allele

def alleletest(sample1,sample2):
    if(sample1==sample2):
	result="MATCH"
    elif(sample1=="HOM1" and sample2=="HOM2" or sample1=="HOM2" and sample2=="HOM1"):
	result="FAIL"
    elif(sample1=="HOM1" and sample2=="HET" or sample1=="HET" and sample2=="HOM1" or \
	 sample1=="HOM2" and sample2=="HET" or sample1=="HET" and sample2=="HOM2"):
	result="MISMATCH"
    return result


## Execute 
if __name__ == '__main__':
    main()
