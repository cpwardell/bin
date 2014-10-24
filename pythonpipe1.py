#!/usr/bin/env python2.7

## Primary sequencing pipeline to perform the following six tasks:
## 1.) Use cutadapt to remove adapters and quality trimming
## 2.) Alignment using BWA-MEM
## 4.) Merge and sort bams
## 5.) Deduplicate
## 6.) Produce metrics


## Import modules
import sys
import logging
import argparse
import os
import random
import subprocess

## Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="Full path to input FASTQ files.  May be specified multiple times",action="append",required=False)
parser.add_argument("-a", type=str, help="RG ID tag (data ID) - e.g. unique ID",required=False)
parser.add_argument("-p", type=str, help="1=cutadapt 2=split 3=align 4=sort 5=dedup 6=metrics",default="123456",required=False)
parser.add_argument("--wgs", help="Use for whole genome sequencing (WGS)",action="store_true")
parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")

## Global variables
uid = str(random.random())[2:8] # A 6-digit random integer
scratchDir="/home/chris_w/tmp/"
previousjob = "no_previous_job"
java = "/home/chris_w/apps/jdk1.8.0_20/bin/java"
python = "/home/chris_w/apps/Python-2.7.8/python"
picard = "/home/chris_w/apps/picard-tools-1.119/"
cutadaptbin="/home/chris_w/apps/cutadapt-1.5/bin/cutadapt"
fwdadapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
revadapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
bwa="/home/chris_w/apps/bwa-0.7.10/bwa"
refgenome="/home/chris_w/resources/b37/human_g1k_v37.fasta.gz"
intersectbed = "/home/chris_w/apps/bedtools2/bin/intersectBed"
coveragebed = "/home/chris_w/apps/bedtools2/bin/coverageBed"
rdepth = "/home/chris_w/bin/median_depth.R"
## For whole-genome runs, fake having a target so median-depth calculation using intersectBed still works
if args.wgs:
    exome = "/home/chris_w/resources/bedfiles/genomic_features/whole_genome_b37d5.bed"
else:
    #exome = "/home/chris_w/resources/bedfiles/agilent/sureselect/SureSelectV5_target_only.bed"
    exome = "/home/chris_w/resources/bedfiles/agilent/sureselect/SureSelectV5_target_only.ext100.bed"

## MAIN ##
def main():
    try:
	## Create directory for output and descend into it
	dircreate(args.a)
	
	## Execute the steps of the pipeline...
	logging.debug("Command string is: "+args.p)
	if "1" in args.p:
	    splitter()
	if "2" in args.p:
	    cutadapt()
	if "3" in args.p:
	    align()
	if "4" in args.p:
	    sort()
	if "5" in args.p:
	    dedup()
	if "6" in args.p:
	    metrics()

    except:
	logging.exception("Error in main method")
	sys.exit()
## MAIN END ##


## SPLIT AND DECOMPRESS ##
## Split only if WGS!!
## Otherwise, decompress the files...
def splitter():
    dircreate("split")
    for directory in args.i:
        fastqs=thesefiles(directory,".bz2")
	for fastq in fastqs:
	    if args.wgs:
                splitcommand = "bzcat "+directory+"/"+fastq+" | split -l 200000000 - "+fastq+"\n"+\
                "FILES=* ; for FILE in $FILES ; do mv $FILE $FILE.fastq ; done"
                logging.debug(splitcommand)
                jobsubmit(splitcommand,"split.sh","nowait")
            else:
                splitcommand = "bzcat "+directory+"/"+fastq+" > "+fastq[0:len(fastq)-4]
                logging.debug(splitcommand)
                jobsubmit(splitcommand,"split.sh","nowait")

    if args.wgs:
        sys.exit("WHOLE GENOME SAMPLE - wait until splitting complete...")

    os.chdir("..")
## SPLITTER END ##

## CUTADAPT ##
def cutadapt():
    dircreate("cutadapt")
    decompressedfiles=[]
    if args.wgs:
	decompressedfiles = os.listdir("../split/")
    else:
	for directory in args.i:
	    decompressedfiles += thesefiles(directory,".bz2")
    ones=[]
    twos=[]

    import sys

    for file in decompressedfiles:
        if "1.fastq" in file:
	    ones.append(file)
	else:
	    twos.append(file)
    ones.sort()
    twos.sort()
    
    ones=[one.replace(".bz2","") for one in ones]
    twos=[two.replace(".bz2","") for two in twos]

    for idx,one in enumerate(ones):
        logging.debug(one+"\t"+twos[idx])
	cutadaptcommand=python+" "+cutadaptbin+" -q 10 -a "+fwdadapter+\
	" --minimum-length 70 -f fastq -o "+one+\
	" -p "+twos[idx]+" ../split/"+one+" ../split/"+twos[idx]+"\n"+\
	python+" "+cutadaptbin+" -q 10 -a "+revadapter+\
	" --minimum-length 70 -f fastq -o "+one+".trimmed.gz -p "+\
	twos[idx]+".trimmed.gz "+one+" "+twos[idx]+"\n"+\
	"rm "+one+" "+twos[idx] 
	logging.debug(cutadaptcommand)
	## Submit job
	jobsubmit(cutadaptcommand,"cutadapt.sh",args.a+"split.sh*")

    ## Return to base dir
    os.chdir("..")
## CUTADAPT END ##

## ALIGNMENT ##
def align():
    dircreate("align")
    trimmedfiles=[]
    #trimmedfiles = os.listdir("../cutadapt/")
    if args.wgs:
	trimmedfiles = os.listdir("../split/")
    else:
	for directory in args.i:
	    trimmedfiles += thesefiles(directory,".bz2")
    ones=[]
    twos=[]
    for file in trimmedfiles:
	if "1.fastq" in file:
	    ones.append(file+".trimmed.gz")
        else:
	    twos.append(file+".trimmed.gz")
	ones.sort()
	twos.sort()

	ones=[one.replace(".bz2","") for one in ones]
	twos=[two.replace(".bz2","") for two in twos]

    for idx,one in enumerate(ones):
	logging.debug(one+"\t"+twos[idx])
	aligncommand=bwa+" mem -R \"@RG\\tID:"+args.a+"\\tSM:"+args.a+"\\t:ILLUMINA\" "+\
	refgenome+" ../cutadapt/"+one+" ../cutadapt/"+twos[idx]+" > "+str(idx)+".sam\n"
	logging.debug(aligncommand)
	## Submit job
	jobsubmit(aligncommand,"align.sh",args.a+"*cutadapt.sh*")
    ## Return to base dir
    os.chdir("..")

    #$BWA mem -R "@RG\tID:$OUTPUT\tSM:$OUTPUT\tPL:ILLUMINA" $INDEX $FORWARD $REVERSE > $OUTPUT.sam
## ALIGNMENT END ##

## SORT ##
def sort():
    ## Create directory to work in and descend into it
    dircreate("sort")
    sams=[]
    ## Get list of files to merge and or sort
    if args.wgs:
	sams = thesefiles("../align/",".sam")
    else:
	for directory in args.i:
	    fastqs = thesefiles(directory,".bz2")
	countofsams=len(fastqs)/2
	for count in range(countofsams):
	    sams.append(str(count)+".sam")
    
    ## Merge and sort the sams
    inputblock=""
    for sam in sams:
	inputblock += " INPUT=../align/"+sam

    ## WE NEED ONE FILE PER PAIR OF INPUTS

    sortcommand = java+" -Djava.io.tmpdir="+scratchDir+\
    " -Xmx4g -jar "+picard+\
    "MergeSamFiles.jar OUTPUT=merged.sorted.bam SORT_ORDER=coordinate "+\
    "CREATE_INDEX=TRUE USE_THREADING=false MAX_RECORDS_IN_RAM=1000000 "+\
    "VALIDATION_STRINGENCY=SILENT TMP_DIR="+scratchDir+inputblock
    logging.debug("Sort command: "+sortcommand)
    
    ## Submit job
    jobsubmit(sortcommand,"sort.sh",args.a+"align.sh*")

    ## Return to base dir
    os.chdir("..")	
## SORT END##

## DEDUPLICATE ##
def dedup():
    ## Create directory to work in and descend into it
    dircreate("dedup")

    ## Deduplicate
    dedupcommand = java+" -Djava.io.tmpdir="+scratchDir+\
    " -Xmx4g -jar "+picard+\
    "MarkDuplicates.jar METRICS_FILE=picard.metrics "+\
    "REMOVE_DUPLICATES=false ASSUME_SORTED=true CREATE_INDEX=true "+\
    "MAX_RECORDS_IN_RAM=1000000 VALIDATION_STRINGENCY=SILENT "+\
    "TMP_DIR="+scratchDir+" INPUT=../sort/merged.sorted.bam OUTPUT=dedup.bam"

    ## Submit job
    jobsubmit(dedupcommand,"dedup.sh")

    ## Return to base dir
    os.chdir("..")

## DEDUPLICATE END ##

## METRICS ##
def metrics():
    ## Create directory to work in and descend into it
    dircreate("metrics")

    ## Median depth    
    mediandepthcommand = coveragebed+" -hist -abam "+\
    "../dedup/dedup.bam"+\
    " -b "+exome+" | grep ^all > depth.ALL.txt\n"+\
    "Rscript "+rdepth+" depth.ALL.txt median"
    jobsubmit(mediandepthcommand,"mediandepth.sh")
    
    ## Only perform off-target calculation if there is a target exome
    if not args.wgs:
        ## Off-target sequencing
        offtargetcommand = "NOTALIGNED=$(samtools view ../dedup/dedup.bam | awk '{if($3==\"*\"){print}}' | wc -l)\n"+\
        "ONTARGET=$("+intersectbed+\
        " -sorted -abam ../dedup/dedup.bam -b "+exome+" | samtools view - | wc -l)\n"+\
        "RAWOFFTARGET=$("+intersectbed+" -sorted -v -abam ../dedup/dedup.bam -b "+exome+\
        " | samtools view - | wc -l)\n"+\
        "OFFTARGET=$(echo \"$RAWOFFTARGET-$NOTALIGNED\" | bc)\n"+\
        "PERCENT=$(echo \"scale=2;$OFFTARGET*100/($ONTARGET+$OFFTARGET)\" | bc)\n"+\
	"echo -e \"Unaligned reads:\t$NOTALIGNED\" > offtarget.txt\n"+\
        "echo -e \"On target reads:\t$ONTARGET\" >> offtarget.txt\n"+\
        "echo -e \"Off target reads:\t$OFFTARGET\" >> offtarget.txt\n"+\
        "echo -e \"Off target percentage:\t$PERCENT\" >> offtarget.txt\n"
        jobsubmit(offtargetcommand,"offtarget.sh")

    ## Return to base dir
    os.chdir("..")
## METRICS END ##

## Writes a command to a shell script and submits it
def jobsubmit(command,scriptname,groupname=None):
    ## Define name of this job
    global previousjob
    thisjob = args.a+scriptname+uid
    ## Write and submit script
    script = open(scriptname,"w")
    script.write("#!/bin/bash\n\n"+command+"\n")
    if(groupname==None):
	submission = "qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 -N "+thisjob+" -hold_jid "+previousjob+" "+scriptname
    else:
	submission = "qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 -N "+thisjob+" -hold_jid	"+groupname+" "+scriptname
    script.close()
    logging.debug(submission)
    subprocess.call(submission,shell=True)

    ## Store the name of this job for the next time...
    previousjob = thisjob
    

## Return list of files with a certain file extension
def thesefiles(location,suffix):
    results=[]
    for file in os.listdir(location):
	if file.endswith(suffix):
	    results.append(file)
    return(results)

## Create a directory if one doesn't already exist and go into it
def dircreate(dirname):
    if not os.path.isdir(dirname):
	os.mkdir(dirname)
    os.chdir(dirname)
    ## Report where we are
    logging.debug(os.getcwd())

## Execute main method
if __name__ == '__main__':
	main()

