#!/usr/bin/env python2.7

## Secondary (temporary?) pipeline to perform variant calls and other 
## "real" analysis

## Import modules
import argparse
import logging
import subprocess
import pysam
import os
import random
import sys

## Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("-n", type=str, help="Full path to input normal bam",required=True)
parser.add_argument("-t", type=str, help="Full path to input normal bam",required=True)
parser.add_argument("-p", type=str, help="Command string details here...",default="123",required=False)
parser.add_argument("--wgs", help="Use for whole genome sequencing (WGS)",action="store_true")
parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")

## Declare some globals
uid = str(random.random())[2:8] # A 6-digit random integer
java="/home/chris_w/apps/jre1.7.0_67/bin/java"
mutectjar = "/home/chris_w/apps/mutect-src/mutect/target/mutect-1.1.7.jar"
reference="/home/chris_w/resources/b37/human_g1k_v37.fasta"
refchroms="/home/chris_w/resources/b37/chroms/"
cosmic="/home/chris_w/resources/cosmic/CosmicCodingMuts.vcf"
dbsnp="/home/chris_w/resources/b37/dbsnp_138.b37.vcf"
tempdir="/home/chris_w/tmp/"
previousjob = "no_previous_job"
tumourrg=""
idtestscript="/home/chris_w/bin/idtest.py"
strelkadir="/home/chris_w/apps/strelka/"

## List of chromosomes:
chroms=map(str,range(1,23))
chroms=chroms+["X","Y"]

# Steps to implement:
# ID test - COMPLETE
# SNV calling with MuTect - COMPLETE
# indel calling with Strelka - COMPLETE
# Tx detection

## Fetch RG tag from bam file using pysam
def getrg(bam):
    try:
	samfile=pysam.Samfile(bam,"rb") # rb = "read bam"
	rg=samfile.header["RG"][0]["ID"]
	return(rg)
    except:
	logging.exception("Error in getrg")

## Create a directory if one doesn't already exist and go into it
def dircreate(dirname):
    if not os.path.isdir(dirname):
	os.mkdir(dirname)
    os.chdir(dirname)
    ## Report where we are
    logging.debug(os.getcwd())

## Writes a command to a shell script and submits it
def jobsubmit(command,scriptname,groupname=None):
    ## Define name of this job
    global previousjob
    thisjob = tumourrg+scriptname+uid
    
    ## Write and submit script
    script = open(scriptname,"w")
    script.write("#!/bin/bash\n\n"+command+"\n")
    ## Jobs in this pipeline don't wait; they just execute ASAP
    ## Strelka jobs are submitted with default RAM, but 6 cores for speed
    if scriptname == "strelka.sh":
	qsub = "qsub -cwd -S /bin/bash -pe def_slot 6 -N "
    else:
	qsub = "qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 -N "
    submission = qsub+thisjob+" "+scriptname
    script.close()
    logging.debug(submission)
    subprocess.call(submission,shell=True)

    ## Store the name of this job for the next time...
    previousjob = thisjob

def mutect():
    try:
	dircreate("mutect")

	## MuTect command.  Note additional environment variable export
	## to override the Java memory limits
	## Note additional command which filters for only passing results
	## This command also annotates using VEP
	## Also, it executes individually for each chromosome
	for chrom in chroms:
	    mutectcommand = "export JAVA_TOOL_OPTIONS=\"-XX:+UseSerialGC -Xmx4g -Djava=\""+\
	    tempdir+" ; "+java+\
	    " -Xmx4g -Djava.io.tmpdir="+tempdir+" -jar "+mutectjar+\
	    " --analysis_type MuTect "+\
	    " --reference_sequence "+reference+\
	    " --cosmic "+cosmic+\
	    " --dbsnp "+dbsnp+\
	    " --input_file:normal "+args.n+\
	    " --input_file:tumor "+args.t+\
	    " -vcf mutect."+chrom+".vcf"+\
	    " --out call_stats."+chrom+".out "+\
	    " --only_passing_calls "+\
	    " -L "+chrom+" ; "+\
	    "/home/chris_w/apps/ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl -i mutect."+chrom+".vcf --cache --offline --everything --force_overwrite -o mutect."+chrom+".vep.txt ; "+\
	    " source /home/chris_w/apps/virtualpythonenvironment/bin/activate ; "+\
	    " /home/chris_w/bin/metalfox.py -f1 call_stats."+chrom+".out -f3 "+args.t+" > call_stats."+chrom+".postfox.out ; "+\
	    " /home/chris_w/bin/vep_txt_parser.py -i mutect."+chrom+".vep.txt > mutect."+chrom+".vep.parsed.txt ; "+\
	    " deactivate "
	    logging.debug(mutectcommand)
	    jobsubmit(mutectcommand,"mutect.sh")
	os.chdir("..")
    except:
	logging.exception("Error in mutect")


def idTest():
    try:
	## Note that we load a Python virtual environment!
	idtestcommand="source /home/chris_w/apps/virtualpythonenvironment/bin/activate ; "+\
	idtestscript+" -n "+args.n+" -t "+args.t+" ; "+\
	"deactivate"
	logging.debug(idtestcommand)
	jobsubmit(idtestcommand,"idtest.sh")

    except:
	logging.exception("Error in ID test")
	sys.exit()

def strelka():
    try:
	## No need to create output directory, Strelka will do this for us
	strelkacommand=strelkadir+"/bin/configureStrelkaWorkflow.pl"+\
        " --normal="+args.n+\
	" --tumor="+args.t+\
        " --ref="+reference+\
	" --config="+strelkadir+"etc/strelka_config_bwa_exomes.ini"+\
        " --output-dir=strelka ; "+\
	" cd strelka ; make -j 6 ; "+\
        " mv results/* . ; rm -rf Makefile chromosomes config results ; "+\
        " /home/chris_w/apps/ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl "+\
	" -i passed.somatic.indels.vcf "+\
	" --cache --offline --everything -o passed.somatic.indels.vep.txt ; "+\
	" /home/chris_w/bin/vep_txt_parser.py -i passed.somatic.indels.vep.txt > passed.somatic.indels.vep.parsed.txt ; "+\
	" /home/chris_w/bin/woodfox.py -n1 args.n -t1 args.t -t2 ../strelka/passed.somatic.indels.vep.parsed.txt "+\
	" > passed.somatic.indels.vep.parsed.postfox.txt " 
	logging.debug(strelkacommand)
        jobsubmit(strelkacommand,"strelka.sh")
    except:
	logging.exception("Error in strelka")

def main():
    global tummourrg
    try:
    	logging.debug("Doing ID test for "+args.t+" vs "+args.n)
	tumourrg=getrg(args.t)
	## Create directory for output and descend into it
	dircreate(tumourrg)
	
	## Execute the steps of the pipeline...
	logging.debug("Command string is: "+args.p)
	if "1" in args.p:
	    idTest()
	if "2" in args.p:
	    mutect()
	if "3" in args.p:
	    strelka()
	
    except:
	logging.exception("Error in main")
	sys.exit()

## Execute 
if __name__ == '__main__':
	main()
