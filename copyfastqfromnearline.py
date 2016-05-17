#!/usr/bin/env python

# We're going to take a file as input, then hunt for the appropriate files on nearline disk
# Then we're going to copy them to a local directory

## The input text file should have 3 columns:
## 1st is sample name, 2nd is string in fastq file, 3rd is flowcell
## HK02_c  HK_Exome_Lib_pool_Set_2_c_20131010_05   BC2FBFACXX


## Import modules
import logging
import argparse
import csv
import os
import subprocess

## Gather command line args
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Turn logging on",action="store_true")
parser.add_argument("-i", type=str, help="Path to input file",required=True)
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")
    logging.debug("Input is: "+args.i)

## Read input file
for row in csv.reader(open(args.i),delimiter="\t"):
    ## Attempt to create directory
    try:
	os.mkdir(row[0])
    except:
	pass
    
    ## Go hunting for the files we want
    potentials=[]
    fastqs=[]
    #for root, dirs, files in os.walk("/near/data/icgc_cgm/data/solexa"):
    for root, dirs, files in os.walk("/archive/data/icgc_cgm/premium_nearline/data/solexa"):
	for file in files:
	    if row[1] in file and "checksum" not in root:
		potentials.append(os.path.join(root,file))

    ## We need to traverse the results by flowcell, which will be the 8th element of the split
    ## If a flowcell directory contains "exome_dual", we keep only those results
    ## Otherwise, we keep everything
    flowcellswithdual=[]
    for potential in potentials:
	if "dual" in potential.split("/")[8]:
	    flowcellswithdual.append(potential.split("/")[7])

    ## Make the list unique
    flowcellswithdual = list(set(flowcellswithdual))

    ## We also want to get rid of "mismatch" data
    for potential in potentials:
	if potential.split("/")[7] in flowcellswithdual:
	    if "dual" in potential:
		fastqs.append(potential)
	elif "mis" not in potential:
	    fastqs.append(potential)

    ## Make the list unique
    fastqs = list(set(fastqs))
    
    ## Use rsync to copy the files to the new local directory
    for idx,fastq in enumerate(fastqs):
	filename=fastq.replace("/",".")
        rsynccommand = "rsync --progress -arv "+fastq+" "+row[0]+"/"+filename[1:len(filename)]
	logging.debug(rsynccommand)
	subprocess.call(rsynccommand,shell=True)


