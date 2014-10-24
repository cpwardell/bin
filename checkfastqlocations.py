#!/usr/bin/env python

# We're going to take a file as input, then hunt for the appropriate files on nearline disk
# Then, print if the input is correct or not

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

    sample=row[0]
    id=row[1]
    flowcell=row[2]

    ## Go hunting for the files we want
    potentials=[]
    fastqs=[]
    for root, dirs, files in os.walk("/near/data/icgc_cgm/data/solexa"):
	for file in files:
	    if flowcell in root and id in file and "checksum" not in root:
		potentials.append(os.path.join(root,file))

    print sample+"\t"+flowcell+"\t"+str(len(potentials))

