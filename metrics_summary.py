#!/usr/bin/env python2.7

## Pass this script an input directory and it will hunt for various sequencing metrics
## Note that it only works with a specific directory structure

import argparse
import logging
import os
import csv

## Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="Target input directory",required=True)
parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")

def main():
    ## Fetch/create results
    samplename=os.path.basename(args.i)
    logging.debug("Processing "+samplename)
    median=getmedian()
    duplication=getduplication()
    total=gettotal()
    mapping=getmapping()
    ontarget=getontarget()

    ## Print results
    print "\t".join([samplename,median,duplication,total,mapping,ontarget])

def getmedian():
    for row in csv.reader(open(os.path.join(args.i,"metrics/median.depth.txt")),delimiter="\t"):
	depth=str(row[0])
    return depth

def getduplication():
    rownumber=0
    for row in csv.reader(open(os.path.join(args.i,"dedup/picard.metrics")),delimiter="\t"):
	if rownumber==7:
	    duplication=str(row[7])
	    return duplication
	rownumber+=1

def gettotal():
    total=0
    for row in csv.reader(open(os.path.join(args.i,"metrics/offtarget.txt")),delimiter="\t"):
	if "reads" in row[0]:
	    total+=int(row[1])
    return str(total)

def getmapping():
    total=int(gettotal())
    for row in csv.reader(open(os.path.join(args.i,"metrics/offtarget.txt")),delimiter="\t"):
	if "Unaligned" in row[0]:
	    unaligned = int(row[1])
    return str(1-float(unaligned)/total)

def getontarget():
    for row in csv.reader(open(os.path.join(args.i,"metrics/offtarget.txt")),delimiter="\t"):
	if "percentage" in row[0]:
	    percentage = 1-(float(row[1]))/100
	    return str(percentage)

## Execute main method
if __name__ == '__main__':
    main()

