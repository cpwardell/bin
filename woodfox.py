#!/usr/bin/env python2.7

## If the ASCII art isn't readable below, you're not viewing this file in the intended way; 
## you MUST use a monospaced font
# ........... .. . .  . ........ . ... ...........................................
# ..........................................?MMMM$................................
# .............       . ................MMMMMMMMMMMMMM............................
# ..................... .............,MMMMMMMMMMMMMMMMMM..........................
# .................................NMMMMMMMMMMMMMMMMMMMMMM........................
# ...............................MMMMMMMMMMMMMMMMMMMMMMMMMMN... .. .. ... . ... ..
# .........M...................:MMMMMMMMMMMMMMMMMMMMMMMMMMMMN ....................
# ........ MM.................MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMI....................
# ........~MM$..............:MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM......... ..........
# ........MMMM............ MMMMMMMMMMMM  MMMMMMMMMMMMMMMMMMMMMD...................
# ........M.=M............MMMMMMMMMM. OMMMMM  +MMMMMMMMMMMMMMMM...................
# .......M ..MM   ...... MMMMMMMMM.............. OMMMMMMMMMMMMMM..................
# .......,... MMMMMMM8...MMMMMMMM ................ MMMMMMMMMMMMM .................
# ......M, MMMMMMMMMMMMMMMMMMMMM.......... N.........NMMMMMMMMMMMM,...............
# .....MMNMMMMMMMMMMMMMMMMMMMMM......... MMMM,........MMMMMMMMMMDOMMMMMNN.. ......
# ....MMMMMMMMMMMMMMMMMMMMMMMM8........ZMM .?M.........MMMMMMMMMM... ZMMMMI.......
# ....MMMMMMMMMMMMMMMMMMMMMMMM........MMO....NN........$M MMMMMMM........OMN .....
# ....MMMMMMMMMMMMM MMMMMMMMMM...... MM........NM.......D.8MMMMMM..........MD ....
# ....MMMMI..MMMMM.?MMMMMMMMMM.... ~M7..........~MM.... ...MMMMMM....M. ...MM ....
# ....MMMM.........MMMMMMMMMMM$...MM..............MM... ..,MMMMMM....IM ....MM....
# ...NMMMM..........MMMMMMMMMMMM MM ...............MM.....MMMMMMM.....M..... M....
# ...MMMM ..........MMMMMMMMMM? MM................. MM....MMMMMM:.....MM.....M~...
# ..MM...:MMMMMMD. MMMMMMMMMN MZ...................MM....MMMMMM .....MM.....N....
# ..MM...~M$....$M..MMMMMMMMM N=....................MM... MMMMM$,.....MMMM........
# ..M... MM.........:MMMMMMMM8M.....................MM...MMMMMMOM+....DMMM$...M...
# ...DNMM=...........MMMMMMMMM......................MM...MMMMMMMMMM...MMMMZ...M ..
# . MMMMM............ MMMMMMMM......................MM...MMMMM+,MM....MMMM ...M8..
# DMMM...............M=MMMMMMM......................MM  ..MMMMM $ ...IMMMMMO.:M...
# .MM?...............M MMMMMMM ......................MM, .MMMMM~.... MMMMMM..MM...
# .. .. . . .DM8? ...MM.MMMMMMZ........................MM .MMMMM... MMMMMMMOMM$...
# ............MMMMMN.MMD MMMMMM......................NM:.MMMMMMM:.=MMMMMMMMMMM....
# .............MMMMMMMMM.MMMMMM,.MM.................MMMMM MMMMMMMMMMMMMMMMMMM.....
# ............. MMMMMMMMM MMMMMMMMM ..............MMMMMMMMMMMMMMMMMMMMMMMMMM:.....
# ..... . .......MMMMMMMM.MMMMMMMMMO...........M?. $MMMMMMMMMMMMMMMMMMMMMMM ......
# ................MMMMMMMM MMMMMMMMM.............8MMMMNMMMMMMMMMMMMMMMMMM$........
# .....  .   ..... MM. . ..MMMMMMMMM .......M MMMMMMMMMMMMMMMMMMMMMMMMMM .........
# .................MM8......MMMMMMMMM..... MMMMMMMMMMMMMMMMMMMMMMMMMMZ ...........
# ....... . ....... MM......=MMM.MMMM.....MMMMM.:MMM...$MMMMMMMMMM$ ..............
# ..................MM.......MM..MMMM.....MM ..,MMM....MMMMM......................
# ....... .  ... . ..MM ..........MMMM....$... MMM....MMMMM.......................
# ...................MM ..........MMMM........MMM....MMMMM:MI.....................
# .......... .  .  . ?M,... .......MMM .... .MMM... MMMMMMMD......................
# ............   . ...ZZ...........MMMM.....NM.....:MMMMMMM.......................
# ......... ... .. .... ............MMM.....M .....MM MMMM........................
# ..................................MMM~.... ........MMMM.........................
# ....... ...            .. .... . ..MMM... .........MMZ..........................
# ..................... .............,MM........... MM ...........................
# ....... .             ........ .... MM ..........$M ............................
# ...........            .. .... . .. .M .. ...... N..............................
# ...........    ..  .  ........ . .....M.........................................
# ...........            ..  ... . ..     .  .....................................
# .....   .              ..  .   . .    ...  . ...................................

## Purpose: woodFox is our all-purpose indel filtering tool.
## Accepts either an exome directory as input, or full paths to the Indelocator vcf and bam files
## It calculates the following metrics with these limits:

## Mapping quality; x >  50 : calculated using Pysam - COMPLETE
## Base quality; > 26 : calculated using Pysam - COMPLETE
## Not somatic (i.e. in normal); no more than 2 indel containing reads allowed : calculated with Pysam - COMPLETE
## Alignability at site must be 1 - using Pybedtools - COMPLETE

import sys # so we can exit the program
import os # to check if directories exist
import argparse # command line args
import csv # line-by-line operations
import pysam # for native bam operations
import re # for regular expressions
import pybedtools # for intersectBed alignability
import socket # Required to find hostname
import logging # For debugging purposes

## Mean function so we don't have to import numpy
def mean(numbers):
    x=float(sum(numbers))/len(numbers)
    return(x)

## Gather command line args
## ONLY THE FIRST ARGUMENT (args.t) ; BY DEFAULT THE PROGRAM INFERS args.t1,args.t2
## Supply args.t OR args.t1,args.t2, NOT both
parser = argparse.ArgumentParser()
parser.add_argument("-t", type=str, help="full path to tumour exome directory; e.g. \"/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_c\"",required=False)
parser.add_argument("-t1", type=str, help="full path to tumour bam; e.g. \"/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_c/dedup/dedup.bam\"",required=False)
parser.add_argument("-t2", type=str, help="full path to parsed Strelka VEP file; e.g.\"/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_c/strelka/passed.somatic.indels.vep.parsed.txt\"",required=False)
parser.add_argument("-n", type=str, help="full path to normal exome directory; e.g. \"/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_L\"",required=False)
parser.add_argument("-n1", type=str, help="full path to normal bam; e.g. \"/home/chris_w/project_bile_duct_cancer/analysed_data/HK02_L/dedup/dedup.bam\"",required=False)
parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()


## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    ## Use statements like this to print to STDOUT
    logging.debug("Debugging mode enabled")

## EXAMPLE INDELS TO FILTER
## /mnt/HPC_haem/cwardell/analysis/exomes/My9-3_My9-CD138+-1/indels/indels.filtered.vcf
#args.t = "/home/chris_w/project_bile_duct_cancer/analysed_data/HK66_c"
#args.n = "/home/chris_w/project_bile_duct_cancer/analysed_data/HK66_L"

## If args.t is supplied, infer args.t1,args.t2
if(args.t!=None):
	## If the input argument doesn't end in a forward slash, add one
	if(args.t[-1] != "/"):
		args.t = args.t+"/"
	args.t1 = args.t + "dedup/dedup.bam"
	args.t2 = args.t + "strelka/passed.somatic.indels.vep.parsed.txt"					

## If args.n is supplied, infer args.n1
if(args.n!=None):
	## If the input argument doesn't end in a forward slash, add one
	if(args.n[-1] != "/"):
		args.n = args.n+"/"
	args.n1 = args.n + "dedup/dedup.bam"

## If args.t1,args.t2 don't exist by this point, they haven't been explicitly supplied 
## or inferred, so warn the user and exit
if(args.t1==None):
#if(args.t1==None and args.t2==None):
	print "You must supply EITHER the -t OR -t1,-t2 arguments"
	sys.exit()

## Check that the args.t1 directory exists.  If it doesn't, exit
if not os.path.isfile(args.t1):
	print args.t1+" does not exist; are you certain that this is a tumour sample and that Strelka has been run?"
	sys.exit()

## Echo arguments for debugging purposes
logging.debug("args.n: "+args.n)
logging.debug("args.t: "+args.t)


## Get coordinates of all indels and store them in an object - a list of tuples
indels = []
inOrDel=[]
headRE=re.compile("^#") # skip header
for row in csv.reader(open(args.t2),delimiter="\t"):
	if not headRE.match(row[0]):
		## Determine if the indel is an insertion or deletion
		type="insertion"
		if row[2]=="-":
		    type="deletion"
		## Append data to relevant object
		inOrDel.append(type)
		indels.append((row))






## Look for reads at the same position in the normal sample; how many indel containing reads are there?
logging.debug("Checking normal sample for indels")
somaticKeepers=[]
somaticindels=[]
## Iterate through every indel - this loop ONLY considers reads with indels in them
## We use enumerate() to create a nice index for us to use
for idx,indel in enumerate(indels):
	## Set properties of indel
	CHROM=indel[1].split(":")[0]
	POS=int(indel[1].split(":")[1].split("-")[0])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1

	## Count of indel-containing reads in normal
	indelReads = 0	
		
	## Open the normal bam file
	samfile=pysam.Samfile(args.n1,"rb") # rb = "read bam"
	
	for alignedread in samfile.fetch(CHROM,POS,POS+1):
		## Note that "is_proper_pair" excludes reads that map to different chromosomes (i.e. involved in interchromosomal translocations)
		if(alignedread.is_proper_pair):
			
			## Only detects deletions
			if(any(pair[0] is None for pair in alignedread.aligned_pairs) and inOrDel[idx]=="deletion"):
				indelReads+=1

			## Only detects insertions
			if(any(pair[1] is None for pair in alignedread.aligned_pairs) and inOrDel[idx]=="insertion"):
				indelReads+=1
			
	## No more than 2 indel containing reads allowed		
	if(indelReads < 3):
		somaticKeepers.append(idx)

	somaticindels.append(indelReads)

	
## Define a list in which to store the indices of all the lines we want to keep
logging.debug("Calcuating mapping and base quality scores")
mapbaseKeepers=[]
mapqualities=[]
basequalities=[]
## Loop through the list of indels again, this time considering base and mapping quality scores
## We use enumerate() to create a nice index for us to use
for idx,indel in enumerate(indels):
	
	## Set properties of indel
	CHROM=indel[1].split(":")[0]
	POS=int(indel[1].split(":")[1].split("-")[0])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1

	## Open the bam file
	samfile=pysam.Samfile(args.t1,"rb") # rb = "read bam"
	
	## Lists in which to store mapping and base quality data
	mapq=[]
	baseq=[]
	
	## Get values at the SNV location
	for alignedread in samfile.fetch(CHROM,POS,POS+1):
		if(alignedread.is_proper_pair):
			## Which base in the read is at the position we want?  Use the "aligned_pairs" list of tuples to determine this
			offset = [item for item in alignedread.aligned_pairs if item[1] == POS][0][0] 

			## offset == None when there is an indel at the site of the SNV
			if(offset!=None):			
				mapq.append(alignedread.mapq)
				baseq.append(ord(alignedread.qual[offset])-33) ## Subtract 33 because SAM specification tells us to
	
	if(mean(mapq) >= 50 and mean(baseq) > 26 ):
		mapbaseKeepers.append(idx)

	mapqualities.append(mean(mapq))
	basequalities.append(mean(baseq))


## Use Pybedtools to investigate alignability
logging.debug("Calculating alignability")

alignKeepers=[]
## We put this in a try clause, as some chromosome names (MT, GL, etc) are
## not in the alignability file and generate errors
alignfile=pysam.Tabixfile("/home/chris_w/resources/bedfiles/genomic_features/alignability/wgEncodeCrgMapabilityAlign100mer.bedGraph.gz")
alignabilities=[]
for idx,indel in enumerate(indels):
    ## Set properties of indel
    CHROM=indel[1].split(":")[0]
    POS=int(indel[1].split(":")[1].split("-")[0])-1 # Pysam coordinates are ZERO-based, so we MUST subtract 1
    try:
	for record in alignfile.fetch(CHROM, POS, POS+1):
	    alignability=float(record.split("\t")[3])
	    print alignability
	    alignabilities.append(alignability)
    except:
	alignabilities.append(0)
	pass

# Now we can write the output
logging.debug("Writing output")

## Iterate through every indel
iter=0
for row in csv.reader(open(args.t2),delimiter="\t"):
    tabrow = "\t".join(row)
    # Skip all header rows
    if(row[0].startswith("#Uploaded")):
	print "\t".join([tabrow,"somaticindels","mean_mapq","mean_baseq","alignability"])
    if(row[0].startswith("##")):
	print tabrow
    
    if not row[0].startswith("#"):
	print "\t".join([tabrow,str(somaticindels[iter]),str(mapqualities[iter]),str(basequalities[iter]),str(alignabilities[iter])])
	iter+=1

print somaticindels
print mapqualities
print basequalities
print alignabilities

