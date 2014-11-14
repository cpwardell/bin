#!/home/chris_w/apps/Python-2.7.8/python

## Import modules
import logging
import argparse
import vcf
import sys
from operator import itemgetter

## Gather command line args
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Turn logging on",action="store_true")
parser.add_argument("--print_all", help="Output ALL annotations",action="store_true")
parser.add_argument("-v", type=str, help="full path to input vcf",required=True)
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
	logging.basicConfig(level=logging.DEBUG)
	logging.debug("Debugging mode enabled")
	logging.debug("Input VCF is: "+args.v)

## Define functions

# Casts lists to ints, with empty strings becoming zero
def mk_int(s):
    s = s.strip()
    return int(s) if s else 0

# Takes a list of transcripts and outputs appropriately
# If args.print_all is true, it will output every transcript as a new line
# Otherwise, will select a single "best" transcript and output one per mutation
def transcript_picker(transcripts):
	if(args.print_all):
		for transcript in transcripts:
			pipeprint(transcript)
		return
	## Only one transcript?  Write and return
	if numtranscripts == 1:
		pipeprint(transcripts[0])
		return
	## Multiple transcripts; make a choice	
	## First, collect which variants are WITHIN transcripts, which affect
	## canonical transcripts and the distances from transcripts
	nodistance=[]
	canonical=[]
	distances=[]
	for idx,transcript in enumerate(transcripts):
		pipesplit=pipesplitter(transcript)
		distances.append(mk_int(pipesplit[11]))
		if pipesplit[11] == "":
			nodistance.append(idx)
		if pipesplit[17] == "YES":
			canonical.append(idx)
	## Output the first annotation set that is within a canonical transcript
	## Otherwise, output the first annotation set that is within a transcript
	## Otherwise, output the annotation set that is closest to a transcript
	nodistance_and_canonical=list(set(nodistance).intersection(set(canonical)))
	if len(nodistance_and_canonical) > 0:
		pipeprint(transcripts[nodistance_and_canonical[0]])
	elif len(nodistance) > 0:
		pipeprint(transcripts[nodistance[0]])
	else:
		minimumdistanceindex = min(enumerate(distances), key=itemgetter(1))[0]
		pipeprint(transcripts[minimumdistanceindex])

# Takes a transcript string and splits it using the pipe symbol
# Returns a list
def pipesplitter(transcript):
	pipesplit=str.split(transcript,"|")
	return pipesplit
# Takes list and returns a string joined with tab characters
def pipejoiner(pipesplit):
	pipejoin="\t".join(pipesplit)
	return pipejoin
# Take a transcript and prints a full set of output
def pipeprint(transcript):
	pipesplit=str.split(transcript,"|")
	pipejoin="\t".join(pipesplit)
	print chr+"\t"+pos+"\t"+ref+"\t"+alt+"\t"+deptht+"\t"+refcountt+"\t"+altcountt+"\t"+\
	depthn+"\t"+refcountn+"\t"+altcountn+"\t"+pipejoin

## Open VCF file for reading
## Pipes (|) are not automatically parsed, so we do that ourselves
## First, open the file for reading and parse out the vep fields
## Join these into a tab-delimited string for easy printing
logging.debug("Accessing VCF header")
vcf_reader = vcf.Reader(open(args.v, 'r'))
vepfields = str.split(vcf_reader.infos['CSQ'][3])[-1]
vepfieldssplit = str.split(vepfields,"|")
vepfieldsheader = "\t".join(vepfieldssplit)
print "chr\tpos\tref\talt\tdeptht\trefcountt\taltcountt\tdepthn\trefcountn\taltcountn\t"+vepfieldsheader
logging.debug("Parsing VCF records")
for record in vcf_reader:
	chr=str(record.CHROM)
	pos=str(record.POS)
	ref=record.REF
	alt=str(record.ALT[0])
	
	## Raw VCF files from ensembl site might have NO format line
	## Only SNV files have AD in the format
	## Only INDEL files have a TAR in the format
	if record.FORMAT==None:
	    refcountt="0"
	    altcountt="0"
	    refcountn="0"
	    altcountn="0"
	elif "AD" in record.FORMAT:
	    refcountt=str(record.samples[1]['AD'][0])
	    altcountt=str(record.samples[1]['AD'][1])
	    refcountn=str(record.samples[0]['AD'][0])
	    altcountn=str(record.samples[0]['AD'][1])
	elif "TAR" in record.FORMAT:
	    refcountt=str(record.samples[1]['TAR'][0])
	    altcountt=str(record.samples[1]['TIR'][0])
	    refcountn=str(record.samples[0]['TAR'][0])
	    altcountn=str(record.samples[0]['TIR'][0])
	else:
	    print "Unknown file format; exiting"
	    sys.exit()

	## Calculate depths from the ALT/REF counts
	deptht=str(int(int(refcountt)+int(altcountt)))
	depthn=str(int(int(refcountn)+int(altcountn)))
	
	numtranscripts = len(record.INFO["CSQ"])
	transcripts=list()	
	for transcript in record.INFO["CSQ"]:
		pipesplit = pipesplitter(transcript)
		pipejoin = pipejoiner(pipesplit)
		transcripts.append(transcript)
	transcript_picker(transcripts)

logging.debug("All records parsed")


