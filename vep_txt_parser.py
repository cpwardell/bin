#!/usr/bin/env python2.7

## Purpose; we want to ingest the text output of the VEP and output a file with a SINGLE
## row per variant.  We select the canonical transcript(s) first and then output the one
## with the worst outcome, according to ENSEMBL.  The process is repeated for variants with
## no canonical transcript

## Import modules
import logging
import argparse
import csv

## Gather command line args
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Turn logging on",action="store_true")
parser.add_argument("-i", type=str, help="full path to input txt",required=True)
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
	logging.basicConfig(level=logging.DEBUG)
	logging.debug("Debugging mode enabled")

### Define functions #################################################################

## Decide to append annotation to canonical or noncanonical list 
def appendmethod(row):
    if("CANONICAL=YES" in row[-1]):
	canonicals.append(row)
    else:
	noncanonicals.append(row)

## Decide between printing an annotation from the canonical or noncanonical list
## If the canonical list is empty, the noncanonical list is chosen
## Note that this method resets the canonical and noncanonical lists
def decider():
    global canonicals,noncanonicals
    ## Chose canonicals or noncanonicals
    if(len(canonicals)>=1):
	chosen=canonicals
    else:
	chosen=noncanonicals

    ## If there's only one option, print it
    ## Otherwise, decide which option is best by using the ordered ensembl annotation
    ## Stored here as the consequences dict
    if(len(chosen)==1):
	print "\t".join(chosen[0]) ## THIS IS WHERE THE ANNOTATION IS PRINTED 1 OF 2
	pass
    else:
	selector(chosen) 

    ## Reset the lists of objects
    canonicals=list()
    noncanonicals=list()

## Compares each annotation to the ensembl VEP ordering.  Prints highest-scoring (most 
## deleterious) annotation
def selector(chosen):
    scores=dict()
    for row in chosen:
	for consequence in consequences.keys():
	    #print row[6]
	    if(consequence in row[6]):
		scores[consequences[consequence]]=row
    #print("SCORE LEN IS: "+str(len(scores)))
    topscore = sorted(scores.keys())[0]
    print "\t".join(scores[topscore]) ## THIS IS WHERE THE ANNOTATION IS PRINTED 2 OF 2



### END OF FUNCTIONS, BEGIN MAIN SCRIPT ###########################################

## Open txt file for reading and create an empty dictionary to store all annotations
logging.debug("Opening text file")
rows=dict()
## Allow header to pass through and print it unmolested.  Keep track of 
## number of header rows so our iterator later starts at zero
headercount=0
for idx,row in enumerate(csv.reader(open(args.i),delimiter="\t")):
    if("#" in row[0]):
	print "\t".join(row)
	headercount+=1
    else:
	rows[idx-headercount]=row
logging.debug("Text file reading complete")
logging.debug(str(idx)+" rows in input file")

## Define some variables, including a dictionary used to score which effect is worst
consequences={"transcript_ablation":1,
       "splice_donor_variant":2,
       "splice_acceptor_variant":3,
       "stop_gained":4,
       "frameshift_variant":5,
       "stop_lost":6,
       "initiator_codon_variant":7,
       "transcript_amplification":8,
       "inframe_insertion":9,
       "inframe_deletion":10,
       "missense_variant":11,
       "splice_region_variant":12,
       "incomplete_terminal_codon_variant":13,
       "stop_retained_variant":14,
       "synonymous_variant":15,
       "coding_sequence_variant":16,
       "mature_miRNA_variant":17,
       "5_prime_UTR_variant":18,
       "3_prime_UTR_variant":19,
       "non_coding_transcript_exon_variant":20,
       "non_coding_exon_variant":20,
       "intron_variant":21,
       "NMD_transcript_variant":22,
       "non_coding_transcript_variant":23,
       "nc_transcript_variant":23,
       "upstream_gene_variant":24,
       "downstream_gene_variant":25,
       "TFBS_ablation":26,
       "TFBS_amplification":27,
       "TF_binding_site_variant":28,
       "regulatory_region_ablation":29,
       "regulatory_region_amplification":30,
       "regulatory_region_variant":31,
       "feature_elongation":32,
       "feature_truncation":33,
       "intergenic_variant":34}

thatid=None
canonicals=list()
noncanonicals=list()
## Now we have the results in a dict, loop through them and print one annotation per variant
logging.debug("Parsing "+str(len(rows.keys()))+" rows")
for key in rows.keys():
    row=rows[key]

    ## Assign id of this row
    thisid=row[0]

    ## If it's not the first row and the id is a new one, decide which row to print
    if(key!=0 and thisid != thatid):
	logging.debug("Deciding on row "+str(key)+" and ID "+str(thatid))
	decider()

    ## Append annotation to appropriate list
    appendmethod(row)
    
    ## If final row, execute the test
    if(key==len(rows)-1):
    	decider()

    ## Assign previous row id
    thatid=row[0]

logging.debug("All records parsed")


