#!/usr/bin/env python2.7

## Secondary (temporary?) pipeline to perform variant calls and other 
## "real" analysis

## Import modules
import argparse
import logging

## Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("-n", type=str, help="Full path to input normal bam",required=False)
parser.add_argument("-t", type=str, help="Full path to input normal bam",required=False)
parser.add_argument("-p", type=str, help="Command string details here...",default="123456",required=False)
parser.add_argument("--debug", help="Turn debugging mode on",action="store_true")
args = parser.parse_args()

## Turn logging on if desired
if(args.debug):
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Debugging mode enabled")

## Declare some globals


# Steps to implement:
# ID test
# SNV calling with MuTect
# indel calling with Strelka
# Tx detection

def idTest():
    try:
	pass

    except:
	logging.exception("Error in ID test")

def main():
    try:
	logging.debug("Doing ID test for "+args.t+" vs "+args.n)
	idTest()
	pass
    except:
	logging.exception("Error in main")
	sys.exit()

## Execute 
if __name__ == '__main__':
	main()
