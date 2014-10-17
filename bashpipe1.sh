#!/bin/bash

## Set up variables 
FWD=$1
REV=$2
RG=$3
FIRST=$4 # Forward fastq file
SECOND=$5 # Reverse fastq file

## Create a unique string to append to things:
USTRING=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 6 | head -n 1)


## Parse file names to produce output file names
FIRSTBASE=$(basename $FIRST .bz2)
SECONDBASE=$(basename $SECOND .bz2)
## VARIABLES
FORWARD=../trimmeddata/$FIRSTBASE.trimmed.gz
REVERSE=../trimmeddata/$SECONDBASE.trimmed.gz
BWA="/home/chris_w/apps/bwa-0.7.10/bwa"
INDEX="/home/chris_w/resources/b37/human_g1k_v37.fasta.gz"

INPUT=../align/$RG.sam
BASEINPUT=$(basename $INPUT .sam)


## Binary locations
PYTHON="/home/chris_w/apps/Python-2.7.8/python" # Python location
CUTADAPT="/home/chris_w/apps/cutadapt-1.5/bin/cutadapt" # cutadapt location
JAVA="/home/chris_w/apps/jre1.7.0_67/bin/java"
PICARD="/home/chris_w/apps/picard-tools-1.119"


##################
#### Cutadapt ####
##################
#FWD="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" # forward adapter
#REV="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" # reverse adapter

## NEXTERA TRANSPOSASE SEQUENCES:
#FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
#REV="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"


#mkdir -p trimmeddata
#cd trimmeddata
#cat >$PWD/cutadapt$RG.sh <<EOL
#!/bin/bash

## BZ2 files don't work, so we create local uncompressed fastq files
#bzcat $FIRST > $FIRSTBASE
#bzcat $SECOND > $SECONDBASE

## Perform filtering and delete temporary files
#$PYTHON $CUTADAPT -q 10 -a $FWD --minimum-length 70 -o $FIRSTBASE.temp.fastq -p $SECONDBASE.temp.fastq $FIRSTBASE $SECONDBASE
#$PYTHON $CUTADAPT -q 10 -a $REV --minimum-length 70 -o $SECONDBASE.trimmed.gz -p $FIRSTBASE.trimmed.gz $SECONDBASE.temp.fastq $FIRSTBASE.temp.fastq
#rm $FIRSTBASE.temp.fastq $SECONDBASE.temp.fastq $SECONDBASE $FIRSTBASE

#EOL

## Submit to queue
#qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 cutadapt$RG.sh
#cd ..

######################
#### Cutadapt END ####
######################

#############
#### BWA ####
#############

mkdir -p align
cd align

cat >$PWD/bwa$RG$USTRING.sh <<EOL
#!/bin/bash

$BWA mem -R "@RG\tID:$RG\tSM:$RG\tPL:ILLUMINA" $INDEX $FORWARD $REVERSE > $RG$USTRING.sam

EOL

## Submit to queue
#qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 -hold_jid cutadapt$RG.sh bwa$RG_$USTRING.sh 
qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 bwa$RG$USTRING.sh 
cd ..

#################
#### BWA END ####
#################

########################
#### Sort and dedup ####
########################

#mkdir -p sort
#cd sort
#cat >$PWD/sort$RG.sh <<EOL
#!/bin/bash

## Unset default Java options
#export JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -Xmx4G -Djava=/home/chris_w/tmp/ "

## Perform sort
#$JAVA -jar $PICARD/SortSam.jar INPUT=$INPUT OUTPUT=$BASEINPUT.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true

## Mark duplicates
#$JAVA -jar $PICARD/MarkDuplicates.jar INPUT=$BASEINPUT.sorted.bam OUTPUT=$BASEINPUT.dedup.bam CREATE_INDEX=true METRICS_FILE=$BASEINPUT.metrics REMOVE_DUPLICATES=false

## Tidy up
#rm $BASEINPUT.sorted.ba*

#EOL
 
## Submit to queue
#qsub -cwd -S /bin/bash -l s_vmem=8G -l mem_req=8 -hold_jid bwa$RG.sh sort$RG.sh 
#cd ..
############################
#### Sort and dedup END ####
############################
