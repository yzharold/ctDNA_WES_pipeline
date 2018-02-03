echo 'STARTING JOG'     # prints to your output file

## This script is to run convert SAM to BAM files using Picard


## Requirements: Picard v1.86 installed in $HOME/applications/picard-tools-1.119 
## Usage example: ./Picard_SAM2BAM.sh  sample1_name

## Assign the command line arguments to correct variables, assuming paired ends data
samplename=$1

## Define some directories for Picard 
PICARDDIR=$HOME/applications/picard-tools-2\.5\.0

## Define some output files for Picard
outputsam=$samplename\.sam
outputbam=$samplename\.bam


## step 2: SAM to BAM conversion using Picard
echo "####MESS Step 2: SAM to BAM conversion using Picard"
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDDIR/picard.jar SortSam SO=coordinate INPUT=$outputsam OUTPUT=$outputbam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
date
