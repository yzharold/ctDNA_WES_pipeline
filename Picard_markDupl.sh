echo 'STARTING JOG'     # prints to your output file

## This script is to mark PCR duplicates in BAM files using Picard


## Requirements: Picard v1.86 installed in $HOME/applications/picard-tools-1.119 
## Usage example: ./Picard_markDupl.sh  sample1_name

## Assign the command line arguments to correct variables, assuming paired ends data
samplename=$1

## Define some directories for Picard 
PICARDDIR=$HOME/applications/picard-tools-2\.5\.0

## Define some output files for Picard
outputsam=$samplename\.sam
outputbam=$samplename\.bam
outputbammarked=$samplename\.marked\.bam


## step 3: Marking PCR duplicates
echo "####MESS Step 3: Marking PCR duplicates using Picard"
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDDIR/picard.jar MarkDuplicates INPUT=$outputbam OUTPUT=$outputbammarked METRICS_FILE=$samplename\.DuplicationMetrics\.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
rm $outputsam
date
