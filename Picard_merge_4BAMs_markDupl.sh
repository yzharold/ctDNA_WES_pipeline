echo 'STARTING JOG'     # prints to your output file

## This script is to merge 4 BAM files (from same sample) and mark PCR duplicates using Picard


## Requirements: Picard v1.86 installed in $HOME/applications/picard-tools-1.119 
## Usage example: ./Picard_merge_4BAM_markDupl.sh  sample1_name bam_1 bam_2 bam_3 bam_4

## Assign the command line arguments to correct variables, assuming paired ends data
samplename=$1
bam_1=$2
bam_2=$3
bam_3=$4
bam_4=$5


## Define some directories for Picard 
PICARDDIR=$HOME/applications/picard-tools-2\.5\.0

## Define some output files for Picard
outputbammarked=$samplename\.merged\.marked\.bam


## step 3: Marking PCR duplicates
echo "####MESS Step 3: Marking PCR duplicates using Picard"
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $PICARDDIR/picard.jar MarkDuplicates INPUT=$bam_1 INPUT=$bam_2 INPUT=$bam_3 INPUT=$bam_4 OUTPUT=$outputbammarked METRICS_FILE=$samplename\.merged\.DuplicationMetrics\.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
date
