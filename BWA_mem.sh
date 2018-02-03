echo 'STARTING JOG'     # prints to your output file

## This script is to run BWA 0.7.15-r1140 to align reads generated from the Illumina platform using 'mem' algorithm


## Requirements: BWA 0.7.15-r1140 installed in $HOME/bin 
## Usage example: ./BWA_mem.sh  sample1_name sample1_exome1_read_fq.1 sample1_exome1_read_fq.2

## Assign the command line arguments to correct variables, assuming paired ends data
#referencename=$1
#referenceindex=$2
samplename=$1
read1=$2
read2=$3

## Define some directories for BWA 
BWADIR=$HOME/bin

## Define some output files for BWA
read1name=$read1
read2name=$read2
outputsam=$samplename\.sam
referenceindex=/data/BCI-BioInformatics/Jun/reference_hg38/index_bwa_0\.7\.15/hg38\.fa

## step 1: BWA alignment
echo "####MESS Step 1:BWA mem alignment: read individual alignment"
## ATTENTION here, choose a unique ID name to integrate into the SAM output,e.g. Exome1
## this step is very important to not to cause any problems later on!!!
$BWADIR/bwa mem -t 4 -M -R "@RG\tID:$samplename\tLB:$samplename\tSM:$samplename\tPL:Illumina" $referenceindex $read1name $read2name > $outputsam
date
