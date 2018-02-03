echo 'STARTING JOG'     # prints to your output file

## This script is to process a bam file to determine coverage for defined list of regions (capture list). Coverage is analysed per locus and per interval specified in the capture list.


## Requirements: GATK v3.6 in $HOME/applications
## Usage example: ./GATK_coverage.sh  sample1_name

## Assign the command line arguments to correct variables, assuming paired ends data
samplename=$1

## Define some directories for GATK
GATKDIR=$HOME/applications

## Define some output files for GATK
reference=/data/BCI-BioInformatics/Jun/reference_hg38/hg38\.fa
capture_list=/data/BCI-BioInformatics/PC_ctDNA/WES_data/Agilent_Human_Exon_V6/S07604514_Covered_hg38_clean\.bed
outputbammarked=$samplename\.marked\.bam
coverage=$samplename\.bam\.coverage


## Determine coverage for defined list of regions
echo "####MESS Determine coverage for defined list of regions"
java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T DepthOfCoverage -R $reference -o $coverage -I $outputbammarked -ct 20 -ct 50 -ct 80 -ct 100 -ct 150 -ct 200 -L $capture_list

date
