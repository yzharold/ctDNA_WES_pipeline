echo 'STARTING JOG'     # prints to your output file

## This script is to perform base quality score recalibration


## Requirements: GATK v3.6 in $HOME/applications
## Usage example: ./GATK_baseQrecalib.sh  sample1_name

## Assign the command line arguments to correct variables, assuming paired ends data
samplename=$1

## Define some directories for Picard and GATK
GATKDIR=$HOME/applications

## Define some output files for Picard and GATK
reference=/data/BCI-BioInformatics/Jun/reference_hg38/hg38\.fa
dbSNP=/data/BCI-BioInformatics/Jun/dbSNP/common_all_20161122_chr_sorted\.vcf\.gz
baserecaldata=$samplename\.recal_data\.grp
realignmentfixbam=$samplename\.marked\.realigned\.fixed\.bam
recalioutbam=$samplename\.recalib\.bam


## step 5: base quality score recalibration
echo "####MESS Step 5: base quality score recalibration"
java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T BaseRecalibrator -I $realignmentfixbam -R $reference -knownSites $dbSNP -o $baserecaldata
echo "####MESS Step 5: print recalibrated reads into BAM"
java -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -R $reference -I $realignmentfixbam -BQSR $baserecaldata -o $recalioutbam
date
date
