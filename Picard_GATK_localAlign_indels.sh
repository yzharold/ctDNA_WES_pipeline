echo 'STARTING JOG'     # prints to your output file

## This script is to perform local alignment around indels


## Requirements: Picard v1.86 installed in $HOME/applications/picard-tools-1.119 and GATK v3.6 in $HOME/applications
## Usage example: ./Picard_GATK_localAlign_indels.sh  sample1_name

## Assign the command line arguments to correct variables, assuming paired ends data
samplename=$1

## Define some directories for Picard and GATK
PICARDDIR=$HOME/applications/picard-tools-2\.5\.0
GATKDIR=$HOME/applications

## Create temporary folder for picard (at this step many intermediate files are produced)
TEMPDIR=`pwd`/tmp\.$samplename

mkdir $TEMPDIR


## Define some output files for Picard and GATK
reference=/data/BCI-BioInformatics/Jun/reference_hg38/hg38\.fa
outputbammarked=$samplename\.marked\.bam
realignmentlist=$samplename\.bam\.list
realignmentbam=$samplename\.marked\.realigned\.bam
realignmentfixbam=$samplename\.marked\.realigned\.fixed\.bam


## step 4: local alignment around indels
echo "####MESS Step 4: local alignment around indels"
echo "####MESS Step 4: first create a table of possible indels"
java -Xmx4g -jar $GATKDIR/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -o $realignmentlist -I $outputbammarked
echo "####MESS Step 4: realign reads around those targets"
java -Xmx4g -Djava.io.tmpdir=/tmp -jar $GATKDIR/GenomeAnalysisTK.jar -I $outputbammarked -R $reference -T IndelRealigner -targetIntervals $realignmentlist -o $realignmentbam
echo "####MESS Step 4: fix paired end mate information using Picard"
java -Djava.io.tmpdir=$TEMPDIR -jar $PICARDDIR/picard.jar FixMateInformation INPUT=$realignmentbam OUTPUT=$realignmentfixbam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true TMP_DIR=$TEMPDIR


## Remove temporary folder for picard
rm -rf $TEMPDIR

date