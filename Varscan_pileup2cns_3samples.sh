echo 'STARTING JOG'     # prints to your output file

##  Description: Call variants (SNPs and indels) using Varscan mpileup2cns for 3 samples. NOTE: the RG tags will be ignored and all reads in one BAM file will be treated as one sample.

## Requirements: Samtools and VarScan2 in $HOME/applications
## Usage example:  nohup ./Varscan_pileup2snp_pileup2indel.sh samples_name sample_1_bam sample_2_bam sample_3_bam  > sample_name_Varscan_somatic.log &


## Assign the command line arguments to correct variables
samplename=$1
input1bam=$2
input2bam=$3
input3bam=$4


## Define directories for Samtools
SAMTOOLSDIR=$HOME/applications/samtools-1.3.1
VARSCAN2DIR=$HOME/applications


## Define input and output files for Samtools and Varscan
reference=/data/BCI-BioInformatics/Jun/reference_hg38/hg38\.fa
vcfcns=$samplename\.cns.vcf


echo "Start samtools mpileup and varscan mpileup2cns " $1
$SAMTOOLSDIR/samtools mpileup -B -q 20 -Q 15 -R -f $reference $input1bam $input2bam $input3bam  | java -jar $VARSCAN2DIR/VarScan.v2.4.3.jar mpileup2cns --min-coverage 10 --min-reads2 1 --min-avg-qual 15 --min-var-freq 0.001 --p-value 0.99 --strand-filter 0 --output-vcf 1 --variants 1 > $vcfcns
echo "End samtools mpileup and varscan mpileup2cns " $1

date
