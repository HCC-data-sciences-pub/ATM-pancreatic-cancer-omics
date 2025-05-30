
echo "hostname = $HOSTNAME"
hostname=htc ## stats or hpc; if hpc for some reason samtools cannot be found; only use stats even for hpc jobs

now=$(date +"%m-%d-%Y_%H:%M:%S")
	
## --------------------------------------
## command line options 
## --------------------------------------

project=ATM_PAAD ## MSLNpilot
sample=$1 ## MTS-1
genome=$2 ## grch38 grch37 hg19
part=$3 ## I or II
aligner=$4 ## bwamem
caller=$5 ## mutect2
tumor=$6
normal=$7

echo "Sample = $sample; part = $part"
if [ -z $sample ] || [ -z $genome ] || [ -z $aligner ]; then echo -e "Input missing. \n$0 <sample> <genome>\n"; exit; fi 

## --------------------------------------
## Tools
## --------------------------------------

## stats
if [ $hostname == 'stats' ]; then 

	## 08/15/2019: switch to gatk/4.1.3.0, R/3.5.3
	modules=' bam-readcount/0.7.4 bcftools/1.2 bedops/2.4.12 bedtools/2.23.0 bowtie2/2.2.6 bwa/0.7.10 cufflinks/2.2.1 fastqc/0.11.5 freebayes/0.9.20 gatk/4.1.3.0 igvtools/2.3.32 java/jdk1.8.0_45 kallisto/0.43.1 novocraft/3.02.08-free picard/2.9.3 R/3.5.3 rsem/1.2.19 rseqc/2.6 sambamba/0.6.3 samtools/1.2 SeqPrep/1.1 trimmomatic/0.36 UCSCtools python/2.7.8-bio pigz/2.3.1 salmon/0.8.2 star/2.5.3a subread/1.5.2 tabix/1.2.1 vcflib vcfsorter vcftools/0.1.12b '
	modules=' '
fi 

## gardner 
## 07/15/2017: do NOT use star/2.5.2a - parameter change between star/2.5.2a and 2.5.3a
## erre msg: EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "genomeFileSizes" in input "genomeParameters.txt"; SOLUTION: use correct parameter name (check the manual)
if [ $hostname == 'gardner' ]; then 
	## 08/11/2019: switch to gatk/4.1.3.0
	## 8/21/2019: switch to R/3.5.0
	modules=' java-jdk/1.8.0_92 gcc/6.2.0 zlib/1.2.8 bam-readcount/0.8.0 picard/2.17.3 bwa/0.7.17 fastqc/0.11.5 bedtools/2.26.0 mosdepth/0.2.6 bedops/2.4.26 bzip2/1.0.6 xz/5.2.2 perl/5.24.0 R/3.5.0 sambamba/0.6.7 UCSCtools pigz/2.3.4 tabix/1.2.1 trimmomatic/0.36 flash/1.2.11 gatk/4.1.3.0 abra/2.12 igvtools/2.3.91 annovar/20180416 lancet/1.0.6'

	## in order to run gatk4-scnv workflow:
	## pre-install all R packages in R 3.2.5 or higher
	## https://gatkforums.broadinstitute.org/gatk/discussion/11682
	## The plotting tools require particular R components. Options are to install these or to use the broadinstitute/gatk Docker. In particular, to match versions, use the broadinstitute/gatk:4.0.1.1 version.
    ## Install R v3.2.5 or above from https://www.r-project.org/, then install the components using the install_R_packages.R script with Rscript install_R_packages.R.
    ## https://github.com/broadinstitute/gatk/blob/4.0.1.1/scripts/docker/gatkbase/install_R_packages.R

    ## 08/21/2019 (temperary solution ....)
	export PATH=$PATH:/group/bioinformatics/software/GATK/4.1.3.0:/group/bioinformatics/software/Control-FREEC/11.5/src

fi 

if [ $hostname == 'htc' ]; then 
	modules='gcc/8.2.0 bedtools/2.29.0 samtools/1.9 r/4.0.0'
	export PATH=$PATH:Software/freec/11.6/src
fi 





## --------------------------------------
## Settings
## --------------------------------------

threads=8

## --------------------------------------
## Paths and directories
## --------------------------------------

projPath=Projects/HCC-CBS-011-IPM-FCouch
outputPath="$projPath/results/exome/${project}_samples_$genome/$sample"
dirlist=" scnv_calls scnv_calls/gatkscnv"

## --------------------------------------
## Commands
## --------------------------------------

for module in $modules; do module load $module; done 
module list 

echo "## --------------------------------------"

cd $projPath
pwd=$PWD
echo "Current Path = $pwd"

## http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
# /group/bioinformatics/software/UCSCtools/faSplit byname /group/referenceFiles/Homo_sapiens/GDC/GRCh38.d1.vd1/GRCh38.d1.vd1.fa chromosomes/

# gunzip -c /group/referenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0/dbsnp_146.hg38.vcf.gz | awk -F"\t" '$1~/^chr/{print $1"\t"$2-1"\t"$2}' > /group/referenceFiles/Homo_sapiens/GDC/broad-references/hg38/v0/dbsnp_146.hg38.bed

echo -e "Running Control-FREEC: tumor = $tumor; normal = $normal"
config=$projPath/scripts/batches/freec.$tumor.cfg

if [ ! -d $outputPath/scnv_calls/freec ]; then mkdir -p $outputPath/scnv_calls/freec; fi 

# module load gcc/6.2.0 samtools/1.6.0 bedtools/2.27.1 R/3.5.0
# module load dmd/2.072.1 sambamba/0.6.5 ## this will inactivate all other modules...
which samtools
which bedtools
which R

# freec -conf $config

echo -e "Running Control-FREEC: calculating p-values and making plots"
cnv=$outputPath/scnv_calls/freec/$tumor.bwamem.merged.rmdup.bam_CNVs
ratio=$outputPath/scnv_calls/freec/$tumor.bwamem.merged.rmdup.bam_ratio.txt
baf=$outputPath/scnv_calls/freec/$tumor.bwamem.merged.rmdup.bam_BAF.txt
info=$outputPath/scnv_calls/freec/$tumor.bwamem.merged.rmdup.bam_info.txt
ploidy=`awk -F"\t" '$1=="Output_Ploidy" {print $2}' $info`

cat Software/freec/11.6//scripts/assess_significance.R | R --slave --args $cnv $ratio

echo -e "ploidy = $ploidy"
cat Software/freec/11.6/scripts/makeGraph.R | R --slave --args $ploidy $ratio $baf

echo "## --------------------------------------"

echo "[ " `date` " ] Program finished!"




