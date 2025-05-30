
echo "hostname = $HOSTNAME"

now=$(date +"%m-%d-%Y_%H:%M:%S")
alignRawReads=1 ## whether to directly align raw reads or trimmed reads
# strand=reverse ## reverse, forward, unstranded
strand=reverse ## reverse, forward, unstranded

## --------------------------------------
## command line options 
## --------------------------------------

project=ATM_PAAD
sample=$1 ## MTS-1
genome=$2 ## grch38 grch37 hg19
aligner=$3 ## star salmon kallisto star2pass

echo "Sample = $sample"
if [ -z $sample ] || [ -z $genome ] || [ -z $aligner ]; then echo -e "Input missing. \n$0 <sample> <genome>\n"; exit; fi 

# genome=$1 ## USA300_FPR3757
# sample=$2 ## mrsaA_F1
# run=$3 ## single multi
# aligner=$4 ## bwamem novo 
# caller=$5 ## mpileup freebayes gatkhc 

## --------------------------------------
## Tools
## --------------------------------------

# ## tarbell
# modules=' bam-readcount/0.7.4 bcftools/1.2 bedops/2.4.12 bedtools/2.23.0 bowtie2/2.2.6 bwa/0.7.10 cufflinks/2.2.1 fastqc/0.11.5 freebayes/0.9.20 gatk/3.4.0 igvtools/2.3.32 java/jdk1.8.0_45 kallisto/0.43.1 novocraft/3.02.08-free picard/2.9.3 R/3.2.1 rsem/1.2.19 rseqc/2.6 sambamba/0.6.3 samtools/1.2 SeqPrep/1.1 trimmomatic/0.36 UCSCtools python/2.7.8-bio pigz/2.3.1 salmon/0.8.2 star/2.5.3a subread/1.5.2 tabix/1.2.1 vcflib vcfsorter vcftools/0.1.12b '

## gardner 
## 07/15/2017: do NOT use star/2.5.2a - parameter change between star/2.5.2a and 2.5.3a
## erre msg: EXITING: FATAL INPUT ERROR: unrecoginzed parameter name "genomeFileSizes" in input "genomeParameters.txt"; SOLUTION: use correct parameter name (check the manual)
## 04/28/2018: switched to kallisto/0.44.0 star/2.6.0a
## 05/03/2018: changed from python/2.7.8-bio to python/2.7.13 (on gardner)
## 05/04/2018: changed from gatk/3.4.0 to gatk/4.0.4.0
## 05/05/2018: added gcc/6.2.0
## 05/15/2018: switched to salmon/0.9.1
# modules=' gcc/6.2.0 java-jdk/1.8.0_92 bam-readcount/0.7.4 bcftools/1.2 bedops/2.4.12 bedtools/2.23.0 bowtie2/2.2.6 bwa/0.7.10 cufflinks/2.2.1 fastqc/0.11.5 freebayes/0.9.20 gatk/4.0.4.0 igvtools/2.3.32 kallisto/0.44.0 novocraft/3.02.08-free picard/2.8.1 R/3.2.1 rsem/1.2.19 rseqc/2.6 sambamba/0.6.3 samtools/1.2 SeqPrep/1.1 trimmomatic/0.36 UCSCtools python/2.7.13 pigz/2.3.1 salmon/0.9.1 star/2.6.0a subread/1.5.0-p2 tabix/1.2.1 vcflib vcfsorter vcftools/0.1.12b '
## switch to pitt crc htc
modules=' gcc/8.2.0 bcftools/1.9-40 bedops/2.4.35 bedtools/2.29.0 bowtie2/2.3.4.2 bwa/0.7.17 cufflinks/2.2.1 fastqc/0.11.7 freebayes/1.2.0 gatk/4.1.2.0 r/3.6.0 rseqc/2.6.6 sambamba/0.6.8 samtools/1.9 subread/1.6.2 '

# checkQC=/group/bioinformatics/Pipelines/Development/CRI-GMD-Pipeline/gmd-EC2-v1.0.0/BDS-ExScaliburGMD-test/util/Check_QC.pl
# pear=/group/bioinformatics/software/pear/0.9.6/pear-0.9.6-bin-64
# flash=/group/bioinformatics/software/FLASH/1.2.11/flash

export PATH=Software/UCSCtools:Software/star/2.7.3a/bin/Linux_x86_64_static:Software/kallisto/0.46.1:Software/salmon/1.0.0/bin:$PATH
PICARD=/ihome/crc/install/picard/2.18.12/picard.jar
TRIMMOMATIC=/ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar

## fastqc returned error:
## /usr/bin/perl: symbol lookup error: /ihome/crc/install/gcc-8.2.0/perl/5.28.0/lib/5.28.0/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined symbol: Perl_xs_handshake
## solution: https://github.com/bioconda/bioconda-recipes/issues/4390
export PERL5LIB=""

## not loaded ...
# bam-readcount
# igvtools
# SeqPrep
# pigz
# tabix
# vcfsorter
# vcflib/1.0.0 vcftools/0.1.16 rsem/1.3.1 

## --------------------------------------
## Settings
## --------------------------------------

threads=8 ## 12 for star, change mem to 36Gb; 4 for other parts, change mem to 10Gb 
## "/ihome/crc/install/trimmomatic/Trimmomatic-0.38/adapters/"
adpFasta="TruSeq3-PE.fa" ## TruSeq3-SE.fa NexteraPE-PE.fa
qcMatrics="QUAL:30:50,NFD:0:10,GC:40:60,N:0:2"
trimEncodingParam=" -phred33"
clipEncodingParam=""
pearEncodingParam=" -b 33"
flashEncodingParam=" -p 33"
if [ $adpFasta == 'TruSeq3-PE.fa' ]; then 
	adp1="TACACTCTTTCCCTACACGACGCTCTTCCGATCT" ## TruSeq3-PE.fa R1 
	adp2="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" ## TruSeq3-PE.fa R2
fi 
if [ $adpFasta == 'TruSeq3-SE.fa' ]; then 
	adp1="AGATCcheckQCGGAAGAGCACACGTCTGAACTCCAGTCAC"  
	adp2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" 
fi 
if [ $adpFasta == 'NexteraPE-PE.fa' ]; then 
	adp1="AGATGTGTATAAGAGACAG" 
	adp2="AGATGTGTATAAGAGACAG" 
fi 

picardStrand="NONE"
featureCountsStrand=0
if [ $strand == "reverse" ]; then 
	picardStrand="SECOND_READ_TRANSCRIPTION_STRAND"
	featureCountsStrand=2
fi 
if [ $strand == "forward" ]; then 
	picardStrand="FIRST_READ_TRANSCRIPTION_STRAND"
	featureCountsStrand=1
fi 

# insertSizeMean=200
# insertSizeStd=100
## 11/27/2017
insertSizeMean=250
insertSizeStd=30 

## --------------------------------------
## Paths and directories
## --------------------------------------

projPath=
outputPath="results/rnaseq/${project}_samples_${genome}_${strand}/$sample"
# outputPath="results/rnaseq/${project}_samples_$genome/$sample"
tmpDir=$outputPath/tmp
dirlist=' qc_reports clean_reads alignment alignment_stats alignment_stats/picard alignment_stats/rseqc alignment_stats/bedtools read_counts read_counts/raw read_counts/salmon read_counts/kallisto tmp '

# if [ $aligner == 'star' ] || [ $aligner == 'star2pass' ]; then outputPath="results/rnaseq/${project}_samples_$genome/$sample"; fi 

## --------------------------------------
## Files
## --------------------------------------

refGenome=
refDict=
refGTF=
starIndex=
# rsemIndex=Homo_sapiens/RSEMgenome/GRCh37.p13_Gencode19/bt2/human_ref
chromSize=
# metadata=$projPath/sampleinfo/$project.rnaseq.metadata.txt
metadata=$projPath/sampleinfo/$project.rnaseq.metadata.20200520.txt

## --------------------------------------
## Commands
## --------------------------------------

for module in $modules; do module load $module; done 
module list 

echo "## --------------------------------------"

cd $projPath
pwd=$PWD
echo "Current Path = $pwd"

for dir in $dirlist; do if [ ! -d $outputPath/$dir ]; then mkdir -p $outputPath/$dir; fi; done 

echo "## --------------------------------------"

echo "Aligner = $aligner"

echo "## --------------------------------------"

## get attributes from the metadata table
rgTotal=`awk -F"\t" -v s=$sample '$1==s' $metadata | cut -f 1- | wc -l` 
echo "readgroup total = $rgTotal"

## =========================== PART I ================================ ## 

## loop through readgroups per sample 
bamString=""
readgroups=""
r1files=""
r2files=""
for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 

	echo "## --------------------------------------"
	echo "readgroup line = $rgline"

	values=( `awk -F"\t" -v s=$sample '$1==s' $metadata | cut -f 1- | head -$rgline | tail -1 ` )
	header=`head -1 $metadata`

	# echo $header
	smIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Sample" {print $1-1}'`
	lbIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Library" {print $1-1}'`
	rgIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="ReadGroup" {print $1-1}'`
	plIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Platform" {print $1-1}'`
	cnIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="SequencingCenter" {print $1-1}'`
	dtIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Date" {print $1-1}'`
	lnIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Lane" {print $1-1}'`
	puIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Unit" {print $1-1}'`
	flIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Flavor" {print $1-1}'`
	ecIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Encoding" {print $1-1}'`
	locIndex=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Location" {print $1-1}'`
	r1Index=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Seqfile1" {print $1-1}'`
	r2Index=`echo $header | sed 's/ /\n/g' | cat -n | awk '$2=="Seqfile2" {print $1-1}'`

	## note shell array index starts with offset 0
	sample_id=${values[smIndex]}
	library=${values[lbIndex]}
	readgroup=${values[rgIndex]}
	platform=${values[plIndex]}
	center=${values[cnIndex]}
	date=${values[dtIndex]}
	lane=${values[lnIndex]}
	unit=${values[puIndex]}
	flavor=${values[flIndex]}
	encoding=${values[ecIndex]}
	loc=${values[locIndex]}
	r1=${values[r1Index]}
	# echo "r2Index = $r2Index"
	if [ ! -z $r2Index ] && [ $r2Index != "" ]; then r2=${values[r2Index]}; fi 
	# echo "r2 = $r2"

	## this should be set for each readgroup NOT over all readgroups
	libraryType='SE'
	if [ ! -z $r2 ] && [ $r2 != "" ];then libraryType='PE'; fi 
	if [ -z $libraryType ] || [ $libraryType == "" ]; then echo "libraryType is not valid. Program terminates."; exit; fi 

	## set adapter sequence accordingly
	if [ $adpFasta == "NexteraPE-PE.fa" ]; then 
		if [ $libraryType != 'PE' ]; then 
			echo "adpFasta is wrong!"
			echo $sample_id $library $readgroup $platform $center $date $lane $unit $flavor $encoding $loc $r1 $r2 $libraryType $adpFasta
			exit
		fi 
	else 
		if [ $libraryType == 'SE' ]; then adpFasta="TruSeq3-SE.fa";fi 
		if [ $libraryType == 'PE' ]; then adpFasta="TruSeq3-PE.fa";fi 
	fi 

	if [ $adpFasta == 'TruSeq3-PE.fa' ]; then 
		adp1="TACACTCTTTCCCTACACGACGCTCTTCCGATCT" 
		adp2="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" 
	fi 
	if [ $adpFasta == 'TruSeq3-SE.fa' ]; then 
		adp1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" 
		adp2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" 
	fi 
	if [ $adpFasta == 'NexteraPE-PE.fa' ]; then 
		adp1="AGATGTGTATAAGAGACAG" 
		adp2="AGATGTGTATAAGAGACAG" 
	fi 

	echo $sample_id $library $readgroup $platform $center $date $lane $unit $flavor $encoding $loc $r1 $r2 $libraryType $adpFasta

	## modify star index for different seq length 
	if [ $genome == 'grch37' ]; then 
		if [ $flavor == "50x1" ] || [ $flavor == "49x1" ] || [ $flavor == "50x2" ] || [ $flavor == "49x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/GRCh37.p13_Gencode19_maskPAR_50bp

			## 04/28/2018: old index is not compatible with star2pass....
			## rebuilt the index using new star versions
			if [ $aligner == 'star2pass' ]; then 
				starIndex=Homo_sapiens/STARgenome/GRCh37.p13_Gencode19_maskPAR_50bp_2.6.0a
			fi 
		fi 
		if [ $flavor == "75x2" ] || [ $flavor == "76x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/GRCh37.p13_Gencode19_maskPAR_75bp
			
			## 04/28/2018: old index is not compatible with star2pass....
			## rebuilt the index using new star versions
			if [ $aligner == 'star2pass' ]; then 
				starIndex=Homo_sapiens/STARgenome/GRCh37.p13_Gencode19_maskPAR_75bp_2.6.0a
			fi 
		fi 
		if [ $flavor == "100x2" ] ||  [ $flavor == "101x2" ] || [ $flavor == "95x2" ] || [ $flavor == "103x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/GRCh37.p13_Gencode19_maskPAR_100bp

			## 04/28/2018: old index is not compatible with star2pass....
			## rebuilt the index using new star versions
			if [ $aligner == 'star2pass' ]; then 
				starIndex=Homo_sapiens/STARgenome/GRCh37.p13_Gencode19_maskPAR_100bp_2.6.0a
			fi 
		fi 
		refGTF=$starIndex/gencode.v19.annotation.maskPAR.gtf
		refGenome=$starIndex/GRCh37.p13.genome.fa
		refBED=$starIndex/gencode.v19.annotation.maskPAR.bed
		refBED12=$starIndex/gencode.v19.annotation.maskPAR.bed12
		refFLAT=$starIndex/gencode.v19.annotation.maskPAR.refFlat.txt
		refRibosomeRNABED=$starIndex/gencode.v19.annotation.rRNA.bed
		refRibosomeRNAINTERVAL=$starIndex/gencode.v19.annotation.rRNA.interval_list
		chromSize=$starIndex/GRCh37.p13.genome.chrom.size
	fi 
	if [ $genome == 'grch38' ]; then 
		if [ $flavor == "50x1" ] || [ $flavor == "49x1" ] || [ $flavor == "50x2" ] || [ $flavor == "49x2" ]; then  
			starIndex=Homo_sapiens/STARgenome/GRCh38.primary_Gencode32_maskPAR_50bp
		fi
		if [ $flavor == "75x2" ] || [ $flavor == "76x2" ]; then
			starIndex=Homo_sapiens/STARgenome/GRCh38.primary_Gencode32_maskPAR_75bp
		fi 
		if [ $flavor == "100x2" ] ||  [ $flavor == "101x2" ] || [ $flavor == "95x2" ] || [ $flavor == "103x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/GRCh38.primary_Gencode32_maskPAR_100bp
		fi 
		if [ $flavor == "150x2" ] ||  [ $flavor == "151x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/GRCh38.primary_Gencode32_maskPAR_150bp
		fi 
		refGTF=$starIndex/gencode.v32.primary_assembly.annotation.maskPAR.gtf
		refGenome=$starIndex/GRCh38.primary_assembly.genome.fa
		refBED=$starIndex/gencode.v32.primary_assembly.annotation.maskPAR.bed
		refBED12=$starIndex/gencode.v32.primary_assembly.annotation.maskPAR.bed12
		refFLAT=$starIndex/gencode.v32.primary_assembly.annotation.maskPAR.refFlat.txt
		# refRibosomeRNABED=$starIndex/gencode.v32.primary_assembly.annotation.maskPAR.rRNA.bed
		# refRibosomeRNAINTERVAL=$starIndex/gencode.v32.primary_assembly.annotation.maskPAR.rRNA.interval_list
		refRibosomeRNABED=$starIndex/GRCh38_rRNA.bed
		refRibosomeRNAINTERVAL=$starIndex/GRCh38_rRNA.interval_list
		chromSize=$starIndex/GRCh38.primary_assembly.genome.chrom.size
	fi 
	if [ $genome == 'hg19' ]; then 
		if [ $flavor == "75x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/hg19_UCSC_75bp
		fi 
		if [ $flavor == "100x2" ] ||  [ $flavor == "101x2" ] || [ $flavor == "95x2" ]; then 
			starIndex=Homo_sapiens/STARgenome/hg19_UCSC_100bp
		fi 
		refGTF=Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
		refGenome=Homo_sapiens/UCSC/hg19/hg19.GATKbundle.2.8/hg19/ucsc.hg19.fasta
		refBED=Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.bed
		# refBED12=$starIndex/gencode.v19.annotation.bed12
		refFLAT=Homo_sapiens/UCSC/hg19/Annotation/Genes/refFlat.txt
		refRibosomeRNABED=Homo_sapiens/UCSC/hg19/Annotation/Genes/rRNA.bed
		refRibosomeRNAINTERVAL=Homo_sapiens/UCSC/hg19/Annotation/Genes/rRNA.bed.locations.txt
		chromSize=Homo_sapiens/UCSC/hg19/hg19.GATKbundle.2.8/hg19/ucsc.hg19.chrom.size
	fi 

	echo "star index = $starIndex"

	## sample sanity check 
	if [ $sample != $sample_id ]; then 
		echo "Sample is not consistent with metadata table. $sample vs $sample_id"
		exit
	fi 

	## check fastqc format (default assuming it is 33)
	if [ $encoding == "64" ] || [ $encoding == "phred64" ] || [ $encoding == "Phred64" ]; then
		echo "Fastq Encoding = $encoding"
		trimEncodingParam=" -phred64"
		clipEncodingParam=" -6"
		pearEncodingParam=" -b 64"
		flashEncodingParam=" -p 64"

		# ## convert 64 to 33
		# echo "[ " `date` " ] Convert FastQ Phred33 to 64: Trimmomatic"
		# for in in $loc/$r1 $loc/$r2; do 
		# 	echo $in 
		# 	out=$in.phred33.fq.gz
		# 	java -Xmx6G -jar ${TRIMMOMATIC} SE -threads $threads $trimEncodingParam $in $out TOPHRED33

		# 	# mv $in $in.bak
		# 	# mv $out $in
		# done 
	fi 

	## readgroup header string (for aligners)
	RGstring="@RG\tID:$readgroup\tPL:$platform\tLB:$library\tSM:$sample\tCN:$center\tDT:$date\tPU:$unit"
	echo "RGstring = $RGstring"

	RGstringSTAR="ID:$readgroup PL:$platform LB:$library SM:$sample CN:$center DT:$date PU:$unit"

	if [ $rgline -eq 1 ]; then RGstringSTARstring="$RGstringSTAR"; fi
	if [ $rgline -gt 1 ]; then RGstringSTARstring="$RGstringSTARstring , $RGstringSTAR"; fi 
	

	# ## =========================== PART II ================================ ## 

	# if [ $aligner == 'star' ]; then 

	# 	echo "## --------------------------------------"

	# 	## fastqc (raw reads)
	# 	echo "[ " `date` " ] Collect raw reads qc metrics: fastqc"
	# 	inString="$loc/$r1"
	# 	if [ $libraryType == 'PE' ]; then inString="$inString $loc/$r2"; fi 
	# 	echo "$inString"
	# 	# fastqc --extract -o $outputPath/qc_reports -t $threads --nogroup $inString

	# 	# echo "## --------------------------------------"

	# 	# ## check qc results 
	# 	# echo "[ " `date` " ] Check raw reads qc metrics: checkQC.pl"
	# 	# inString="$outputPath/qc_reports/$r1:$outputPath/qc_reports/$r2"
	# 	# inString=${inString//.fastq.gz/_fastqc.zip}
	# 	# $checkQC -f $inString -rg $readgroup -s $sample -o $outputPath -m $qcMatrics

	# 	echo "## --------------------------------------"

	# 	if [ $alignRawReads -eq 0 ]; then 
	# 		## trimmomatic: trim low quality bases, clip adapters
	# 		if [ $libraryType == 'SE' ]; then 
	# 			echo "[ " `date` " ] Trim low quality bases (and clip adapters): Trimmomatic"
	# 			in1=$loc/$r1
	# 			out1=$outputPath/clean_reads/$readgroup.R1.trim.fq.gz
	# 			java -Xmx6G -jar ${TRIMMOMATIC} SE -threads $threads $trimEncodingParam $in1 $out1 ILLUMINACLIP:${ADAPTERS}/$adpFasta:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

	# 			## join all SE reads
	# 			echo "Double check all adapters were removed ... [adapter = $adp1 $adp2]"
	# 			for file in $out1 ; do 
	# 				echo $file
	# 				zcat $file | grep $adp1 | wc -l
	# 				zcat $file | grep $adp2 | wc -l
	# 			done 
	# 		fi 
	# 		if [ $libraryType == 'PE' ]; then 
	# 			echo "[ " `date` " ] Trim low quality bases (and clip adapters): Trimmomatic"
	# 			in1=$loc/$r1
	# 			in2=$loc/$r2
	# 			out1=$outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz
	# 			out2=$outputPath/clean_reads/$readgroup.R1.trim.unpaired.fq.gz
	# 			out3=$outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz
	# 			out4=$outputPath/clean_reads/$readgroup.R2.trim.unpaired.fq.gz
	# 			out5=$outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz
	# 			java -Xmx6G -jar ${TRIMMOMATIC} PE -threads $threads $trimEncodingParam $in1 $in2 $out1 $out2 $out3 $out4 ILLUMINACLIP:${ADAPTERS}/$adpFasta:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

	# 			## join all SE reads
	# 			tmpFile=$outputPath/clean_reads/fq.tmp
	# 			if [ -e $tmpFile ]; then rm $tmpFile; fi 
	# 			for file in $out2 $out4; do echo $file; zcat $file >> $tmpFile; done 
	# 			gzip -c $outputPath/clean_reads/fq.tmp > $out5
	# 			rm $tmpFile
	# 			echo "Double check all adapters were removed ... [adapter = $adp1 $adp2]"
	# 			for file in $out1 $out2 $out3 $out4 $out5; do 
	# 				echo $file
	# 				zcat $file | grep $adp1 | wc -l
	# 				zcat $file | grep $adp2 | wc -l
	# 			done 
	# 		fi 

	# 		echo "## --------------------------------------"

	# 		## fastqc (clean reads )
	# 		echo "[ " `date` " ] Collect clean reads qc metrics: fastqc"
	# 		inString=""
	# 		for file in $outputPath/clean_reads/$readgroup.*.gz; do 
	# 			echo $file
	# 			inString="$inString $file"
	# 		done 
	# 		fastqc --extract -o $outputPath/qc_reports -t $threads --nogroup $inString
	# 	fi 

	# 	## =========================== PART III ================================ ## 

	# 	echo "## --------------------------------------"

	# 	## star: alignment; Note it uses 48 Gb mem with 12 threeads but FAST! Change the mem and ppn and threads at this point! also add qos=biocorex in PBS header
	# 	echo "[ " `date` " ] Map to reference genome: STAR"
	# 	if [ $alignRawReads -eq 1 ]; then 
	# 		if [ $libraryType == 'SE' ]; then 
	# 			in1=$loc/$r1
	# 			out1=$outputPath/alignment/$readgroup.$aligner.

	# 			echo "Align SE reads .... "
	# 			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
	# 		fi 
	# 		if [ $libraryType == 'PE' ]; then 
	# 			in1=$loc/$r1
	# 			in2=$loc/$r2
	# 			out1=$outputPath/alignment/$readgroup.$aligner.

	# 			echo "Align PE reads .... "
	# 			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 $in2 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
	# 		fi 
	# 	fi 
	# 	if [ $alignRawReads -eq 0 ]; then 
	# 		if [ $libraryType == 'SE' ]; then 
	# 			in1=$outputPath/clean_reads/$readgroup.R1.trim.fq.gz
	# 			out1=$outputPath/alignment/$readgroup.$aligner.

	# 			echo "Align SE reads .... "
	# 			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
	# 		fi 
	# 		if [ $libraryType == 'PE' ]; then 
	# 			in1=$outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz
	# 			in2=$outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz
	# 			in3=$outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz
	# 			out1=$outputPath/alignment/$readgroup.$aligner.pe.
	# 			out2=$outputPath/alignment/$readgroup.$aligner.se.
	# 			out3=$outputPath/alignment/$readgroup.$aligner.Aligned.sortedByCoord.out.bam

	# 			echo "Align PE reads .... "
	# 			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 $in2 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
	# 			echo "Align SE reads .... "
	# 			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out2 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in3 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR

	# 			echo "Merge PE and SE read alignment ..."
	# 			sambamba merge -t $threads $out3 ${out1}Aligned.sortedByCoord.out.bam ${out2}Aligned.sortedByCoord.out.bam
	# 		fi 
	# 	fi

	# fi 

	echo "## --------------------------------------"

	## record which bam file output 
	bamString="$bamString $outputPath/alignment/$readgroup.$aligner.Aligned.sortedByCoord.out.bam"

	## record readgroups
	readgroups="$readgroups $readgroup"

	## record R1 fastq files
	r1files="$r1files $loc/$r1"
	r2files="$r2files $loc/$r2"
done ; 

echo "## --------------------------------------"

if [ $aligner == 'star' ]; then 

	## combined readgroup BAM files per sample 
	echo "[ " `date` " ] Merge readgroup alignment into sample level: sambamba"
	echo $bamString
	echo $readgroups
	out=$outputPath/alignment/$sample.$aligner.merged.bam
	if [ $rgTotal -eq 1 ]; then 
		echo "readgroup total = $rgTotal. No need to merge!"
		link=${bamString/*\/}
		if [ ! -e $out ]; then ln -s $link $out; fi 
		if [ ! -e $out.bai ]; then sambamba index -t $threads $out; fi 
	fi 
	if [ $rgTotal -gt 1 ]; then 
		echo "readgroup total = $rgTotal. Merge bamfiles!"
		sambamba merge -t $threads $out $bamString
	fi

	# echo "## --------------------------------------"

	# ## output MapQ distribution to make sure featureCounts behave as expected 
	# ## (eg. counts only primary & unique alignment)
	# echo "[ " `date` " ] Count MapQ distribution of the bam file ..."
	# sambamba view -t $threads $out | cut -f 5 | sort | uniq -c > $out.mapq_hist

	echo "## --------------------------------------"

	## freaturecounts: count raw reads  
	## featurecount works with BAM of both PE and SE reads (https://groups.google.com/forum/#!topic/subread/qImVtrKtZsw)
	echo "[ " `date` " ] Count raw reads: featureCounts"
	in=$outputPath/alignment/$sample.$aligner.merged.bam 
	out1=$outputPath/read_counts/raw/$sample.$aligner.merged.bam.primary.raw_counts.txt
	out2=$outputPath/read_counts/raw/$sample.$aligner.merged.bam.primary_uniq.raw_counts.txt
	out3=$outputPath/read_counts/raw/$sample.$aligner.merged.bam.mapq1.raw_counts.txt

	## Note that "Assigned + Unassigned_Ambiguity + Unassigned_NoFeatures == number of reads with MapQ 255"
	## STAR manual: uniquely mapped reads are assigned with mapq 255
	echo "Count primary alignments and uniquely mapped reads only ..."
	featureCounts -s $featureCountsStrand -a $refGTF -o $out2 -p -t exon -g gene_id -Q 255 -T $threads -J -G $refGenome --primary -C $in

	echo "Count primary alignments only ..."
	featureCounts -s $featureCountsStrand -a $refGTF -o $out1 -p -t exon -g gene_id -Q 0 -T $threads -J -G $refGenome --primary -C $in

	## Kyle's commands 
	echo "Count reads with mapping quality at least 1; no primary or unique filters ..."
	featureCounts -s $featureCountsStrand -a $refGTF -o $out3 -p -t exon -g gene_id -Q 1 -T $threads -G $refGenome $in

	## =========================== PART IV ================================ ## 

	echo "## --------------------------------------"

	## picard: collect rnaseq matrics
	echo "[ " `date` " ] Collect alignment metrics: picard CollectMultipleMetrics"

	echo "Collect multiple metrics ..." ## PROGRAM=RnaSeqMetrics will throw nullpoints in java! 
	## cannot process in mixed PE & SE bam files ... need to run on each
	if [ $alignRawReads -eq 1 ]; then 
		in=$outputPath/alignment/$sample.$aligner.merged.bam
		out=$outputPath/alignment_stats/picard/$sample.$aligner.merged.bam

		if [ $libraryType == 'SE' ]; then 
			java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
		fi 
		if [ $libraryType == 'PE' ]; then 
			java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
		fi 
	fi 

	if [ $alignRawReads -eq 0 ]; then 
		echo "readgroups = $readgroups"
		
		if [ $libraryType == 'SE' ]; then 
			for readgroup in $readgroups; do 
				echo "SE read bam ... $readgroup"
				in=$outputPath/alignment/$readgroup.$aligner.Aligned.sortedByCoord.out.bam
				out=$outputPath/alignment_stats/picard/$readgroup.$aligner.bam
				java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
			done 
		fi 
		if [ $libraryType == 'PE' ]; then 
			for readgroup in $readgroups; do 
				echo "PE read bam ... $readgroup"
				in=$outputPath/alignment/$readgroup.$aligner.pe.Aligned.sortedByCoord.out.bam
				out=$outputPath/alignment_stats/picard/$readgroup.$aligner.pe.bam
				java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000

				echo "SE read bam ... $readgroup"
				in=$outputPath/alignment/$readgroup.$aligner.se.Aligned.sortedByCoord.out.bam
				out=$outputPath/alignment_stats/picard/$readgroup.$aligner.se.bam
				java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
			done 
		fi 
	fi 

	echo "Collect rnaseq metrics ..." 
	in=$outputPath/alignment/$sample.$aligner.merged.bam
	out=$outputPath/alignment_stats/picard/$sample.$aligner.merged.bam
	java -Xmx6G -jar ${PICARD} CollectRnaSeqMetrics I=$in O=$out.rnaseq_metrics REF_FLAT=$refFLAT RIBOSOMAL_INTERVALS=$refRibosomeRNAINTERVAL STRAND=$picardStrand CHART=$out.rnaseq.pdf METRIC_ACCUMULATION_LEVEL=SAMPLE TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 

	echo "Estimate library complexity ..." 
	java -Xmx6G -jar ${PICARD} EstimateLibraryComplexity I=$in O=$out.lib_complex_metrics TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 

	# echo "## --------------------------------------"

	# ## rseqc: collect rnaseq matrics 
	# ## essential functions
	# echo "[ " `date` " ] Collect alignment metrics: RSeQC essential functions "
	# in=$outputPath/alignment/$sample.$aligner.merged.bam
	# out=$outputPath/alignment_stats/rseqc/$sample.$aligner.merged.bam

	# echo "Clipping_profile ... "
	# clipping_profile.py -i $in -o $out

	# echo "Infer_experiment ... "
	# infer_experiment.py -i $in -r $refBED12 -s 200000 > $out.infer_experiment

	# # echo "GeneBody_coverage ..." ## needs LOTs of mem (> 10Gb); job killed at this step!
	# # geneBody_coverage.py -i $in -r $refBED12 -o $out 

	# # echo "RPKM_saturation ..." 
	# # RPKM_saturation.py -i $in -o $out -r $refBED12 -c 0.01 

	# echo "## --------------------------------------"

	# ## bedtools: coverage bigwig
	# echo "[ " `date` " ] Calculating read coverage in bigwig: bedtools "
	# in=$outputPath/alignment/$sample.$aligner.merged.bam
	# out=$outputPath/alignment_stats/bedtools/$sample.$aligner.merged.bam.cov

	# scaled=50000000
	# scaleFactor=`awk -v s=$scaled '$2==255 {print s/$1}' $outputPath/alignment/$sample.$aligner.merged.bam.mapq_hist`
	# echo -e "scaleFactor = $scaleFactor"

	# echo "Calculating genome coverage ... "
	# bedtools genomecov -bga -split -ibam $in -scale $scaleFactor > $out.bdg

	# echo "Converting bedgraph to bigwig ... "
	# bedGraphToBigWig $out.bdg $chromSize $out.bw

fi 

## =========================== PART V ================================ ## 

# echo "## --------------------------------------"

# ## cufflinks
# echo "[ " `date` " ] Count cufflinks reads: cuffquant"
# in=$outputPath/alignment/$sample.$aligner.merged.bam
# out=$outputPath/read_counts/cufflinks
# cuffquant -o $out -b $refGenome -u -p $threads $refGTF $in 

# echo "## --------------------------------------"

# ## rsem 
# echo "[ " `date` " ] Count rsem reads: rsem"
# in1=$outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz
# in2=$outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz
# out=$outputPath/read_counts/rsem2/$readgroup.bt2.rsem.pe ## bt2 aligner, rsem reference

# echo "Align PE reads to RSEM reference transcriptome using Bowtie2 aligner and calculate expression values.... "
# rsem-calculate-expression -p $threads --paired-end --bowtie2 --estimate-rspd --output-genome-bam --time $in1 $in2 $rsemIndex $out

# # echo "Generate expression values .... "
# # rsem-calculate-expression -p $threads --paired-end --bam --estimate-rspd --append-names --output-genome-bam --time $out.bam $rsemIndex $out.bam

echo "## --------------------------------------"


## kallisto 
## Note that if you use --pseudobam, only 1 thread is allowd "Error: pseudobam is not compatible with running on many threads"

if [ $aligner == 'kallisto' ]; then 
	echo "Aligner = $aligner"

	echo `date` "Map reads to human transcriptome: kallisto"
	# kallistoIndex=Homo_sapiens/KALLISTOtranscriptome/GRCh38.primary_Gencode32/gencode.v32.transcripts.idx
	## if you want to run kallisto [--fusion] at the same time ...
	kallistoIndex=Homo_sapiens/KALLISTOtranscriptome/GRCh38.primary_Gencode32_slim_maskPAR/gencode.v32.transcripts.slim_header.maskPAR.idx
	if [ $genome == 'grch37' ]; then 
		kallistoIndex=Homo_sapiens/KALLISTOtranscriptome/GRCh37.p13_Gencode19_slim_maskPAR/gencode.v19.pc_transcripts.slim_header.maskPAR.idx
	fi 
	echo "kallistoIndex = $kallistoIndex"

	strandness=" "
	in=

	if [ $strand == 'reverse' ];  then strandness=" --rf-stranded"; fi 
	if [ $strand == 'forward' ];  then strandness=" --fr-stranded"; fi 
	echo "kallisto strandness = $strandness"

	if [ $alignRawReads -eq 1 ]; then 
		if [ $libraryType == 'SE' ]; then
			in=''
			out=$outputPath/read_counts/kallisto/$sample.raw

			## process multiple readgroups together!! coz cannot merge counts afterwards
			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nreadgroup = $readgroup"

				in="$in $r1"

				echo "## --------------------------------------"

				## fastqc (raw reads)
				echo "[ " `date` " ] Collect raw reads qc metrics: fastqc"
				fastqc --extract -o $outputPath/qc_reports -t $threads --nogroup $r1

			done 

			echo "in = $in"
			kallisto quant --seed=10 $strandness -i $kallistoIndex -t $threads --bias -b 100 -o $out --single -l $insertSizeMean -s $insertSizeStd $in 
		fi 

		if [ $libraryType == 'PE' ]; then
			in=''
			out=$outputPath/read_counts/kallisto/$sample.raw

			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				r2=`echo $r2files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nr2 = $r2\nreadgroup = $readgroup"

				in="$in $r1 $r2"

				echo "## --------------------------------------"

				## fastqc (raw reads)
				echo "[ " `date` " ] Collect raw reads qc metrics: fastqc"
				fastqc --extract -o $outputPath/qc_reports -t $threads --nogroup $r1 $r2

			done 

			echo "in = $in"
			kallisto quant --seed=10 $strandness -i $kallistoIndex -t $threads --fusion --bias -b 100 -o $out $in 
		fi 
	fi 

	if [ $alignRawReads -eq 0 ]; then 
		if [ $libraryType == 'SE' ]; then 
			in=''
			out=$outputPath/read_counts/kallisto/$sample.trim

			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nreadgroup = $readgroup"

				in="$in $outputPath/clean_reads/$readgroup.R1.trim.fq.gz"
			done 

			echo "Align SE reads .... "
			echo "in = $in"
			kallisto quant --seed=10 $strandness -i $kallistoIndex -t $threads -bias -b 100 -o $out --single -l $insertSizeMean -s $insertSizeStd $in
		fi 

		if [ $libraryType == 'PE' ]; then 
			in1=''
			in2=''
			out1=$outputPath/read_counts/kallisto/$sample.trim.pe
			out2=$outputPath/read_counts/kallisto/$sample.trim.se

			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				r2=`echo $r2files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nr2 = $r2\nreadgroup = $readgroup"

				in1="$in1 $outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz $outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz"
				in2="$in2 $outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz"
			done 

			echo "Align PE reads ... "
			echo "in1 = $in1"
			kallisto quant --seed=10 $strandness -i $kallistoIndex -t $threads --fusion --bias -b 100 -o $out1 $in1 

			echo "Align SE reads ... "
			echo "in2 = $in2"
			kallisto quant --seed=10 $strandness -i $kallistoIndex -t $threads -bias -b 100 -o $out2 --single -l $insertSizeMean -s $insertSizeStd $in2

		fi 
	fi
fi 

echo "## --------------------------------------"

## salmon 
## map to human transcriptome
if [ $aligner == 'salmon' ]; then 
	echo "Aligner = $aligner"

	echo `date` "Map raw reads to human transcriptome: salmon"
	salmonIndex=Homo_sapiens/SALMONtranscriptome/GRCh38.primary_Gencode32_slim_maskPAR/gencode.v32.transcripts.slim_header.maskPAR.index
	if [ $genome == 'grch37' ]; then 
		salmonIndex=
	fi 
	echo "salmonIndex = $salmonIndex"

	# strandness=A
	strandness=' '
	aux=aux
	in=

	if [ $libraryType == 'SE' ]; then 
		strandness=' U'
		if [ $strand == 'reverse' ]; then strandness=" SR"; fi 
		if [ $strand == 'forward' ]; then strandness=" SF"; fi 
	fi 
	if [ $libraryType == 'PE' ]; then 
		strandness=' IU'
		if [ $strand == 'reverse' ]; then strandness=" ISR"; fi 
		if [ $strand == 'forward' ]; then strandness=" ISF"; fi 
	fi 
	echo "salmon strandness = $strandness"

	if [ $alignRawReads -eq 1 ]; then 
		in1=''
		in2=''
		in=''
		out=$outputPath/read_counts/salmon/$sample.raw

		## process multiple readgroups together!! coz cannot merge counts afterwards
		for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
			r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
			r2=`echo $r2files | awk -v i=$rgline '{print $i}'`
			readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

			echo -e "r1 = $r1\nr2 = $r2\nreadgroup = $readgroup"

			if [ $libraryType == 'SE' ]; then in1="$in1 $r1"; fi 
			if [ $libraryType == 'PE' ]; then in1="$in1 $r1"; in2="$in2 $r2"; fi 
		done 

		if [ $libraryType == 'SE' ]; then in=" -r $in1"; fi 
		if [ $libraryType == 'PE' ]; then in=" -1 $in1 -2 $in2"; fi 
		
		echo "in = $in"
		salmon quant -i $salmonIndex -l $strandness -o $out -p $threads -g $refGTF --auxDir $aux --writeUnmappedNames --dumpEq --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias --biasSpeedSamp 10 $in

	fi 

	if [ $alignRawReads -eq 0 ]; then 
		if [ $libraryType == 'SE' ]; then 
			in1=''
			in=''
			out=$outputPath/read_counts/salmon/$sample.trim

			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nreadgroup = $readgroup"

				in1="$in1 $outputPath/clean_reads/$readgroup.R1.trim.fq.gz"
			done 

			in=" -r $in1"
			
			echo "Align SE reads .... "
			echo "in = $in"
			salmon quant -i $salmonIndex -l $strandness -o $out -p $threads -g $refGTF --auxDir $aux --writeUnmappedNames --dumpEq --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias --biasSpeedSamp 10 $in
		fi 
		if [ $libraryType == 'PE' ]; then 
			in1=''
			in2=''
			in3=''
			in=''
			out1=$outputPath/read_counts/salmon/$sample.trim.pe
			out2=$outputPath/read_counts/salmon/$sample.trim.se

			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				r2=`echo $r2files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nr2 = $r2\nreadgroup = $readgroup"

				in1="$in1 $outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz"
				in2="$in2 $outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz"
				in3="$in3 $outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz"
			done 

			in=" -1 $in1 -2 $in2"

			echo "Align PE reads .... "
			echo "in = $in"
			salmon quant -i $salmonIndex -l $strandness -o $out1 -p $threads -g $refGTF --auxDir $aux --writeUnmappedNames --dumpEq --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias --biasSpeedSamp 10 $in

			in=" -r $in3"

			echo "Align SE reads .... "
			echo "in = $in"
			salmon quant -i $salmonIndex -l $strandness -o $out2 -p $threads -g $refGTF --auxDir $aux --writeUnmappedNames --dumpEq --useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias --biasSpeedSamp 10 $in
		fi 
	fi

fi 

echo "## --------------------------------------"

## star2pass
## for variant calling

if [ $aligner == 'star2pass' ]; then 
	echo "Aligner = $aligner"

	echo "RGstringSTARstring = $RGstringSTARstring"
	echo "libraryType = $libraryType"

	echo "[ " `date` " ] Map to reference genome: STAR 2-pass"
	if [ $alignRawReads -eq 1 ]; then 
		if [ $libraryType == 'SE' ]; then 
			in1=$loc/$r1
			out1=$outputPath/alignment/$readgroup.$aligner.

			echo "Align SE reads .... "
			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
		fi 
		if [ $libraryType == 'PE' ]; then 
			in1= ## r1,r1,r1
			in2= ## r2,r2,r2
			out1=$outputPath/alignment/$sample.$aligner.

			## process multiple readgroups together!! coz cannot merge counts afterwards
			for (( rgline=1; rgline<=$rgTotal; rgline++ )); do 
				r1=`echo $r1files | awk -v i=$rgline '{print $i}'`
				r2=`echo $r2files | awk -v i=$rgline '{print $i}'`
				readgroup=`echo $readgroups | awk -v i=$rgline '{print $i}'`

				echo -e "r1 = $r1\nr2 = $r2\nreadgroup = $readgroup"

				if [ $rgline -eq 1 ]; then in1=$r1; in2=$r2; fi 
				if [ $rgline -gt 1 ]; then in1="$in1,$r1"; in2="$in2,$r2"; fi 

			done 

			echo -e "in1 = $in1\nin2 = $in2"

			echo "Align PE reads .... "
			STAR --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 --chimOutType WithinBAM SoftClip --chimSegmentMin 15 --genomeDir $starIndex --genomeLoad NoSharedMemory --limitSjdbInsertNsj 1200000 --outFileNamePrefix $out1 --outFilterIntronMotifs None --outFilterMatchNminOverLread 0.33 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 --outFilterType BySJout --outSAMattributes NH HI AS nM NM ch --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat --runThreadN $threads --twopassMode Basic --readFilesIn $in1 $in2 --outSAMattrRGline $RGstringSTARstring

			echo "Sort bam file ... "
			out2=$outputPath/alignment/$sample.$aligner.srt.bam
			sambamba sort --tmpdir $tmpDir -o $out2 -t $threads ${out1}Aligned.out.bam
			sambamba index $out2
		fi 
	fi 
	if [ $alignRawReads -eq 0 ]; then 
		if [ $libraryType == 'SE' ]; then 
			in1=$outputPath/clean_reads/$readgroup.R1.trim.fq.gz
			out1=$outputPath/alignment/$readgroup.$aligner.

			echo "Align SE reads .... "
			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
		fi 
		if [ $libraryType == 'PE' ]; then 
			in1=$outputPath/clean_reads/$readgroup.R1.trim.pe.fq.gz
			in2=$outputPath/clean_reads/$readgroup.R2.trim.pe.fq.gz
			in3=$outputPath/clean_reads/$readgroup.R1R2.trim.unpaired.fq.gz
			out1=$outputPath/alignment/$readgroup.$aligner.pe.
			out2=$outputPath/alignment/$readgroup.$aligner.se.
			out3=$outputPath/alignment/$readgroup.$aligner.Aligned.sortedByCoord.out.bam

			echo "Align PE reads .... "
			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out1 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in1 $in2 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR
			echo "Align SE reads .... "
			STAR --runMode alignReads --genomeLoad NoSharedMemory --outFileNamePrefix $out2 --readFilesCommand zcat --genomeDir $starIndex --readFilesIn $in3 --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $RGstringSTAR

			echo "Merge PE and SE read alignment ..."
			sambamba merge -t $threads $out3 ${out1}Aligned.sortedByCoord.out.bam ${out2}Aligned.sortedByCoord.out.bam
		fi 
	fi
	
	# echo "## --------------------------------------"

	# targetBED=$projPath/results/rnaseq/checkmate_multisample_grch38/gene_variants/FCGR/FCGR.gencode27.grch38.bed
	# gene=FCGR

	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Subset bam by gene : sambamba"
	# in=$outputPath/alignment/$sample.$aligner.srt.bam
	# out=$outputPath/alignment/$sample.$aligner.srt.$gene.bam
	# sambamba view -f bam -L $targetBED -t $threads -o $out $in 

	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Mark duplicates: sambamba"
	# in=$outputPath/alignment/$sample.$aligner.srt.$gene.bam
	# out=$outputPath/alignment/$sample.$aligner.srt.$gene.markdup.bam
	# sambamba markdup -t $threads --tmpdir=$tmpDir $in $out
	# sambamba index $out 

	# echo "## --------------------------------------"

	# ## No GATK tools should run in between star and vardict. SKIP!!
	# ## https://github.com/chapmanb/bcbio-nextgen/issues/389
	# # echo "[ " `date` " ] Deal with Ns in CIGAR, reassign MapQ: Split'N'Trim"
	# # in=$outputPath/alignment/$sample.$aligner.srt.$gene.markdup.bam
	# # out=$outputPath/alignment/$sample.$aligner.srt.$gene.markdup.split.bam
	# # java -jar -Xmx6G ${GATK} -T SplitNCigarReads -R $refGenome -I $in -o $out -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	# echo "## --------------------------------------"

	# echo "[ " `date` " ] Variant calling: VarDict"
	# ## minimum allele frequency
	# minAF=0.10
	# minAltReads=2
	# minMapQ=30
	# minBaseQ=30
	# in=$outputPath/alignment/$sample.$aligner.srt.$gene.markdup.bam
	# out1=$outputPath/alignment/$sample.$aligner.srt.$gene.markdup.vardict.txt
	# out2=${out1/.txt/.vcf}

	# export PATH=/group/bioinformatics/software/VarDict/1.5.1/VarDict:$PATH
	# vardict -G $refGenome -f $minAF -N $sample -b $in -c 1 -S 2 -E 3 -g 4 -F 0x500 -k 1 -r $minAltReads -Q $minMapQ -q $minBaseQ $targetBED | teststrandbias.R > $out1
	# var2vcf_valid.pl -N $sample -E -f $minAF $out1 > $out2 

	echo "## --------------------------------------"

	## picard: collect rnaseq matrics
	echo "[ " `date` " ] Collect alignment metrics: picard CollectMultipleMetrics"

	echo "Collect multiple metrics ..." ## PROGRAM=RnaSeqMetrics will throw nullpoints in java! 
	## cannot process in mixed PE & SE bam files ... need to run on each
	if [ $alignRawReads -eq 1 ]; then 
		in=$outputPath/alignment/$sample.$aligner.srt.bam
		out=$outputPath/alignment_stats/picard/$sample.$aligner.srt.bam

		if [ $libraryType == 'PE' ]; then 
			java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
		fi 
	fi 

	if [ $alignRawReads -eq 0 ]; then 
		echo "readgroups = $readgroups"
		
		if [ $libraryType == 'PE' ]; then 
			for readgroup in $readgroups; do 
				echo "PE read bam ... $readgroup"
				in=$outputPath/alignment/$readgroup.$aligner.pe.Aligned.sortedByCoord.out.bam
				out=$outputPath/alignment_stats/picard/$readgroup.$aligner.pe.bam
				java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000

				echo "SE read bam ... $readgroup"
				in=$outputPath/alignment/$readgroup.$aligner.se.Aligned.sortedByCoord.out.bam
				out=$outputPath/alignment_stats/picard/$readgroup.$aligner.se.bam
				java -Xmx6G -jar ${PICARD} CollectMultipleMetrics I=$in O=$out R=$refGenome PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000
			done 
		fi 
	fi 

	echo "Collect rnaseq metrics ..." 
	in=$outputPath/alignment/$sample.$aligner.srt.bam
	out=$outputPath/alignment_stats/picard/$sample.$aligner.srt.bam
	java -Xmx6G -jar ${PICARD} CollectRnaSeqMetrics I=$in O=$out.rnaseq_metrics REF_FLAT=$refFLAT RIBOSOMAL_INTERVALS=$refRibosomeRNAINTERVAL STRAND=$picardStrand CHART=$out.rnaseq.pdf METRIC_ACCUMULATION_LEVEL=SAMPLE TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 

	echo "Estimate library complexity ..." 
	java -Xmx6G -jar ${PICARD} EstimateLibraryComplexity I=$in O=$out.lib_complex_metrics TMP_DIR=$tmpDir VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 

	echo "## --------------------------------------"

	echo "[ " `date` " ] Mark duplicates: sambamba"
	in=$outputPath/alignment/$sample.$aligner.srt.bam
	out=$outputPath/alignment/$sample.$aligner.srt.dedup.bam
	sambamba markdup -t $threads --tmpdir=$tmpDir $in $out
	sambamba index $out 

	# echo "## --------------------------------------"

	## freaturecounts: count raw reads  
	## featurecount works with BAM of both PE and SE reads (https://groups.google.com/forum/#!topic/subread/qImVtrKtZsw)
	echo "[ " `date` " ] Count raw reads: featureCounts"
	in=$outputPath/alignment/$sample.$aligner.srt.bam 
	out1=$outputPath/read_counts/raw/$sample.$aligner.srt.bam.primary.raw_counts.txt
	out2=$outputPath/read_counts/raw/$sample.$aligner.srt.bam.primary_uniq.raw_counts.txt
	out3=$outputPath/read_counts/raw/$sample.$aligner.srt.bam.mapq1.raw_counts.txt

	## Note that "Assigned + Unassigned_Ambiguity + Unassigned_NoFeatures == number of reads with MapQ 255"
	## STAR manual: uniquely mapped reads are assigned with mapq 255
	echo "Count primary alignments and uniquely mapped reads only ..."
	featureCounts -s $featureCountsStrand -a $refGTF -o $out2 -p -t exon -g gene_id -Q 255 -T $threads -J -G $refGenome --primary -C $in

	echo "Count primary alignments only ..."
	featureCounts -s $featureCountsStrand -a $refGTF -o $out1 -p -t exon -g gene_id -Q 0 -T $threads -J -G $refGenome --primary -C $in

	## Kyle's commands 
	echo "Count reads with mapping quality at least 1; no primary or unique filters ..."
	featureCounts -s $featureCountsStrand -a $refGTF -o $out3 -p -t exon -g gene_id -Q 1 -T $threads -G $refGenome $in

	# echo "## --------------------------------------"

	# ## for gdc samples  only ....
	# echo "[ " `date` " ] Count raw reads: htseq-count"
	# in=$outputPath/alignment/$sample.$aligner.srt.bam 
	# out=$outputPath/read_counts/raw/$sample.$aligner.srt.bam.htseq.raw_counts.txt
	# samtools view -F 4 $in | htseq-count -m intersection-nonempty -i gene_id -r pos -s no - $refGTF > $out
	
fi

echo "## --------------------------------------"

echo "[ " `date` " ] Program finished!"

## CRC HTC: print job log ....
Software/htc_utilities/crc-job-stats.py


