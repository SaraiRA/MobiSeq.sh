#! /bin/bash

## Load the modules required to run this pipeline on the high performance cluster
module load python/v2.7.12
module load perl/v5.24.0
module load R/v3.4.1
module load java/v1.8.0_131
module load mapDamage/v2.0.6
module load AdapterRemoval/v2.2.2
module load bwa/v0.7.15
module load htslib/v1.6
module load samtools/v1.6
module load paleomix/v1.2.12
module load preseq/v2.0
module load bedtools/v2.26.0
module load kentTools/v22032018
module load picard/v2.13.2
module load agplus/v1.0
module load cutadapt/v1.11
module load angsd/v0.921
module load ngsTools
module load fastme/v2.1.5
module load RAxML-NG/v0.5.1b

## Define the variables
PROJECT=/groups/hologenomics/sarai/data/wolfProject/fecalProject
FASTQ=$PROJECT/fastqs
STATS=$PROJECT/stats
MAP=$PROJECT/mapping
ADAPTER=$PROJECT/adapRem
ANGSD=$PROJECT/angsd
WOLFGENOME=/groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta
ALBA=$PROJECT/albaFiles
ALBAMAP=/groups/hologenomics/ariglesi/data/MobiSeq_wolf_popgen/mapping/
FREBAM=/groups/hologenomics/fmfava/data/wolfproject/subset_highCoverage/

DEPTHS=$PROJECT/depths

DEPTHLOCI=$DEPTHS/loci
DEPTH=$DEPTHS/depth
DEPTHTISSUE=$DEPTHS/depth_tissue
DEPTHFRE=$DEPTHLOCI/depth_frederikke

DEPTHSNP=$DEPTHS/snp
DEPTHSNPFECAL=$DEPTHSNP/fecal
DEPTHSNPTISSUE=$DEPTHSNP/tissue

## Define file variables



##### MERGING #####
#Merging files of the same lane 
#Dont run with the rest of the pipeline

echo "Merge file of the same line"
# -p: no error if existing 
mkdir -p $FASTQ && cd $FASTQ

## LINE
for f in $FASTQ/LINE/*_L001_R1_001.fastq.gz
	do
		echo "Merging"
    		bn=$(basename $f _L001_R1_001.fastq.gz)
    		cat "$bn"_L00*_R1_001.fastq.gz > "$bn"_ME_R1_001.fastq.gz
		cat "$bn"_L00*_R2_001.fastq.gz > "$bn"_ME_R2_001.fastq.gz
	done
rm *L00*

## SINE
for f in $FASTQ/SINE/*_L001_R1_001.fastq.gz
	do
		echo "Merging"
    		bn=$(basename $f _L001_R1_001.fastq.gz)
    		cat "$bn"_L00*_R1_001.fastq.gz > "$bn"_ME_R1_001.fastq.gz
		cat "$bn"_L00*_R2_001.fastq.gz > "$bn"_ME_R2_001.fastq.gz
	done
rm *L00*

## SINE1
for f in $FASTQ/SINE1/*_L001_R1_001.fastq.gz
	do
		echo "Merging"
    		bn=$(basename $f _L001_R1_001.fastq.gz)
    		cat "$bn"_L00*_R1_001.fastq.gz > "$bn"_ME_R1_001.fastq.gz
		cat "$bn"_L00*_R2_001.fastq.gz > "$bn"_ME_R2_001.fastq.gz
	done
rm *L00*


##### ADAPTERS #####

## Verify that all the sequences have the TE target (cutadapt).
# Filter all the sequences that does not end in the 5' of R2 with the TE target, allowing 5% mismatches.
# IMPORTANT: We are not removing the primer sequence, just the sequences that does not have the primer ;)

## Remove the illumina adapters P5 and P7 (AdapterRemoval). 
# Yeah, here we are removing the adapters

#INPUT: fastq files

echo "Filter for TE target and remove adapters"
# -p: no error if existing 
mkdir -p $ADAPTER && cd $ADAPTER

## LINE
LINEADAP=$ADAPTER/LINE
mkdir -p $LINEADAP && cd $LINEADAP
# Check if the adapter removal is done, if yes then this file must exist.
# If it does, skip this step, if not then perform cutadapt steps.
if [ ! -e .adap.done ]; then
	# For each fastq pair, do cutadapt
	for f in $FASTQ/LINE/*LINE*_R1_001.fastq.gz
		do
		# Figure out the name of the sample.
		bn=$(basename $f _R1_001.fastq.gz)
    		
		## Run cutadapt to verify that all the sequences have the TE target. 
    		# -a : check for sequence on the 3' end of read1. It is N, so any base will do.
		# -G : check for sequence on the 5' end of read2. It is the primer for LINE, ^ means it should be in the star. 
  		# --no-trim: do not remove the primer sequences.
    		# --no-indels: do not account for short insertions and deletions.
    		# -e : mismatch rate, here it is 5%, 1 mismatch in the primer sequence.
    		# -o : output name for read 1
  		# -p : output name for read 2
    		# $f and ${f/R1/R2}: read 1 and read 2 input files respectively.
   
 		## Run Adapter removal upon the successful completion of cutadapt
		# collapse: Combined into a single consensus sequence pair aligments
		# Output: output_paired.collapsed containing merged reads,
		# output_paired.collapsed.truncated containing merged reads that have been trimmed due to the --trimns or --trimqualities options.
		# The sequencing primers are not specified, since AdapterRemoval knows the Illumina sequencing primers.
   		echo "cutadapt -a N$ -G ^GATAGCCAAACTGTGGAAGG --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    		AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
		--collapse --basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J LINEadap --
 	touch .adap.done
fi

## SINE
SINEADAP=$ADAPTER/SINE
mkdir -p $SINEADAP && cd $SINEADAP
if [ ! -e .adap.done ]; then
	for f in $FASTQ/SINE/*SINE*_R1_001.fastq.gz
  		do
    		bn=$(basename $f _R1_001.fastq.gz)
    		echo "cutadapt -a N$ -G ^GAGACCCGGGATCGAATCCC --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    		AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 \
		--collapse --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
    		--basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
  		done | xsbatch -c 1 --mem-per-cpu=2G -R -J SINEadap --
  		touch .adap.done
fi

## SINE1
SINE1ADAP=$ADAPTER/SINE1
mkdir -p $SINE1ADAP && cd $SINE1ADAP
if [ ! -e .adap.done ]; then
  	for f in $FASTQ/SINE1/*SINE*_R1_001.fastq.gz
  		do
   	 	bn=$(basename $f _R1_001.fastq.gz)
    		echo "cutadapt -a N$ -G ^CAGAGACCCGGGATCGAATCCC --discard-untrimmed --no-trim --no-indels -e 0.05 -o ${bn}_primer_R1.fastq.gz -p ${bn}_primer_R2.fastq.gz $f ${f/R1/R2} &&\
    		AdapterRemoval --qualitybase 33 --qualitybase-output 33 --qualitymax 45 \
    		--gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
		--collapse --basename ${bn}_primer_noAdap --file1 ${bn}_primer_R1.fastq.gz --file2 ${bn}_primer_R2.fastq.gz"
 		done | xsbatch -c 1 --mem-per-cpu=2G -R -J SINE1adap --
  		touch .adap.done
fi


##### MAPPING AND MARK DUPLICATES #####
# Map to the appropriate reference genomes using BWA mem. 
# Use mem so that we can soft clip the primer sequences in the beginning of the read.

# Sort the mapping by coordinates using 
echo "Map to reference genome"
# -p: no error if existing 
mkdir -p $MAP && cd $MAP

### COLLAPSED READS 

#INPUT: *_primer_noAdap.collapsed.gz and  *_primer_noAdap.collapsed.truncated.gz

## LINE
LINEBAM=$MAP/LINE
mkdir -p $LINEBAM
cd $LINEBAM
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/LINE/*_primer_noAdap.collapsed.gz
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta <(cat $f $ADAPTER/LINE/${bn}_primer_noAdap.collapsed.truncated.gz) | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}.Wolf_noHets_collapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .map.done
fi

## SINE
SINEBAM=$MAP/SINE
mkdir -p $SINEBAM
cd $SINEBAM
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/SINE/*_primer_noAdap.collapsed.gz
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta <(cat $f $ADAPTER/LINE/${bn}_primer_noAdap.collapsed.truncated.gz) | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}.Wolf_noHets_collapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J sine -R --max-array-jobs=10 --
  touch .map.done
fi

## SINE1
SINE1BAM=$MAP/SINE1
mkdir -p $SINE1BAM
cd $SINE1BAM
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/SINE1/*_primer_noAdap.collapsed.gz
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta <(cat $f $ADAPTER/LINE/${bn}_primer_noAdap.collapsed.truncated.gz) | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}.Wolf_noHets_collapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J sine1 -R --max-array-jobs=10 --
  touch .map.done
fi

### UNCOLLAPSED READS  

#INPUT: **_primer_noAdap.pair1.truncated.gz *_primer_noAdap.pair2.truncated.gz

## LINE
cd $LINEBAM
if [ ! -e .map.done.uncollapsed ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/LINE/*_primer_noAdap.pair1.truncated.gz
  		do
    		bn=$(basename $f _primer_noAdap.pair1.truncated.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta $f $ADAPTER/LINE/${bn}_primer_noAdap.pair2.truncated.gz | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}.Wolf_noHets_uncollapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .map.done.uncollapsed
fi

## SINE
cd $SINEBAM
if [ ! -e .map.done.uncollapsed ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/SINE/*_primer_noAdap.pair1.truncated.gz
  		do
    		bn=$(basename $f _primer_noAdap.pair1.truncated.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta $f $ADAPTER/SINE/${bn}_primer_noAdap.pair2.truncated.gz | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}.Wolf_noHets_uncollapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J sine -R --max-array-jobs=10 --
  touch .map.done.uncollapsed
fi

## SINE1
cd $SINE1BAM
if [ ! -e .map.done.uncollapsed ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $ADAPTER/SINE1/*_primer_noAdap.pair1.truncated.gz
  		do
    		bn=$(basename $f _primer_noAdap.pair1.truncated.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M /groups/hologenomics/data/genomes/wolf/L.Dalen_14_wolf.scf.noHets.fasta $f $ADAPTER/SINE1/${bn}_primer_noAdap.pair2.truncated.gz | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}.Wolf_noHets_uncollapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J sine1 -R --max-array-jobs=10 --
  touch .map.done.uncollapsed
fi


### MERGE COLLAPSED AND UNCOLLAPSED READS 

#INPUT: *_ME.Wolf_noHets_collapsed.markdup.bam  *_ME.Wolf_noHets_uncollapsed.markdup.bam 

##LINE
cd $LINEBAM
if [ ! -e .map.done.merge ]; then
	# Run samtools to merge bam files
  	for f in $LINEBAM/*_ME.Wolf_noHets_collapsed.markdup.bam
  		do
    		bn=$(basename $f _ME.Wolf_noHets_collapsed.markdup.bam)
    		# Run samtools to merge bam files
    		echo "(samtools merge ${bn}_ME.Wolf_noHets_MEcollapsed.markdup.bam $f $LINEBAM/${bn}_ME.Wolf_noHets_uncollapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .map.done.merge
fi

## SINE
cd $SINEBAM
if [ ! -e .map.done.merge ]; then
	# Run samtools to merge bam files
  	for f in $SINEBAM/*_ME.Wolf_noHets_collapsed.markdup.bam
  		do
    		bn=$(basename $f _ME.Wolf_noHets_collapsed.markdup.bam)
    		# Run samtools to merge bam files
    		echo "(samtools merge ${bn}_ME.Wolf_noHets_MEcollapsed.markdup.bam $f $SINEBAM/${bn}_ME.Wolf_noHets_uncollapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J sine -R --max-array-jobs=10 --
  touch .map.done.merge
fi

## SINE1
cd $SINE1BAM
if [ ! -e .map.done.merge ]; then
	# Run samtools to merge bam files
  	for f in $SINE1BAM/*_ME.Wolf_noHets_collapsed.markdup.bam
  		do
    		bn=$(basename $f _ME.Wolf_noHets_collapsed.markdup.bam)
    		# Run samtools to merge bam files
    		echo "(samtools merge ${bn}_ME.Wolf_noHets_MEcollapsed.markdup.bam $f $SINE1BAM/${bn}_ME.Wolf_noHets_uncollapsed.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J sine1 -R --max-array-jobs=10 --
  touch .map.done
fi



############################### QUALITY CONTROL ################################ 

######## DEPTH OF LOCI #########
###DEPTH FECAL SAMPLES 
# Quality control, to ensure that we are getting reasonable coverage and data out of the fecal samples.
#Index files first by script
#./index.sh *.bam

#Depth of Coverage
echo "Depth of coverage"
# -p: no error if existing 
mkdir -p $DEPTH && cd $DEPTH

## SINE
#short bed: awk '{print $1 "\t" $2 "\t" $3}' SINE.90pct.collapsed.bed > SINE.90pct.collapsed_correct.bed
SINEQUAL=$DEPTH/SINE
mkdir -p $SINEQUAL
cd $SINEQUAL
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $MAP/SINE/*_MEcollapsed.markdup.bam
  		do
    		bn=$(basename $f _MEcollapsed.markdup.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SINE.90pct.collapsed_correct.bed  $f > ${bn}_SUMdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J sine -R --max-array-jobs=10 --
  touch .depth
fi

## SINE1
#Use same bed as SINE
#short bed: awk '{print $1 "\t" $2 "\t" $3}' SINE.90pct.collapsed.bed > SINE.90pct.collapsed_correct.bed
SINE1QUAL=$DEPTH/SINE1;
mkdir -p $SINE1QUAL;
cd $SINE1QUAL;
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $MAP/SINE1/*_MEcollapsed.markdup.bam;
  		do
    		bn=$(basename $f _MEcollapsed.markdup.bam);
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SINE.90pct.collapsed_correct.bed  $f > ${bn}_SUMdepth.txt)"
		done | xsbatch -c 1 --mem-per-cpu=2G -J sine1 -R --max-array-jobs=10 --
	touch .depth
fi

## LINE
#short bed: awk '{print $1 "\t" $2 "\t" $3}' LINE.90pct.collapsed.bed > LINE.90pct.collapsed_correct.bed
LINEQUAL=$DEPTH/LINE
mkdir -p $LINEQUAL
cd $LINEQUAL
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $MAP/LINE/*_MEcollapsed.markdup.bam
  		do
    		bn=$(basename $f _MEcollapsed.markdup.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/LINE.90pct.collapsed_correct.bed  $f > ${bn}_SUMdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=2G -J line -R --max-array-jobs=10 --
  touch .depth
fi

### DEPTH TISSUE SAMPLES
#*rg.realigned.bam
#Depth of Coverage
echo "Depth of coverage"
# -p: no error if existing 
mkdir -p $DEPTHTISSUE && cd $DEPTHTISSUE

## SINE
#short bed: awk '{print $1 "\t" $2 "\t" $3}' SINE.90pct.collapsed.bed > SINE.90pct.collapsed_correct.bed
SINEDEPTH=$DEPTHTISSUE/SINE
mkdir -p $SINEDEPTH
cd $SINEDEPTH
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $ALBAMAP/SINE_collapsed/*_ME.collapsed.markdup.90pct.rg.realigned.bam
  		do
    		bn=$(basename $f _ME.collapsed.markdup.90pct.rg.realigned.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SINE.90pct.collapsed_correct.bed  $f > ${bn}_SUMdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J sine -R --max-array-jobs=10 --
  touch .depth
fi

## LINE
#short bed: awk '{print $1 "\t" $2 "\t" $3}' LINE.90pct.collapsed.bed > LINE.90pct.collapsed_correct.bed
LINEDEPTH=$DEPTHTISSUE/LINE
mkdir -p $LINEDEPTH
cd $LINEDEPTH
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $ALBAMAP/LINE_collapsed/*_ME.collapsed.markdup.90pct.rg.realigned.bam
  		do
    		bn=$(basename $f _ME.collapsed.markdup.90pct.rg.realigned.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/LINE.90pct.collapsed_correct.bed  $f > ${bn}_SUMdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=2G -J line -R --max-array-jobs=10 --
  touch .depth
fi


### DEPTH SUBSET_FREDERIKA

#No run
#convert sam to bam
for f in $(ls *.sam)
do 
	bn=$(basename $f .sam)
	samtools view -Sb $f > ${bn}.bam
done

#Index bam 
./index.sh *.bam

#Depth of Coverage
echo "Depth of coverage"
# -p: no error if existing 
mkdir -p $DEPTHFRE && cd $DEPTHFRE

## SINE
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $FREBAM/*.bam
  		do
    		bn=$(basename $f .bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
		echo "(samtools bedcov $ALBA/SINE.90pct.collapsed_correct.bed  $f > ${bn}_SUMdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J sine -R --max-array-jobs=10 --
  touch .depth
fi


######## DEPTH OF SNP #########

######Fecal
#Depth of Coverage for the SNPs
echo "Depth of coverage for the SNPs"
# -p: no error if existing 
mkdir -p $DEPTHSNPFECAL && cd $DEPTHSNPFECAL

## SINE
SINEDEPTHSNPFECAL=$DEPTHSNPFECAL/SINE
mkdir -p $SINEDEPTHSNPFECAL && cd $SINEDEPTHSNPFECAL
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $MAP/SINE/*_MEcollapsed.markdup.bam
  		do
    		bn=$(basename $f _MEcollapsed.markdup.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SNPsPos_SINE-merged.mafs.bed  $f > ${bn}_SNPdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J sine -R --max-array-jobs=10 --
  touch .depth
fi

## LINE
LINEDEPTHSNPFECAL=$DEPTHSNPFECAL/LINE
mkdir -p $LINEDEPTHSNPFECAL && cd $LINEDEPTHSNPFECAL
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $MAP/LINE/*_MEcollapsed.markdup.bam
  		do
    		bn=$(basename $f _MEcollapsed.markdup.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SNPsPos_LINE-merged.mafs.bed  $f > ${bn}_SNPdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=2G -J line -R --max-array-jobs=10 --
  touch .depth
fi


######Tissue
#Depth of Coverage for the SNPs
echo "Depth of coverage for the SNPs"
# -p: no error if existing 

mkdir -p $DEPTHSNPTISSUE && cd $DEPTHSNPTISSUE

## SINE
SINEDEPTHSNPTISSUE=$DEPTHSNPTISSUE/SINE
mkdir -p $SINEDEPTHSNPTISSUE && cd $SINEDEPTHSNPTISSUE
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $ALBAMAP/SINE_collapsed/*_ME.collapsed.markdup.90pct.rg.realigned.bam
  		do
    		bn=$(basename $f _ME.collapsed.markdup.90pct.rg.realigned.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SNPsPos_SINE-merged.mafs.bed  $f > ${bn}_SNPdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J sineT -R --max-array-jobs=10 --
  touch .depth
fi

## LINE

LINEDEPTHSNPTISSUE=$DEPTHSNPTISSUE/LINE
mkdir -p $LINEDEPTHSNPTISSUE && cd $LINEDEPTHSNPTISSUE
if [ ! -e .depth ]; then
	## Run samtools bedcov
	# Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. The regions are output as they appear in the BED file and are 0-based. Counts for each alignment file supplied are reported in separate columns. 
  	for f in $ALBAMAP/LINE_collapsed/*_ME.collapsed.markdup.90pct.rg.realigned.bam
  		do
    		bn=$(basename $f _ME.collapsed.markdup.90pct.rg.realigned.bam)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(samtools bedcov $ALBA/SNPsPos_LINE-merged.mafs.bed  $f > ${bn}_SUMdepth.txt)"
  	done | xsbatch -c 1 --mem-per-cpu=2G -J lineT -R --max-array-jobs=10 --
  touch .depth
fi










