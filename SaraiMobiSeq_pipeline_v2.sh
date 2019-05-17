#!/bin/bash

# CALL MODULES
#module avail 
module load AdapterRemoval/v2.2.0 
module load perl/v5.24.0
module load prinseq/v0.20.4 
module load bwa/v0.7.15
module load htslib/v1.6
module load samtools/v1.6 
module load bedtools/v20170719

## Define the variables
PROJECT=/groups/hologenomics/sarai/data/batProject
RAWREADS=$PROJECT/rawReads
FASTQ=$PROJECT/fastqs
ADAPTER=$PROJECT/adapRem
NOHOST=$PROJECT/noHost
MAP=$PROJECT/mapping
PREY=$PROJECT/preys

# references
DIREF=/groups/hologenomics/data/genomes
BATGENOME=/groups/hologenomics/sarai/data/my_referenceGenomes/vampire/vampire_bat_10Aug2015_WXtkA.fasta
COWGENOME=${DIREF}/cow/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fasta
PIGGENOME=${DIREF}/pig/susScr11.fasta
SHEEPGENOME=${DIREF}/sheep/oviAri4.fasta
DONKEYGENOME=${DIREF}/donkey/GCF_001305755.1_ASM130575v1_genomic.fasta
HORSEGENOME=${DIREF}/horse/equCab2.masked.fasta
CHICKENGENOME=${DIREF}/chicken/galGal5.fasta
RHINOGENOME=${DIREF}/white_rhino/cerSim1.fasta


############### MERGING ###################
#Merging files of the same lane 
#Dont run with the rest of the pipeline

echo "Merge file of the same line"

## LINE
for f in $(ls *_L006_R1_001.fastq.gz)
do bn=$(basename $f _L006_R1_001.fastq.gz); cat "$bn"_L006_R2_0*.fastq.gz > "$bn"_L006_R2_ME.fastq.gz
done


##### ADAPTERS #####

## Remove the illumina adapters P5 and P7 (AdapterRemoval). 
# Yeah, here we are removing the adapters

#INPUT: fastq files

echo "Remove adapters"
# -p: no error if existing 
mkdir -p $ADAPTER && cd $ADAPTER

## LINE
if [ ! -e .adap.done ]; then
	# For each fastq pair, do cutadapt
	for f in $FASTQ/*_L006_R1_ME.fastq.gz
		# Figure out the name of the sample.
		do
		bn=$(basename $f _L006_R1_ME.fastq.gz)
    		
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
		echo "AdapterRemoval  --qualitybase 33 --qualitybase-output 33 --qualitymax 45 --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
			--collapse --basename ${bn}_primer_noAdap --file1 $f --file2 ${f/R1/R2}"
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J adap --
 	touch .adap.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy



############### REMOVE DUPLICATES #################
# Remove exact duplicates 

echo "Remove duplicated reads" 
cd $ADAPTER

### COLLAPSED READS 
#INPUT: *_primer_noAdap.collapsed.gz and  *_primer_noAdap.collapsed.truncated.gz

if [ ! -e .RD.done ]; then
	## Run prinseq
	# -derep: Type of duplicates to filter. 1 (exact duplicate)
	# -derep_min: This option specifies the number of allowed duplicates. If you want to remove sequence duplicates that occur more than x times
	# -out_format 3, fastq and names
	# uncompress gz file for prinseq
  	for f in $ADAPTER/*_noAdap.collapsed.gz
  		do
    		bn=$(basename $f _noAdap.collapsed.gz)
    		echo "( gzip -dc <(cat $f $ADAPTER/${bn}_noAdap.collapsed.truncated.gz) |  prinseq-lite -verbose -fastq stdin -derep 1 -derep_min 2 -out_format 3 -out_good ${bn}_noAdapCollapsed_noDupPrinSeq -out_bad ${bn}_noAdapCollapsed_DupPrinSeq; gzip $ADAPTER/${bn}_noAdapCollapsed_noDupPrinSeq.fastq $ADAPTER/${bn}_noAdapCollapsed_DupPrinSeq.fastq)"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J prinseq -R --max-array-jobs=10 --
  touch .RD.done
fi

### UNCOLLAPSED READS 
#INPUT: *_primer_noAdap.collapsed.gz and  *_primer_noAdap.collapsed.truncated.gz

cd $ADAPTER
if [ ! -e .RDun.done ]; then
	## Run prinseq
	# -derep: Type of duplicates to filter. 1 (exact duplicate)
	# -derep_min: This option specifies the number of allowed duplicates. If you want to remove sequence duplicates that occur more than x times
	# -out_format 3, fastq and names
  	for f in $ADAPTER/*_noAdap.pair1.truncated.gz
  		do
    		bn=$(basename $f _noAdap.pair1.truncated.gz)
    		echo "(gunzip $f $ADAPTER/${bn}_noAdap.pair2.truncated.gz; prinseq-lite -fastq $ADAPTER/${bn}_noAdap.pair1.truncated -fastq2 $ADAPTER/${bn}_noAdap.pair2.truncated -derep 1 -derep_min 2 -out_format 3 -out_good ${bn}_noAdapUncollapsed_noDupPrinSeq -out_bad ${bn}_noAdapUncollapsed_DupPrinSeq)"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J prinseq_un -R --max-array-jobs=10 --
  touch .RDun.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy

gzip *noAdapUncollapsed*

#Wait until the jobs above finish 
echo "Waiting"
read dummy


################ REMOVE HOST ###############

# Map to the host reference genomes using BWA mem. 
# Use mem so that we can soft clip the primer sequences in the beginning of the read.

# Sort the mapping by coordinates using 
echo "Map to reference genome"
# -p: no error if existing 
mkdir -p $NOHOST && cd $NOHOST

### MAPPING ####
### COLLAPSED READS 

#INPUT: *_noAdapCollapsed_noDupPrinSeq

cd $NOHOST
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
	for f in $ADAPTER/*_noAdapCollapsed_noDupPrinSeq.fastq.gz  
		do     
		bn=$(basename $f _noAdapCollapsed_noDupPrinSeq.fastq.gz)
		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
		echo "(bwa mem -M $BATGENOME $f | samtools view -b -| samtools sort - -o ${bn}_MapBat_collapsed.markdup.bam ; samtools view -b -f 4 $NOHOST/${bn}_MapBat_collapsed.markdup.bam > $NOHOST/${bn}_UNMap_Bat_collapsed.markdup.bam)" 
	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --   
	touch .map.done
fi

### UNCOLLAPSED READS 

if [ ! -e .UNcollapsed.map.done ]; then   
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
	for f in $ADAPTER/*_primer_noAdapUncollapsed_noDupPrinSeq_1.fastq 
		do     
		bn=$(basename $f _primer_noAdapUncollapsed_noDupPrinSeq_1.fastq)     
		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
		echo "(bwa mem -M $BATGENOME $f $ADAPTER/${bn}_primer_noAdapUncollapsed_noDupPrinSeq_2.fastq | samtools view -b -| samtools sort - -o ${bn}_MapBat_uncollapsed.markdup.bam ; samtools view -b -f 4 $NOHOST/${bn}_MapBat_uncollapsed.markdup.bam > $NOHOST/${bn}_UNMap_Bat_uncollapsed.markdup.bam)" 
	done | xsbatch -c 1 --mem-per-cpu=10G -J Unmap -R --max-array-jobs=10 --   
	touch .UNcollapsed.map.done
fi


#Wait until the jobs above finish 
echo "Waiting"
read dummy

### CONVERTING BAM TO FASTQ ###

### COLLAPSED READS 
cd $NOHOST
if [ ! -e .recoverfq.done ]; then for f in $NOHOST/*_UNMap_Bat_collapsed.markdup.bam; do bn=$(basename $f _UNMap_Bat_collapsed.markdup.bam); echo "(bedtools bamtofastq -i $f -fq $NOHOST/${bn}_UNMap_Bat_collapsed.markdup.fastq; gzip $NOHOST/${bn}_UNMap_Bat_collapsed.markdup.fastq )"; done | xsbatch -c 1 --mem-per-cpu=5G -J fq -R --max-array-jobs=10 --; touch .recoverfq.done; fi

### UNCOLLAPSED READS 
cd $NOHOST
if [ ! -e .ucollap.recoverfq.done ]; then for f in $NOHOST/*_UNMap_Bat_uncollapsed.markdup.bam; do bn=$(basename $f _UNMap_Bat_uncollapsed.markdup.bam); echo "(bamToFastq -i $f -fq $NOHOST/${bn}_UNMap_Bat_uncollapsed.markdup.fastq; gzip $NOHOST/${bn}_UNMap_Bat_uncollapsed.markdup.fastq )"; done | xsbatch -c 1 --mem-per-cpu=5G -J uncollpased_fq -R --max-array-jobs=10 --; touch .uncollap.recoverfq.done; fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy

###MERGING#### 
cd $NOHOST
#Merging files of the same lane 
#Dont run with the rest of the pipeline

echo "Merge file of the same line"

for f in $(ls *_UNMap_Bat_collapsed.markdup.fastq.gz)
do bn=$(basename $f _UNMap_Bat_collapsed.markdup.fastq.gz); cat $f ${bn}*fastq.gz > ${bn}_UNMap_BatME.markdup.fastq.gz
done

#Wait until the jobs above finish 
echo "Waiting"
read dummy

################ MAPPING AND GENOME ASSEMBLIES OF DESCRIBED PREYS #################
#Mapping the rest of the reads from host, to every prey separatly

echo "Assemble the preys"
# -p: no error if existing 
mkdir -p $PREY && cd $PREY


### COW ###
PREYCOW=$PREY/cow
mkdir -p $PREYCOW
cd $PREYCOW
if [ ! -e .map.cow ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $COWGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPcow_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.cow
fi

### PIG ###
PREYPIG=$PREY/pig
mkdir -p $PREYPIG
cd $PREYPIG
if [ ! -e .map.pig ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $PIGGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPpig_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.pig
fi

#####################################################here
### SHEEP ###
PREYSHEEP=$PREY/sheep
mkdir -p $PREYSHEEP
cd $PREYSHEEP
if [ ! -e .map.sheep ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $SHEEPGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPsheep_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.sheep
fi

### DONKEY ###
PREYDONKEY=$PREY/donkey
mkdir -p $PREYDONKEY
cd $PREYDONKEY
if [ ! -e .map.donkey ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $DONKEYGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPsheep_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.donkey
fi

### HORSE ###
PREYHORSE=$PREY/horse
mkdir -p $PREYHORSE
cd $PREYHORSE
if [ ! -e .map.horse ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $HORSEGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPsheep_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.horse
fi

### CHICKEN ###
PREYCHICKEN=$PREY/chicken
mkdir -p $PREYCHICKEN
cd $PREYCHICKEN
if [ ! -e .map.chicken ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $CHICKENGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPsheep_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.chicken
fi

### RHINO ###
PREYRHINO=$PREY/rhino
mkdir -p $PREYRHINO
cd $PREYRHINO
if [ ! -e .map.rhino ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $RHINOGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPsheep_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --
  touch .map.rhino
fi









