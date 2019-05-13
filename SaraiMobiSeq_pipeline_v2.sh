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


# Count reads 













### Im here ####
















































#Wait until the jobs above finish 
echo "Mapping jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy


## MARK DUPLICATES ## CHECK again parameters

# Do the duplicate marking using samtools markdup

#LINE
cd $LINEBAM
if [ ! -e .dup.done ]; then
  for bam in $LINEBAM/*.Wolf_noHets.bam
  do
    bn=$(basename $bam .Wolf_noHets.bam)
    ## sort first by name, so that read1 and read2 are together.
    ## then fixmate, so that the correct md and cigar are present,
    ## then sort it by coordinates, then mark duplicates
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.Wolf_noHets.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J line -R --max-array-jobs=10 --
  touch .dup.done
fi

#SINE
cd $SINEBAM
if [ ! -e .dup.done ]; then
  for bam in $SINEBAM/*.Wolf_noHets.bam
  do
    bn=$(basename $bam .Wolf_noHets.bam)
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.Wolf_noHets.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J Sine -R --max-array-jobs=10 --
  touch .dup.done
fi

#SINE1
cd $SINE1BAM
if [ ! -e .dup.done ]; then
  for bam in $SINE1BAM/*.Wolf_noHets.bam
  do
    bn=$(basename $bam .Wolf_noHets.bam)
    echo "samtools sort -n $bam | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - $bn.Wolf_noHets.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -J Sine1 -R --max-array-jobs=10 --
  touch .dup.done
fiP

#Wait until the jobs above finish 
echo "Mark duplicate jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy














####Im here####









## Merge the rat bams - from different tissues - into one,
## Do this for both the before and after duplicate marked sets of bams
cd $RATBAM
if [ ! -e .merge.done ]; then
  for i in {1..4}
  do
    echo "samtools merge -r -O bam Rat$i.rn6.bam ${i}-*.rn6.bam"
    echo "samtools merge -r -O bam Rat$i.rn6.markdup.bam ${i}-*.rn6.markdup.bam"
  done | xsbatch -c 1 --mem-per-cpu=10G -R --
  touch .merge.done
fi

## Wait for previous commands to be done before moving on.
echo "Rat tissue merge jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy

## This par is only to be done when the bams are ready.
## WAIT FOR BAMS TO BE DONE!!
## For each bam from the markdup file, get the bed with the starting positions of
## read2.
LINEMATCH=$MATCH/wolf/LINE
mkdir -p $LINEMATCH
cd $LINEMATCH
if [ ! -e .match.done ]; then
  for bam in $LINEBAM/*markdup.bam
  do
    ## First remove all the crappy reads - duplicates and secondary alignments (-F 1292)
    ## and restrict to read2 (-f 128), then convert that to a bed file.
    ## and then use the strand information to get the start of read2 in the right orientation
    ## and make a bed with a 20 bp window (will be different depending on the primer
    ## length)
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+20),$4,$5,$6;} $6=="-"{print $1,($3-20),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk  'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  ## Take all the beds, one for each sample, merge them, so that there is at least a 20 bp overlap,
  ## so basically the whole primer overlaps, and make this the bed of primer locations.
  ## After this, using the count of how many samples had this site, make a bed
  ## of locations where at least 9 of the 10 samples contained the site.
  ## the number 8 and 11 are applicable to LINE (90% of 10 samples), but they will differ for the
  ## other primers depending on the number of samples.
  cat *.Wolf_noHets.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -20 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > LINE.allSamples.bed
  awk '$5>8 && $5<11{print $0;}' < LINE.allSamples.bed > LINE.allSamples.90pct.bed
  touch .match.done
fi

SINEMATCH=$MATCH/wolf/SINE
mkdir -p $SINEMATCH
cd $SINEMATCH
if [ ! -e .match.done ]; then
  for bam in $SINEBAM/*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+20),$4,$5,$6;} $6=="-"{print $1,($3-20),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.Wolf_noHets.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -20 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > SINE.allSamples.bed
  awk '$5>8 && $5<11{print $0;}' < SINE.allSamples.bed > SINE.allSamples.90pct.bed
  touch .match.done
fi

DEERMATCH=$MATCH/deer
mkdir -p $DEERMATCH
cd $DEERMATCH
if [ ! -e .match.done ]; then
  for bam in $DEERBAM/*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+24),$4,$5,$6;} $6=="-"{print $1,($3-24),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk  'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.CervusElaphus.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -24 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > BOV2A.allSamples.bed
  ## For the deer, make a bed with the dama dama samples removed.
  cat $(ls *.CervusElaphus.bed | grep -v DD) | bedtools sort | bedtools merge -c 5 -o count -s -d -24 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > BOV2A.onlyCE.bed
  awk '$5>25 && $5<29{print $0;}' < BOV2A.allSamples.bed > BOV2A.allSamples.90pct.bed
  awk '$5>23 && $5<27{print $0;}' < BOV2A.onlyCE.bed > BOV2A.onlyCE.90pct.bed
  touch .match.done
fi

RATMATCH=$MATCH/rats
mkdir -p $RATMATCH
cd $RATMATCH
if [ ! -e .match.done ]; then
  for bam in $RATBAM/R*markdup.bam
  do
    samtools view -O bam -F 1292 -f 128 $bam | bedtools bamtobed -cigar -i /dev/stdin | \
    awk 'BEGIN{OFS="\t";} $6=="+"{ print $1,$2,($2+20),$4,$5,$6;} $6=="-"{print $1,($3-20),$3,$4,$5,$6;}' | \
    bedtools sort | bedtools merge -s -c 5 -o count | awk  'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > $(basename $bam .markdup.bam).bed
  done
  cat *.rn6.bed | bedtools sort | bedtools merge -c 5 -o count -s -d -26 | awk 'BEGIN{OFS="\t";}{print $1,$2,$3,$1"_"$2,$5,$4; }' > L1.allSamples.bed
  awk '$5>3 && $5<5{print $0;}' < L1.allSamples.bed > L1.allSamples.90pct.bed
  touch .match.done
fi
echo "Done generating global matches."

## For each element make a tree with the presence ablsence data, on the full set, and then on the 90pct set.
## This analysis is not included in the manuscript figures.
## using the presence and absence data for each primer, in all the samples,
## make a neighbor joining tree with this information.
cd $LINEMATCH
if [ ! -e .tree.done ]; then
  ## Get the names of the windows in both all loci, and 90% loci.
  cut -f4 LINE.allSamples.bed > LINE.allSamples.presence.txt
  cut -f4 LINE.allSamples.90pct.bed > LINE.allSamples.presence.90pct.txt
  for bed in *.Wolf_noHets.bed; do
    ## For each sample, get the intersection of the sample loci with the
    ## merged for all samples loci.
    bedtools intersect -a LINE.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste LINE.allSamples.presence.txt temp > temp2
    mv temp2 LINE.allSamples.presence.txt
    bedtools intersect -a LINE.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste LINE.allSamples.presence.90pct.txt temp > temp2
    mv temp2 LINE.allSamples.presence.90pct.txt
    rm temp
  done
  touch .tree.done
  ## The matrix with infomration on presence and absence is created here.
  ## Use R and hierarchical clustering to make the NJ tree.
fi

cd $SINEMATCH
if [ ! -e .tree.done ]; then
  cut -f4 SINE.allSamples.bed > SINE.allSamples.presence.txt
  cut -f4 SINE.allSamples.90pct.bed > SINE.allSamples.presence.90pct.txt
  for bed in *.Wolf_noHets.bed; do
    bedtools intersect -a SINE.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste SINE.allSamples.presence.txt temp > temp2
    mv temp2 SINE.allSamples.presence.txt
    bedtools intersect -a SINE.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste SINE.allSamples.presence.90pct.txt temp > temp2
    mv temp2 SINE.allSamples.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

cd $DEERMATCH
if [ ! -e .tree.done ]; then
  cut -f4 BOV2A.allSamples.bed > BOV2A.allSamples.presence.txt
  cut -f4 BOV2A.allSamples.90pct.bed > BOV2A.allSamples.presence.90pct.txt
  for bed in *.CervusElaphus.bed; do
    bedtools intersect -a BOV2A.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.allSamples.presence.txt temp > temp2
    mv temp2 BOV2A.allSamples.presence.txt
    bedtools intersect -a BOV2A.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.allSamples.presence.90pct.txt temp > temp2
    mv temp2 BOV2A.allSamples.presence.90pct.txt
    rm temp
  done
  cut -f4 BOV2A.onlyCE.bed > BOV2A.onlyCE.presence.txt
  cut -f4 BOV2A.onlyCE.90pct.bed > BOV2A.onlyCE.presence.90pct.txt
  for bed in $(ls *.CervusElaphus.bed | grep -v DD); do
    bedtools intersect -a BOV2A.onlyCE.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.onlyCE.presence.txt temp > temp2
    mv temp2 BOV2A.onlyCE.presence.txt
    bedtools intersect -a BOV2A.onlyCE.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste BOV2A.onlyCE.presence.90pct.txt temp > temp2
    mv temp2 BOV2A.onlyCE.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

cd $RATMATCH
if [ ! -e .tree.done ]; then
  cut -f4 L1.allSamples.bed > L1.allSamples.presence.txt
  cut -f4 L1.allSamples.90pct.bed > L1.allSamples.presence.90pct.txt
  for bed in *.rn6.bed; do
    bedtools intersect -a L1.allSamples.bed -b $bed -wa -c -s | cut -f7 > temp
    paste L1.allSamples.presence.txt temp > temp2
    mv temp2 L1.allSamples.presence.txt
    bedtools intersect -a L1.allSamples.90pct.bed -b $bed -wa -c -s | cut -f7 > temp
    paste L1.allSamples.presence.90pct.txt temp > temp2
    mv temp2 L1.allSamples.presence.90pct.txt
    rm temp
  done
  touch .tree.done
fi

## Compute library complexity using preseq.
PRESEQ=$PROJECT/preseq
mkdir -p $PRESEQ
cd $PRESEQ
if [ ! -e .allreads.preseq ]; then
  for bam in $LINEBAM/*Wolf_noHets.bam $SINEBAM/*Wolf_noHets.bam $DEERBAM/*CervusElaphus.bam $RATBAM/R*rn6.bam; do
    bn=$(basename $bam | cut -f1-2 -d ".")
    ## For each bam compute the library complexity using preseq.
    ## -B means bam, -P means paired end, -e is extrapolation range,
    ## and -s is the step size.
    if [ ! -e $bn2.lc ]; then
      echo "preseq lc_extrap -B -P -o $bn.allReads.lc -e 2e9 -s 5e4 $bam"
    fi
  done | xsbatch -c 1 --mem-per-cpu=5G -R --max-array-jobs=156 -J preseq --
  touch .allreads.preseq
fi

## Using the 90% loci, subset each bam so that only the read pairs where read2
## overlaps one of the 90% loci intervals is retained.

## The first step is to create an interval file from the bed file, so that we can use GATK
## This step is not used downstream, so ignore it if you are not going to use GATK.
cd $LINEMATCH
if [ ! -e LINE.allSamples.90pct.intervals ]; then
  java -Xmx2g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=LINE.allSamples.90pct.bed O=LINE.allSamples.90pct.intervals SD=${WOLFGENOME/fasta/dict}
fi
cd $SINEMATCH
if [ ! -e SINE.allSamples.90pct.intervals ]; then
  java -Xmx2g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=SINE.allSamples.90pct.bed O=SINE.allSamples.90pct.intervals SD=${WOLFGENOME/fasta/dict}
fi
cd $DEERMATCH
if [ ! -e BOV2A.allSamples.90pct.intervals ]; then
  java -Xmx4g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=BOV2A.allSamples.90pct.bed O=BOV2A.allSamples.90pct.intervals SD=${DEERGENOME/fasta/dict}
fi
if [ ! -e BOV2A.onlyCE.90pct.intervals ]; then
  java -Xmx4g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=BOV2A.onlyCE.90pct.bed O=BOV2A.onlyCE.90pct.intervals SD=${DEERGENOME/fasta/dict}
fi
cd $RATMATCH
if [ ! -e L1.allSamples.90pct.intervals ]; then
  java -Xmx2g -jar /groups/hologenomics/software/picard/v2.13.2/picard.jar BedToIntervalList I=L1.allSamples.90pct.bed O=L1.allSamples.90pct.intervals SD=${RATGENOME/fasta/dict}
fi

## Use samtools and bedtools to subset the bams to retaine only read pairs that
## have one of the reads overlap one of the intervals in the 90% loci list.
## In the process, also remove the duplicate and secondary alignments
cd $LINEBAM
if [ ! -e .nodupsec.done ]; then
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $LINEMATCH/LINE.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=6g --
  touch .nodupsec.done
fi

cd $SINEBAM
if [ ! -e .nodupsec.done ]; then
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $SINEMATCH/SINE.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=6g --
  touch .nodupsec.done
fi

cd $DEERBAM
if [ ! -e .nodupsec.done ]; then
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $DEERMATCH/BOV2A.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/allSamples.90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=10g --max-array-jobs 12 --
  for bam in *markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $DEERMATCH/BOV2A.onlyCE.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/onlyCE.90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=10g --max-array-jobs 12 --
  touch .nodupsec.done
fi

cd $RATBAM
if [ ! -e .nodupsec.done ]; then
  for bam in Rat*markdup.bam; do
    echo "samtools view -F 1292 -O bam $bam | samtools sort -n -O bam - | bedtools pairtobed -abam /dev/stdin -b $RATMATCH/L1.allSamples.90pct.bed | samtools sort -O bam -o ${bam/markdup.bam/90pct.nodupsec.bam} -"
  done | xsbatch -R -c 1 --mem-per-cpu=6g --
  touch .nodupsec.done
fi

## Wait for these jobs to be finised before continuing.
echo "Bam subsetting jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy

## Generate the bam indices for the recently created bams
cd $LINEBAM
if [ ! -e .bai.done ]; then
  for bam in *bam; do
      echo "samtools index $bam"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
cd $SINEBAM
if [ ! -e .bai.done ]; then
  for bam in *bam; do
      echo "samtools index $bam"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
cd $DEERBAM
if [ ! -e .bai.done ]; then
  for bam in *bam; do
      echo "samtools index $bam"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
cd $RATBAM
if [ ! -e .bai.done ]; then
  for bam in Rat*bam; do
      echo "samtools index $bam"
  done | xsbatch -c 1 --mem-per-cpu=2G -R --
  touch .bai.done
fi
## Again wait for these to get done before moving on.
echo "Generating indices for bams. Press any key to continue, and Ctrl+C to quit."
read dummy

## Make the aggregate plots for coverage around the primer locations, using aggplot.
### Aggregate plot time
AGGDIR=$PROJECT/aggPlot
mkdir -p $AGGDIR
cd $AGGDIR

## Use agplus on the 90pct.nodupsec.bam files to make the reads agplus files,
## which contain the reads starts for the reads that fall along the intervals
## in the 90% loci file. within a +/- 1 kb window.
## do this for each primer.
if [ ! -e .agplus.done ]; then
  for bam in $LINEBAM/*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $LINEMATCH/LINE.allSamples.90pct.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J LINE --
  for bam in $SINEBAM/*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $SINEMATCH/SINE.allSamples.90pct.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J SINE --
  for bam in $DEERBAM/*.allSamples.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.allSamples.90pct.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J BOV2A --
  for bam in $DEERBAM/*.onlyCE.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.onlyCE.90pct.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J BOV2A --
  for bam in $RATBAM/Rat*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam)
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/rn6_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $RATMATCH/L1.allSamples.90pct.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J L1 --
  touch .agplus.done
fi

## Do the agplus computations, but after splitting into + and - strands.
## Split into positive and negative strands
cd $LINEMATCH
if [ ! -e LINE.allSamples.90pct.Plus.bed ]; then
  awk '$6=="+"' < LINE.allSamples.90pct.bed > LINE.allSamples.90pct.Plus.bed
fi
if [ ! -e LINE.allSamples.90pct.Minus.bed ]; then
  awk '$6=="-"' < LINE.allSamples.90pct.bed > LINE.allSamples.90pct.Minus.bed
fi
cd $SINEMATCH
if [ ! -e SINE.allSamples.90pct.Plus.bed ]; then
  awk '$6=="+"' < SINE.allSamples.90pct.bed > SINE.allSamples.90pct.Plus.bed
fi
if [ ! -e SINE.allSamples.90pct.Minus.bed ]; then
  awk '$6=="-"' < SINE.allSamples.90pct.bed > SINE.allSamples.90pct.Minus.bed
fi
cd $DEERMATCH
if [ ! -e BOV2A.allSamples.90pct.Plus.bed ]; then
  awk '$6=="+"' < BOV2A.allSamples.90pct.bed > BOV2A.allSamples.90pct.Plus.bed
fi
if [ ! -e BOV2A.allSamples.90pct.Minus.bed ]; then
  awk '$6=="-"' < BOV2A.allSamples.90pct.bed > BOV2A.allSamples.90pct.Minus.bed
fi
if [ ! -e BOV2A.onlyCE.90pct.Plus.bed ]; then
  awk '$6=="+"' < BOV2A.onlyCE.90pct.bed > BOV2A.onlyCE.90pct.Plus.bed
fi
if [ ! -e BOV2A.onlyCE.90pct.Minus.bed ]; then
  awk '$6=="-"' < BOV2A.onlyCE.90pct.bed > BOV2A.onlyCE.90pct.Minus.bed
fi
cd $RATMATCH
if [ ! -e L1.allSamples.90pct.Plus.bed ]; then
  awk '$6=="+"' < L1.allSamples.90pct.bed > L1.allSamples.90pct.Plus.bed
fi
if [ ! -e L1.allSamples.90pct.Minus.bed ]; then
  awk '$6=="-"' < L1.allSamples.90pct.bed > L1.allSamples.90pct.Minus.bed
fi

## Repeat the aggplot process for the +/- strand bed files.
cd $AGGDIR
if [ ! -e .strand.agplus.done ]; then
  for bam in $LINEBAM/*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam).Plus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $LINEMATCH/LINE.allSamples.90pct.Plus.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
    bn=$(basename $bam .bam).Minus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $LINEMATCH/LINE.allSamples.90pct.Minus.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J LINE --
  for bam in $SINEBAM/*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam).Plus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $SINEMATCH/SINE.allSamples.90pct.Plus.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
    bn=$(basename $bam .bam).Minus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/wolf_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $SINEMATCH/SINE.allSamples.90pct.Minus.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J SINE --
  for bam in $DEERBAM/*.allSamples.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam).Plus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.allSamples.90pct.Plus.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
    bn=$(basename $bam .bam).Minus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.allSamples.90pct.Minus.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J BOV2A --
  for bam in $DEERBAM/*.onlyCE.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam).Plus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.onlyCE.90pct.Plus.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
    bn=$(basename $bam .bam).Minus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/deer_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $DEERMATCH/BOV2A.onlyCE.90pct.Minus.bed -dstart -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J BOV2A --
  for bam in $RATBAM/Rat*.90pct.nodupsec.bam; do
    bn=$(basename $bam .bam).Plus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/rn6_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $RATMATCH/L1.allSamples.90pct.Plus.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
    bn=$(basename $bam .bam).Minus
    if [ ! -e $bn.agplus.txt ]; then
      echo "bam2bwshifted -s 1 -o $bn.bw -g $GENOME/rn6_chrlengths.genome $bam && bigWigToWig $bn.bw $bn.wig && agplus -b $RATMATCH/L1.allSamples.90pct.Minus.bed -d start -o $bn.agplus.txt -r -1000,1000 $bn.wig"
    fi
  done | xsbatch -c 1 --mem-per-cpu=20G -R -J L1 --
  touch .strand.agplus.done
fi
## Wait for these jobs to be done before continuing.
echo "Agplus and agplus +/- jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy

## This section is for computing the statistics for the intervals, like coverage,
## number of samples etc.
mkdir -p $STATS
cd $STATS
## Compute the coverage for the full sets
if [ ! -e .line.done ]; then
  echo "bedtools multicov -bams $(ls $LINEBAM/*markdup.bam | tr -s "\n" " ") -bed $LINEMATCH/LINE.allSamples.bed > LINE.allSamples.bedcov" | xsbatch -c 1 --mem-per-cpu=4G --
  echo "bedtools multicov -bams $(ls $LINEBAM/*nodupsec.bam | tr -s "\n" " ") -bed $LINEMATCH/LINE.allSamples.90pct.bed > LINE.allSamples.90pct.bedcov" | xsbatch -c 1 --mem-per-cpu=4G --
  touch .line.done
fi
if [ ! -e .sine.done ]; then
  echo "bedtools multicov -bams $(ls $SINEBAM/*markdup.bam | tr -s "\n" " ") -bed $SINEMATCH/SINE.allSamples.bed > SINE.allSamples.bedcov" | xsbatch -c 1 --mem-per-cpu=4G --
  echo "bedtools multicov -bams $(ls $SINEBAM/*nodupsec.bam | tr -s "\n" " ") -bed $SINEMATCH/SINE.allSamples.90pct.bed > SINE.allSamples.90pct.bedcov" | xsbatch -c 1 --mem-per-cpu=4G --
  touch .sine.done
fi
if [ ! -e .bov2a.done ]; then
  echo "bedtools multicov -bams $(ls $DEERBAM/*markdup.bam | tr -s "\n" " ") -bed $DEERMATCH/BOV2A.allSamples.bed > BOV2A.allSamples.bedcov" | xsbatch -c 1 --mem-per-cpu=20G --
  echo "bedtools multicov -bams $(ls $DEERBAM/*markdup.bam | grep -v DD | tr -s "\n" " ") -bed $DEERMATCH/BOV2A.onlyCE.bed > BOV2A.onlyCE.bedcov" | xsbatch -c 1 --mem-per-cpu=10G --
  echo "bedtools multicov -bams $(ls $DEERBAM/*allSamples.90pct.nodupsec.bam | tr -s "\n" " ") -bed $DEERMATCH/BOV2A.allSamples.90pct.bed > BOV2A.allSamples.90pct.bedcov" | xsbatch -c 1 --mem-per-cpu=20G --
  echo "bedtools multicov -bams $(ls $DEERBAM/*onlyCE.90pct.nodupsec.bam | grep -v DD | tr -s "\n" " ") -bed $DEERMATCH/BOV2A.onlyCE.90pct.bed > BOV2A.onlyCE.90pct.bedcov" | xsbatch -c 1 --mem-per-cpu=10G --
  touch .bov2a.done
fi
if [ ! -e .l1.done ]; then
  echo "bedtools multicov -bams $(ls $RATBAM/Rat*markdup.bam | tr -s "\n" " ") -bed $RATMATCH/L1.allSamples.bed > L1.allSamples.bedcov" | xsbatch -c 1 --mem-per-cpu=10G --
  echo "bedtools multicov -bams $(ls $RATBAM/Rat*nodupsec.bam | tr -s "\n" " ") -bed $RATMATCH/L1.allSamples.90pct.bed > L1.allSamples.90pct.bedcov" | xsbatch -c 1 --mem-per-cpu=10G --
  touch .l1.done
fi

## Now call variants in the 90% sets using angsd.
## Variant calling
mkdir -p $ANGSD
cd $ANGSD
LINEVAR=$ANGSD/LINE
mkdir -p $LINEVAR
SINEVAR=$ANGSD/SINE
mkdir -p $SINEVAR
DEERVAR=$ANGSD/deer
mkdir -p $DEERVAR
RATVAR=$ANGSD/rats
mkdir -p $RATVAR

cd $LINEVAR
if [ ! -e .angsd.done ]; then
  ## Make a file with the list of bams that should be used for the variant calling
  ls $LINEBAM/*nodupsec.bam > LINE.bamlist
  ## Use angsd, with min qual 30, min map qual 30,
  ## use GL 1 model (samtools genotype likelihood model), compute mafs, filter
  ## for snps with a p-value of 1e-6, get major minor alleles, output a beagle file
  ## and remove individuals with depth less than 3, and remove sites with less than
  ## 50% inds that are retained (in this case 50% of 10 = 5 samples)
  echo "angsd -bam LINE.bamlist -minq 30 -minmapq 30 -GL 1 -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doglf 2 -minind 5 -minIndDepth 3 -out LINE.allSamples.90pct" | xsbatch -c 1 --mem-per-cpu=5G --
  touch .angsd.done
fi
cd $SINEVAR
if [ ! -e .angsd.done ]; then
  ls $SINEBAM/*nodupsec.bam > SINE.bamlist
  echo "angsd -bam SINE.bamlist -minq 30 -minmapq 30 -GL 1 -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doglf 2 -minind 5 -minIndDepth 3 -out SINE.allSamples.90pct" | xsbatch -c 1 --mem-per-cpu=5G --
  touch .angsd.done
fi
cd $DEERVAR
if [ ! -e .angsd.done ]; then
  ls $DEERBAM/*allSamples.90pct.nodupsec.bam > BOV2A.allSamples.bamlist
  echo "angsd -bam BOV2A.allSamples.bamlist -minq 30 -minmapq 30 -GL 1 -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doglf 2 -minind 14 -minIndDepth 3 -out BOV2A.allSamples.90pct" | xsbatch -c 1 --mem-per-cpu=20G --
  ls $DEERBAM/*onlyCE.90pct.nodupsec.bam | grep -v DD > BOV2A.onlyCE.bamlist
  echo "angsd -bam BOV2A.onlyCE.bamlist -minq 30 -minmapq 30 -GL 1 -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doglf 2 -minind 13 -minIndDepth 3 -out BOV2A.onlyCE.90pct" | xsbatch -c 1 --mem-per-cpu=20G --
  touch .angsd.done
fi
cd $RATVAR
if [ ! -e .angsd.done ]; then
  ls $RATBAM/Rat*nodupsec.bam > L1.bamlist
  echo "angsd -bam L1.bamlist -minq 30 -minmapq 30 -GL 1 -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1 -doglf 2 -minind 2 -minIndDepth 3 -out L1.allSamples.90pct" | xsbatch -c 1 --mem-per-cpu=3G --
  touch .angsd.done
fi
## Angsd variant calling in progress. Wait before doing the next steps.
echo "Angsd variant calling jobs launched, wait for them to finish before proceeding. Press Ctrl+C to 1 times to exit. Any other key will continue"
read dummy

## Use the variant calls, beagle files to do ngs analysis.
## Make ngsdist trees
NGSANAL=$PROJECT/ngsAnalysis
mkdir -p $NGSANAL
cd $NGSANAL
module load ngsTools
## Compute the pairwise distances between the samples using the beagle file,
## using ngsDist.
if [ ! -e .ngsdist.done ]; then
  nsites=$(zcat $LINEVAR/LINE.allSamples.90pct.mafs.gz | wc -l)
  let nsites=nsites-1
  while read line; do basename $line .Wolf_noHets.90pct.nodupsec.bam; done < $LINEVAR/LINE.bamlist > LINE.labels
  echo "ngsDist --geno $LINEVAR/LINE.allSamples.90pct.beagle.gz --probs --n_ind 10 --n_sites $nsites --pairwise_del --out LINE.dist --labels LINE.labels --n_boot_rep 100" | xsbatch -c 1 --mem-per-cpu=4G --
  nsites=$(zcat $SINEVAR/SINE.allSamples.90pct.mafs.gz | wc -l)
  let nsites=nsites-1
  while read line; do basename $line .Wolf_noHets.90pct.nodupsec.bam; done < $SINEVAR/SINE.bamlist > SINE.labels
  echo "ngsDist --geno $SINEVAR/SINE.allSamples.90pct.beagle.gz --probs --n_ind 10 --n_sites $nsites --pairwise_del --out SINE.dist --labels SINE.labels --n_boot_rep 100" | xsbatch -c 1 --mem-per-cpu=4G --
  nsites=$(zcat $DEERVAR/BOV2A.allSamples.90pct.mafs.gz | wc -l)
  let nsites=nsites-1
  while read line; do basename $line .CervusElaphus.allSamples.90pct.nodupsec.bam; done < $DEERVAR/BOV2A.allSamples.bamlist > BOV2A.allSamples.labels
  echo "ngsDist --geno $DEERVAR/BOV2A.allSamples.90pct.beagle.gz --probs --n_ind 28 --n_sites $nsites --pairwise_del --out BOV2A.allSamples.dist --labels BOV2A.allSamples.labels --n_boot_rep 100" | xsbatch -c 1 --mem-per-cpu=4G --
  nsites=$(zcat $DEERVAR/BOV2A.onlyCE.90pct.mafs.gz | wc -l)
  let nsites=nsites-1
  while read line; do basename $line .CervusElaphus.onlyCE.90pct.nodupsec.bam; done < $DEERVAR/BOV2A.onlyCE.bamlist > BOV2A.onlyCE.labels
  echo "ngsDist --geno $DEERVAR/BOV2A.onlyCE.90pct.beagle.gz --probs --n_ind 26 --n_sites $nsites --pairwise_del --out BOV2A.onlyCE.dist --labels BOV2A.onlyCE.labels --n_boot_rep 100" | xsbatch -c 1 --mem-per-cpu=4G --
  nsites=$(zcat $RATVAR/L1.allSamples.90pct.mafs.gz | wc -l)
  let nsites=nsites-1
  while read line; do basename $line .rn6.90pct.nodupsec.bam; done < $RATVAR/L1.bamlist > L1.labels
  echo "ngsDist --geno $RATVAR/L1.allSamples.90pct.beagle.gz --probs --n_ind 4 --n_sites $nsites --pairwise_del --out L1.dist --labels L1.labels --n_boot_rep 100" | xsbatch -c 1 --mem-per-cpu=4G --
  touch .ngsdist.done
fi

## using the ngsdistances computed, make neighbor joining tree using fastme and
## raxml-ng to place support on the nodes.
if [ ! -e .tree.done ]; then
  for i in *dist; do
    fastme -i $i -s -D 101 -o $(basename $i .dist).nwk
    grep -v ^$ $(basename $i .dist).nwk > temp
    mv temp $(basename $i .dist).nwk
    head -1 $(basename $i .dist).nwk > $(basename $i .dist).main.nwk
    tail -n +2 $(basename $i .dist).nwk > $(basename $i .dist).boot.nwk
    raxml-ng --support --tree $(basename $i .dist).main.nwk --bs-trees $(basename $i .dist).boot.nwk --prefix $(basename $i .dist)
  done
  touch .tree.done
fi

## Use the maf file from angsd to compute mafs and the number of individuals that
## have data at each site
if [ ! -e .mafinfo.done ]; then
  for i in $ANGSD/*/*mafs.gz; do
    zcat $i | cut -f7 | tail -n +2 > $(basename $i .mafs.gz).inds
    zcat $i | cut -f5 | tail -n +2 > $(basename $i .mafs.gz).mafs
  done
  touch .mafinfo.done
fi

## Compute GC content in a 1kb window around the 90% loci, and then measure coverage
## at those windows
GCANAL=$PROJECT/gcAnalysis
mkdir -p $GCANAL

## LINE processing
## First, extend the interval by 500 bp on both sides, then sort the bed file,
## then use nuc to compute the GC content, then compute coverage at these
## windows.
bedtools slop -i $LINEMATCH/LINE.allSamples.90pct.bed -g $GENOME/wolf_chrlengths.genome -b 500 |
  sort -k1,1 -k2,2n |
  bedtools nuc -fi $WOLFGENOME -bed - | cut -f1-4,8 |
  bedtools coverage -a - -b $LINEBAM/*90pct.nodupsec.bam -mean > LINE.gc_coverage.bed

bedtools slop -i $SINEMATCH/SINE.allSamples.90pct.bed -g $GENOME/wolf_chrlengths.genome -b 500 |
  sort -k1,1 -k2,2n |
  bedtools nuc -fi $WOLFGENOME -bed - | cut -f1-4,8 |
  bedtools coverage -a - -b $SINEBAM/*90pct.nodupsec.bam -mean > SINE.gc_coverage.bed

bedtools slop -i $DEERMATCH/BOV2A.onlyCE.90pct.bed -g $GENOME/deer_chrlengths.genome -b 1000 |
  sort -k1,1 -k2,2n |
  bedtools nuc -fi $DEERGENOME -bed - | cut -f1-4,8 |
  bedtools coverage -a - -b $(ls $DEERBAM/*90pct.nodupsec.bam | grep -v DD) -mean > BOV2A_onlyCE.gc_coverage.bed

bedtools slop -i $RATMATCH/L1.allSamples.90pct.bed -g $GENOME/rn6_chrlengths.genome -b 1000 |
  sort -k1,1 -k2,2n |
  bedtools nuc -fi $RATGENOME -bed - | cut -f1-4,8 |
  bedtools coverage -a - -b $RATBAM/*90pct.nodupsec.bam -mean > L1.gc_coverage.bed

## Finally, compute the distance of these snps, from angsd to the closest primer
## site in the 90% loci list.
cd $ANGSD
zcat LINE/LINE.allSamples.90pct.mafs.gz | awk 'BEGIN{OFS="\t";} NR>1{print $1,$2-1,$2;}' | sort -k1,1 -k2,2n |
  bedtools closest -a - -b $LINEMATCH/LINE.allSamples.90pct.bed -d | cut -f1,3,10 > LINE/LINE.allSamples.90pct.snps.txt
zcat SINE/SINE.allSamples.90pct.mafs.gz | awk 'BEGIN{OFS="\t";} NR>1{print $1,$2-1,$2;}' | sort -k1,1 -k2,2n |
  bedtools closest -a - -b $SINEMATCH/SINE.allSamples.90pct.bed -d | cut -f1,3,10 > SINE/SINE.allSamples.90pct.snps.txt
zcat deer/BOV2A.onlyCE.90pct.mafs.gz    | awk 'BEGIN{OFS="\t";} NR>1{print $1,$2-1,$2;}' | sort -k1,1 -k2,2n |
  bedtools closest -a - -b $DEERMATCH/BOV2A.onlyCE.90pct.bed -d | cut -f1,3,10 > deer/BOV2A.onlyCE.90pct.snps.txt
zcat rats/L1.allSamples.90pct.mafs.gz   | awk 'BEGIN{OFS="\t";} NR>1{print $1,$2-1,$2;}' | sort -k1,1 -k2,2n |
bedtools closest -a - -b $RATMATCH/L1.allSamples.90pct.bed -d | cut -f1,3,10 > rats/L1.allSamples.90pct.snps.txt
