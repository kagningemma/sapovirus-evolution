#!/bin/bash 
echo "GANBATTE! Kagning Tsinda Emmanuel"

# index ref
bwa index 2884_phil_batch2.consensus.fa
echo "reference indexed"

# map reads to reference
bwa mem -t 32  2884_phil_batch2.consensus.fa  2926_S8_L001_R1_001.fastq.gz 2926_S8_L001_R2_001.fastq.gz | samtools view -b -F 4 -F 2048 | samtools sort -o 2926_phil_batch2_a.sorted.bam

tput setaf 1 echo "mapping completed"

#Let's now trim off the primer sequences using ivar. In order to do this we need three files.BED file with primer coordinates. We will generate this by aligning the primer sequences to the reference sequence and then using bedtools to create the BED file. Aligned and sorted BAM file generate in the previous step.


bwa mem -k 5 -T 16 2884_phil_batch2.consensus.fa Primalprimers_SaV_GI.1.fa  | samtools view -b -F 4 > Primalprimers_SaV_GI.1.fa.bam

bedtools bamtobed -i Primalprimers_SaV_GI.1.fa.bam > 2926_phil_batch2_primers.bed

tput setaf 1 echo "bed file with primers locations on ref generated"
#"We will now use this BED file as input to ivar to trim primer sequences"

# Note that this BED file has to be generated only once for a reference sequence.

ivar trim -b 2926_phil_batch2_primers.bed -p 2926_phil_batch2.trimmed -i 2926_phil_batch2_a.sorted.bam

echo "primers trimming complete"

# Sort and index trimmed BAM file.
samtools sort -o 2926_phil_batch2.trimmed.sorted.bam 2926_phil_batch2.trimmed.bam ; samtools index 2926_phil_batch2.trimmed.sorted.bam
echo "Sorted and indexed trimmed BAM file"
#Let's quicky take a look at the depth of the trimmed vs untrimmed BAM file. We'll extract the depth using the samtools depth command.
mkdir depth
samtools depth -a 2926_phil_batch2.trimmed.sorted.bam > depth/2926_phil_batch2.trimmed.sorted.bam.depth ; samtools depth -a 2926_phil_batch2_a.sorted.bam > depth/2926_phil_batch2.sorted.bam.depth
tput setaf 1 echo "please plot the depth per position"

# "now we need to identify primer sequences that might have a mismatch with the consensus sequence to ensure that we remove reads from any amplicon that might bias the iSNV frequency due to varying primer binding effeciency"

#To do this, we do this following steps,

#Call consensus on merged BAM file.

#Align primer sequences to consensus after creating a bwa index from the consensus sequence called.


samtools mpileup -A -d 0 -Q 0 2926_phil_batch2.trimmed.sorted.bam  | ivar consensus -p 2926_phil_batch2.consensus
echo "consensus called on merged BAM file."
bwa index -p 2926_phil_batch2.consensus 2926_phil_batch2.consensus.fa

bwa mem -k 5 -T 16 2926_phil_batch2.consensus Primalprimers_SaV_GI.1.fa  | samtools view -bS -F 4 | samtools sort -o 2926_phil_batch2_primers_consensus.bam

tput setaf 1 echo "primers mapped to consensus"
# "Let's now call iSNVs on this BAM file at a minimum threshold of 3% and the default minimum quality threshold of 20"

samtools mpileup -A -d 0 --reference 2926_phil_batch2.consensus.fa -Q 0 2926_phil_batch2_primers_consensus.bam | ivar variants -p 2926_phil_batch2_primers_consensus_1 -t 0.03

tput setaf 1 echo "iSNVs on this raw BAM (without removal of reads from mismached primers) called at a minimum threshold of 3%"

# echo "Let's now get the indices of primers with mismtaches and their respective pairs. To get the pair information, we need a tsv file with two columns to represent the pairs of primers. This file is at pair_information.tsv"

bedtools bamtobed -i 2926_phil_batch2_primers_consensus.bam > 2926_phil_batch2_primers_consensus.bam.bed ; ivar getmasked -i 2926_phil_batch2_primers_consensus_1.tsv -b 2926_phil_batch2_primers_consensus.bam.bed -f Primal-primers_SaVGI.1pairs.tsv -p primer_mismatchers_indices


tput setaf 1 echo "you got the indices of primers with mismtaches and their respective pairs"
# echo "let us remove reads from mismatched primers"

ivar removereads -i 2926_phil_batch2.trimmed.sorted.bam -p 2926_phil_batch2.bad_are_masked.bam -t primer_mismatchers_indices.txt -b 2926_phil_batch2_primers.bed

tput setaf 1 echo "you removed reads from mismatched primers"
samtools sort -o 2926_phil_batch2.bad_are_masked.sorted.bam 2926_phil_batch2.bad_are_masked.bam 
samtools depth -a 2926_phil_batch2.bad_are_masked.sorted.bam > depth/2926_phil_batch2.bad_are_masked.sorted.depth

#echo "Let's now call iSNVs from the BAMS without reads from the masked amplicons"

samtools mpileup -A -d 0 --reference 2884_phil_batch2.consensus.fa -Q 0 2926_phil_batch2.bad_are_masked.sorted.bam | ivar variants -p 2926_phil_batch2_final -t 0.03 ; tput setaf 1 echo "analysis is completed"
tput setaf 1 echo "please exclude variants with depth <400 and count iSNV"
