#!/bin/bash
# Align reads using minimap2 and get mapping statistics
# Usage: align_reads_minimap2.sh <index.mmi> <graph.gbz> <reads.fastq> <output prefix>
# Example: align_reads_minimap2.sh linear_ref_chr12_asat/chm13.chr12asat.mmi \
#                                  linear_ref_chr12_asat/chm13.chr12asat.giraffe.gbz
#                                  HG00408_pat_and_HG00099_hap1_chr12_hor_array.hifi.fastq \
#                                  diploid_subsets/HG00408_pat_and_HG00099_hap1_chr12_hor_array.hifi.chm13.minimap2

INDEX=$1
GRAPH=$2
READS=$3
OUT=$4

minimap2 -ax map-hifi -t 20 $INDEX $READS > $OUT.sam 2> $OUT.log
# Filter out secondary/supplementary alignments before converting to BAM
samtools view -b -F 256 -F 2048 $OUT.sam > $OUT.bam
rm $OUT.sam

# Pull vg's identity statistic
vg inject -x $GRAPH --add-identity $OUT.bam > $OUT.gam
vg filter --tsv-out "name;identity" $OUT.gam > $OUT.tsv