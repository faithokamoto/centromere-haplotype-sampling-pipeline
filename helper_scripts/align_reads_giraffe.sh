#!/bin/bash
# Align reads using vg giraffe and get mapping statistics
# Usage: align_reads_giraffe.sh <graph.gbz> <reads.fastq> <output prefix>
# Example: align_reads_giraffe.sh linear_ref_chr12_asat/chm13.chr12asat.giraffe.gbz\
#                                 HG00408_pat_and_HG00099_hap1_chr12_hor_array.hifi.fastq \
#                                 diploid_subsets/HG00408_pat_and_HG00099_hap1_chr12_hor_array.hifi.chm13.giraffe

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
GRAPH=$1
READS=$2
OUT=$3

vg giraffe --progress -b hifi --fastq-in $READS --gbz-name $GRAPH -t 20 > $OUT.gam 2> $OUT.log
vg filter --tsv-out "name;identity;nodes" $OUT.gam > $OUT.tsv