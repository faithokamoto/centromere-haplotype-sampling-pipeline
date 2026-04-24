#!/bin/bash
# Set up inputs for alignment_details_fig.py
# Assumes you've run: retain_files_alignments.sh HG01106.1 chr10
# Usage: massage_inputs.sh

HAP_NAME=HG01106.1
CHROM=chr10
SAMPLE_ID=`echo "$HAP_NAME" | cut -f1 -d "." `
HAP_NUM=`echo "$HAP_NAME" | cut -f2 -d "." `
PREFIX=${CHROM}.${HAP_NAME}

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds
SAMPLE_TMP=$TMPDIR/faith_$SAMPLE_ID

READS=$PROJ_DIR/to_align/$PREFIX

GRAPH_DIR=$PROJ_DIR/graph/haploid
ALN_DIR=$PROJ_DIR/alignments/haploid

SAMPLED_GRAPH=$GRAPH_DIR/${PREFIX}.sampled

OWN_HAP_ALN=$ALN_DIR/${PREFIX}.own_hap
NEIGHBOR_ALN=$ALN_DIR/${PREFIX}.neighbor
CHM13_ALN=$ALN_DIR/${PREFIX}.CHM13
SAMPLED_ALN=$ALN_DIR/${PREFIX}.sampled

# ---- depth ----

# BAMs have already had secondary/supplementary alignments filtered out
samtools sort ${OWN_HAP_ALN}.real.minimap2.bam | samtools depth -a - > ${OWN_HAP_ALN}.depth.tsv
samtools sort ${NEIGHBOR_ALN}.real.minimap2.bam | samtools depth -a - > ${NEIGHBOR_ALN}.depth.tsv
samtools sort ${CHM13_ALN}.real.minimap2.bam | samtools depth -a - > ${CHM13_ALN}.depth.tsv

# Get per-node depths for sampled graph
vg convert --gfa-out ${SAMPLED_GRAPH}.real.gbz > ${SAMPLED_GRAPH}.real.gfa
vg pack -x ${SAMPLED_GRAPH}.real.gbz --gam ${SAMPLED_ALN}.real.giraffe.gam \
    --as-table > ${SAMPLED_ALN}.real.giraffe.pos.depth

./hap_focus/annotate_depth.py -g ${SAMPLED_GRAPH}.real.gfa \
    -d ${SAMPLED_ALN}.real.giraffe.pos.depth > ${SAMPLED_GRAPH}.annot.gfa

# ---- truth coordinates ----

# Get original BAM coordinates for each read
reads=`grep "^$SAMPLE_ID," "./input_data/aws_file_locations.csv" | cut -f3 -d ","`
aws s3 --no-sign-request cp "$reads" $SAMPLE_TMP/full.bam
fgrep chr10 ${BED_DIR}/${HAP_NAME}_asat_arrays.bed > $SAMPLE_TMP/chr10.bed
samtools view -@32 -L $SAMPLE_TMP/chr10.bed -h $SAMPLE_TMP/full.bam > $SAMPLE_TMP/chr10.sam

old_path_name=`cut -f1 $SAMPLE_TMP/chr10.bed`
start=`cut -f2 $SAMPLE_TMP/chr10.bed`
end=`cut -f3 $SAMPLE_TMP/chr10.bed`
chrom=`cut -f4 $SAMPLE_TMP/chr10.bed`
new_path_name="${SAMPLE_ID}#${HAP_NUM}#${HAP_NAME}"

samtools view $SAMPLE_TMP/chr10.sam | sed "s/$old_path_name/$new_path_name/" > $SAMPLE_TMP/no_header.sam
# Filter for reads which appear within the BED file's boundaries
# while also editing the coordinates to be graph-friendly
./get_inputs/edit_sam.py --start "$start" --end "$end" \
    $SAMPLE_TMP/no_header.sam > ${READS}.real.truth.sam

rm -rf $SAMPLE_TMP