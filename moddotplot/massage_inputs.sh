#!/bin/bash
# Set up inputs for moddotplot supplement (in practice, I ran interactively)
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

# ---- vg / samtools; use "cenhap-sample" Conda environment ----

# BAMs have already had secondary/supplementary alignments filtered out
samtools view ${OWN_HAP_ALN}.real.minimap2.bam > ${OWN_HAP_ALN}.clean.sam
samtools view ${NEIGHBOR_ALN}.real.minimap2.bam > ${NEIGHBOR_ALN}.clean.sam
samtools view ${CHM13_ALN}.real.minimap2.bam > ${CHM13_ALN}.clean.sam

# Get per-node depths for sampled graph
vg convert --gfa-out ${SAMPLED_GRAPH}.real.gbz > ${SAMPLED_GRAPH}.real.gfa
vg pack -x ${SAMPLED_GRAPH}.real.gbz --gam ${SAMPLED_ALN}.real.giraffe.gam \
    --as-table > ${SAMPLED_ALN}.real.giraffe.pos.depth

# Average depth per node across offsets https://unix.stackexchange.com/a/465815
awk '{ 
  sum[$2] += $4 
  count[$2] += 1
} 
END { 
  for (key in count) { 
    print key, sum[key] / count[key] 
  } 
}' ${SAMPLED_ALN}.real.giraffe.pos.depth > ${SAMPLED_ALN}.real.giraffe.node.depth

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

# ---- ModDotPlot; use "venv" pip environment ----

# Self-identity plot
moddotplot static --fasta $GRAPH_DIR/$PREFIX.fasta \
    --no-plot --output-dir plot_outputs
# vs. CHM13
moddotplot static --fasta $GRAPH_DIR/$PREFIX.fasta $GRAPH_DIR/$CHROM.CHM13.fasta \
    --compare-only --no-plot --output-dir plot_outputs
# vs. neighbor
moddotplot static --fasta $GRAPH_DIR/$PREFIX.fasta $GRAPH_DIR/$CHROM.HG01891.2.fasta \
    --compare-only --no-plot --output-dir plot_outputs
