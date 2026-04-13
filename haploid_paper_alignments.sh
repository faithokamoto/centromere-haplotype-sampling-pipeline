#!/bin/bash
# Run haplotype sampling with a leave-one-out graph
# Usage: haploid_paper_alignments.sh <hap name> <chromosome>
# Example: haploid_paper_alignments.sh HG00099.1 chr12

set -e

# ---- process arguments ----

echo "Top of haploid_paper_alignments.sh with ${1} ${2} input"

git rev-parse HEAD
vg version

HAP_NAME=$1
CHROM=$2
SAMPLE_ID=`echo "$HAP_NAME" | cut -f1 -d "." `
PREFIX=${CHROM}.${HAP_NAME}

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/$CHROM
DISTS=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/${CHROM}_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=/private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/${CHROM}/${CHROM}.cenhap_predictions.tsv

READS=$PROJ_DIR/to_align/$PREFIX

GRAPH_DIR=$PROJ_DIR/graph/haploid
ALN_DIR=$PROJ_DIR/alignments/haploid

KMER_DIR=$TMPDIR/faith_$PREFIX
mkdir $KMER_DIR

nearest_neighbor=`grep "$HAP_NAME" "$DISTS" | sed 's/,/\t/g' | sort -k3 -n | head -1 \
    | cut -f1-2 | tr "\t" "\n" | grep -v "$HAP_NAME"`
if [ "$nearest_neighbor" == "CHM13.0" ]; then
    # Use name of ref graph
    nearest_neighbor=CHM13
fi
echo "Nearest neighbor: $nearest_neighbor"

CHM13_GRAPH=$GRAPH_DIR/${CHROM}.CHM13
OWN_HAP_GRAPH=$GRAPH_DIR/$PREFIX
NEIGHBOR_GRAPH=$GRAPH_DIR/${CHROM}.${nearest_neighbor}
SAMPLED_GRAPH=$GRAPH_DIR/${PREFIX}.sampled

OWN_HAP_ALN=$ALN_DIR/${PREFIX}.own_hap
NEIGHBOR_ALN=$ALN_DIR/${PREFIX}.neighbor
CHM13_ALN=$ALN_DIR/${PREFIX}.CHM13
SAMPLED_ALN=$ALN_DIR/${PREFIX}.sampled

GUESS_LOG=$GRAPH_DIR/${PREFIX}.guess

ABSENT_SCORE=0.05

# get rid of old graphs / alignments
rm -f ${SAMPLED_GRAPH}.* ${GUESS_LOG}.* $ALN_DIR/${PREFIX}.*

# ---- basic alignments (real reads) ----

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh "$OWN_HAP_GRAPH" \
    ${READS}.real.fastq.gz ${OWN_HAP_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${READS}.real.fastq.gz ${OWN_HAP_ALN}.real.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh "$CHM13_GRAPH" \
    ${READS}.real.fastq.gz ${CHM13_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${READS}.real.fastq.gz ${CHM13_ALN}.real.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh "$NEIGHBOR_GRAPH" \
    ${READS}.real.fastq.gz ${NEIGHBOR_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${READS}.real.fastq.gz ${NEIGHBOR_ALN}.real.giraffe

# ---- basic alignments (sim reads) ----

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh "$OWN_HAP_GRAPH" \
    ${READS}.sim.fastq.gz ${OWN_HAP_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${READS}.sim.fastq.gz ${OWN_HAP_ALN}.sim.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh "$CHM13_GRAPH" \
    ${READS}.sim.fastq.gz ${CHM13_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${READS}.sim.fastq.gz ${CHM13_ALN}.sim.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh "$NEIGHBOR_GRAPH" \
    ${READS}.sim.fastq.gz ${NEIGHBOR_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${READS}.sim.fastq.gz ${NEIGHBOR_ALN}.sim.giraffe

# ---- align to to haplotype-sampled graphs ----

echo "Haplotype sampling on real reads"

kmc -k29 -m128 -okff -t16 -hp ${READS}.real.fastq.gz \
    $KMER_DIR/${PREFIX}.real "$KMER_DIR"

# Sample 5 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $KMER_DIR/${PREFIX}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 5 --haploid-scoring -d ${BIG_GRAPH}.dist --absent-score $ABSENT_SCORE \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> ${GUESS_LOG}.real.log

# Use logfile to guess
./guess_n_and_cenhap.py --cenhap-table "$CENHAP_TABLE" \
        --dist-matrix "$DISTS" ${GUESS_LOG}.real.log &>> ${GUESS_LOG}.real.log
cat ${GUESS_LOG}.real.log

n_to_sample=`fgrep Best ${GUESS_LOG}.real.log | cut -d " " -f4`

vg haplotypes -k $KMER_DIR/${PREFIX}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist --absent-score $ABSENT_SCORE \
    -g ${SAMPLED_GRAPH}.real.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix ${SAMPLED_GRAPH}.real --no-guessing \
    --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.real.gbz 2> /dev/null

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.real.gbz ${READS}.real.fastq.gz ${SAMPLED_ALN}.real.giraffe
# Align sim reads to same graph
./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.real.gbz ${READS}.sim.fastq.gz ${SAMPLED_ALN}.sim.giraffe

# ---- get stats! ----

./helper_scripts/calculate_alignment_stats.py -c "$CHROM" -n "$HAP_NAME" \
    -g ${BIG_GRAPH}.gbz.gfa -r $PROJ_DIR/to_align -a "$ALN_DIR" > $ALN_DIR/${PREFIX}.stats.log

cat $ALN_DIR/${PREFIX}.stats.log

# Clean up behind for space reasons
rm $ALN_DIR/${PREFIX}.*.bam $ALN_DIR/${PREFIX}.*.sam $ALN_DIR/${PREFIX}.*.gam $ALN_DIR/${PREFIX}.*.tsv
rm -rf "$KMER_DIR"