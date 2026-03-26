#!/bin/bash
# Run haplotype sampling with a leave-one-out graph
# Usage: haploid_paper_alignments.sh <original path name> <chromosome>
# Example: haploid_paper_alignments.sh HG00099.1 chr12

set -e

# ---- process arguments ----

echo "Top of haploid_paper_alignments.sh with ${1} ${2} input"

ORIG_PATH_NAME=$1
CHROM=$2
SAMPLE_ID=`echo "$ORIG_PATH_NAME" | cut -f1 -d "." `
HAPLO_NUM=`echo "$ORIG_PATH_NAME" | cut -f2 -d "." `
PATH_NAME="${SAMPLE_ID}#${HAPLO_NUM}#${ORIG_PATH_NAME}#0"

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/$CHROM
BAM_LOCS=$PROJ_DIR/to_align/aws_file_locations.csv
MIRA_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2
BED_DIR=$MIRA_DIR/per_smp_asat_beds
DISTS=$MIRA_DIR/all_pairs/distance_matrices/${CHROM}_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=/private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/${CHROM}/${CHROM}.cenhap_predictions.tsv

READS=$PROJ_DIR/to_align/${CHROM}.${ORIG_PATH_NAME}

GRAPH_DIR=$PROJ_DIR/graph/haploid
ALN_DIR=$PROJ_DIR/alignments/haploid
KMER_DIR=$PROJ_DIR/to_align/kmers/$ORIG_PATH_NAME

mkdir -p $KMER_DIR

nearest_neighbor=`grep "$ORIG_PATH_NAME" "$DISTS" | sed 's/,/\t/g' | sort -k3 -n | head -1 \
    | cut -f1-2 | tr "\t" "\n" | grep -v "$ORIG_PATH_NAME"`
echo "Nearest neighbor: $nearest_neighbor"

CHM13_GRAPH=$GRAPH_DIR/${CHROM}.chm13
OWN_HAP_GRAPH=$GRAPH_DIR/${CHROM}.${ORIG_PATH_NAME}
NEIGHBOR_GRAPH=$GRAPH_DIR/${CHROM}.${nearest_neighbor}
SAMPLED_GRAPH=$GRAPH_DIR/${CHROM}.${ORIG_PATH_NAME}.sampled

OWN_HAP_ALN=$ALN_DIR/${CHROM}.${ORIG_PATH_NAME}.own_hap
NEIGHBOR_ALN=$ALN_DIR/${CHROM}.${ORIG_PATH_NAME}.neighbor
CHM13_ALN=$ALN_DIR/${CHROM}.${ORIG_PATH_NAME}.chm13
SAMPLED_ALN=$ALN_DIR/${CHROM}.${ORIG_PATH_NAME}.sampled

GUESS_LOG=$GRAPH_DIR/${CHROM}.${ORIG_PATH_NAME}.guess

# ---- basic alignments (real reads) ----

echo "====="
echo "Linear reference alignments for real reads"

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh "$OWN_HAP_GRAPH" \
    ${READS}.real.fastq ${OWN_HAP_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${READS}.real.fastq ${OWN_HAP_ALN}.real.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh "$CHM13_GRAPH" \
    ${READS}.real.fastq ${CHM13_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${READS}.real.fastq ${CHM13_ALN}.real.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh "$NEIGHBOR_GRAPH" \
    ${READS}.real.fastq ${NEIGHBOR_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${READS}.real.fastq ${NEIGHBOR_ALN}.real.giraffe

# ---- basic alignments (sim reads) ----

echo "====="
echo "Linear reference alignments for sim reads"

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh "$OWN_HAP_GRAPH" \
    ${READS}.sim.fastq ${OWN_HAP_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${READS}.sim.fastq ${OWN_HAP_ALN}.sim.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh "$CHM13_GRAPH" \
    ${READS}.sim.fastq ${CHM13_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${READS}.sim.fastq ${CHM13_ALN}.sim.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh "$NEIGHBOR_GRAPH" \
    ${READS}.sim.fastq ${NEIGHBOR_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${READS}.sim.fastq ${NEIGHBOR_ALN}.sim.giraffe

# ---- align to to haplotype-sampled graphs (real) ----

echo "====="
echo "Haplotype sampling on real reads"

kmc -k29 -m128 -okff -t16 -hp ${READS}.real.fastq \
    $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.real "$KMER_DIR"

# Sample 5 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 5 --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> ${GUESS_LOG}.real.log

# Use logfile to guess
./guess_n_and_cenhap.py --cenhap-table "$CENHAP_TABLE" \
        --dist-matrix "$DISTS" ${GUESS_LOG}.real.log &>> ${GUESS_LOG}.real.log
cat ${GUESS_LOG}.real.log

n_to_sample=`fgrep Best ${GUESS_LOG}.real.log | cut -d " " -f4`

vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${SAMPLED_GRAPH}.real.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix ${SAMPLED_GRAPH}.real --no-guessing \
    --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.real.gbz 2> /dev/null

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.real.gbz ${READS}.real.fastq ${SAMPLED_ALN}.real.giraffe

# ---- align to to haplotype-sampled graphs (sim) ----

echo "====="
echo "Haplotype sampling on sim reads"

kmc -k29 -m128 -okff -t16 -hp ${READS}.sim.fastq \
    $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.sim "$KMER_DIR"

# Sample 5 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.sim.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 5 --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> ${GUESS_LOG}.sim.log

# Use logfile to guess
./guess_n_and_cenhap.py --cenhap-table "$CENHAP_TABLE" \
        --dist-matrix "$DISTS" ${GUESS_LOG}.sim.log &>> ${GUESS_LOG}.sim.log
cat ${GUESS_LOG}.sim.log

n_to_sample=`fgrep Best ${GUESS_LOG}.sim.log | cut -d " " -f4`

vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.sim.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${SAMPLED_GRAPH}.sim.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix ${SAMPLED_GRAPH}.sim --no-guessing \
    --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.sim.gbz 2> /dev/null

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.sim.gbz ${READS}.sim.fastq ${SAMPLED_ALN}.sim.giraffe