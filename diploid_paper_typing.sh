#!/bin/bash
# Run typing experiments with a diploid sample
# Usage: diploid_paper_typing.sh <sample ID> <chromosome>
# Example: diploid_paper_typing.sh HG00099 chr12

set -e

# ---- process arguments ----

echo "Top of diploid_paper_typing.sh with ${1} ${2} input"

SAMPLE_ID=$1
CHROM=$2
PREFIX=${CHROM}.${SAMPLE_ID}

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/$CHROM
MIRA_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2
DISTS=$MIRA_DIR/all_pairs/distance_matrices/${CHROM}_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=/private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/${CHROM}/${CHROM}.cenhap_predictions.tsv

READS_DIR=$PROJ_DIR/to_align
KMER_DIR=$PROJ_DIR/kmers/$PREFIX

mkdir -p $KMER_DIR

GUESS_LOG=$PROJ_DIR/graph/diploid/${PREFIX}.guess

# ---- run typing ----

# Combine haplotype-specific reads
cat $READS_DIR/${PREFIX}.1.real.fastq $READS_DIR/${PREFIX}.2.real.fastq \
    > $READS_DIR/${PREFIX}.diploid.real.fastq

kmc -k29 -m128 -okff -t16 -hp $READS_DIR/${PREFIX}.diploid.real.fastq \
    $KMER_DIR/${PREFIX}.real "$KMER_DIR"

# Sample 10 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $KMER_DIR/${PREFIX}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 10 --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> ${GUESS_LOG}.real.log

# Use logfile to guess
./guess_n_and_cenhap.py --cenhap-table "$CENHAP_TABLE" --ploidy 2 --fall-threshold 2000 \
        --dist-matrix "$DISTS" ${GUESS_LOG}.real.log &>> ${GUESS_LOG}.real.log
cat ${GUESS_LOG}.real.log

# No need to save these
rm $READS_DIR/${PREFIX}.diploid.real.fastq
rm -rf $KMER_DIR