#!/bin/bash
# Run chr4 pangenome personalization & alignments with default parameters
# Usage: default_param_alignments.sh <hap name>
# Example: default_param_alignments.sh.sh HG00099.2

HAP_NAME=$1
CHROM=chr4
SAMPLE_ID=`echo "$HAP_NAME" | cut -f1 -d "." `
PREFIX=${CHROM}.${HAP_NAME}

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/$CHROM
DISTS=./input_data/${CHROM}_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=./input_data/${CHROM}.cenhap_predictions.tsv

READS=$PROJ_DIR/to_align/$PREFIX
KMER_DIR=$TMPDIR/faith_$PREFIX
rm -rf "$KMER_DIR"
mkdir "$KMER_DIR"

SAMPLED_GRAPH=$PROJ_DIR/graph/default/$PREFIX
SAMPLED_ALN=$PROJ_DIR/alignments/default/$PREFIX
GUESS_LOG=${SAMPLED_GRAPH}.guess.log

kmc -k29 -m128 -okff -t16 -hp ${READS}.real.fastq.gz \
    $KMER_DIR/${PREFIX} "$KMER_DIR"

# Sample 5 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $KMER_DIR/${PREFIX}.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 5 --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> "$GUESS_LOG"

# Use logfile to guess
./guess_n_and_cenhap.py --cenhap-table "$CENHAP_TABLE" \
        --dist-matrix "$DISTS" "$GUESS_LOG" &>> "$GUESS_LOG"

n_to_sample=`fgrep Best "$GUESS_LOG" | cut -d " " -f4`

# Perform alignment to personalized graph
vg haplotypes -k $KMER_DIR/${PREFIX}.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${SAMPLED_GRAPH}.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix "$SAMPLED_GRAPH" --no-guessing \
    --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.gbz 2> /dev/null

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.gbz ${READS}.real.fastq.gz "$SAMPLED_ALN"

# More stats!
./param_tests/calculate_private_depths.py -g ${BIG_GRAPH}.gbz.gfa \
    -a ${SAMPLED_ALN}.tsv -l "$GUESS_LOG" &>> "$GUESS_LOG"

cat "$GUESS_LOG"