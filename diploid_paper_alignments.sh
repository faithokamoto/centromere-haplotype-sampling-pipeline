#!/bin/bash
# Run alignment experiments with a diploid sample
# Usage: diploid_paper_alignments.sh <sample ID> <chromosome>
# Example: diploid_paper_alignments.sh HG00099 chr12

set -e

# ---- process arguments ----

echo "Top of diploid_paper_alignments.sh with ${1} ${2} input"

SAMPLE_ID=$1
CHROM=$2
PREFIX=${CHROM}.${SAMPLE_ID}

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/$CHROM
DISTS=./input_data/${CHROM}_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=./input_data/${CHROM}.cenhap_predictions.tsv

GRAPH_DIR=$PROJ_DIR/graph/diploid
ALN_DIR=$PROJ_DIR/alignments/diploid

READS_DIR=$PROJ_DIR/to_align
SAMPLE_TMP=$TMPDIR/faith_$SAMPLE_ID
mkdir -p $SAMPLE_TMP
DIPLOID_READS=$SAMPLE_TMP/${PREFIX}.diploid.real.fastq.gz

nearest_neighbor_1=`grep "$SAMPLE_ID.1" "$DISTS" | sed 's/,/\t/g' | sort -k3 -n | head -1 \
    | cut -f1-2 | tr "\t" "\n" | grep -v "$SAMPLE_ID.1"`
if [ "$nearest_neighbor_1" == "CHM13.0" ]; then
    # Use name of ref graph
    nearest_path_1="CHM13#0#CHM13.0#0"
else
    neighbor_sample_1=`echo "$nearest_neighbor_1" | cut -f1 -d "."`
    neighbor_haplo_1=`echo "$nearest_neighbor_1" | cut -f2 -d "."`
    nearest_path_1="${neighbor_sample_1}#${neighbor_haplo_1}#${neighbor_sample_1}.${neighbor_haplo_1}#0"
fi
nearest_neighbor_2=`grep "$SAMPLE_ID.2" "$DISTS" | sed 's/,/\t/g' | sort -k3 -n | head -1 \
    | cut -f1-2 | tr "\t" "\n" | grep -v "$SAMPLE_ID.2"`
if [ "$nearest_neighbor_2" == "CHM13.0" ]; then
    # Use name of ref graph
    nearest_path_2="CHM13#0#CHM13.0#0"
else
    neighbor_sample_2=`echo "$nearest_neighbor_2" | cut -f1 -d "."`
    neighbor_haplo_2=`echo "$nearest_neighbor_2" | cut -f2 -d "."`
    nearest_path_2="${neighbor_sample_2}#${neighbor_haplo_2}#${neighbor_sample_2}.${neighbor_haplo_2}#0"
fi

CHM13_GRAPH=$PROJ_DIR/graph/haploid/${CHROM}.CHM13
OWN_HAP_GRAPH=$GRAPH_DIR/${PREFIX}.own_hap
NEIGHBOR_GRAPH=$GRAPH_DIR/${PREFIX}.neighbor
SAMPLED_GRAPH=$GRAPH_DIR/${PREFIX}.sampled

OWN_HAP_ALN=$ALN_DIR/${PREFIX}.own_hap
NEIGHBOR_ALN=$ALN_DIR/${PREFIX}.neighbor
CHM13_ALN=$ALN_DIR/${PREFIX}.CHM13
SAMPLED_ALN=$ALN_DIR/${PREFIX}.sampled

GUESS_LOG=$GRAPH_DIR/${PREFIX}.guess

ABSENT_SCORE=0.05

# Combine haplotype-specific reads
cat $READS_DIR/${PREFIX}.1.real.fastq.gz $READS_DIR/${PREFIX}.2.real.fastq.gz > $DIPLOID_READS

# ---- basic alignments ----

# align to own diplotype
./get_inputs/subset_ref_by_path_name.sh ${BIG_GRAPH}.pg \
    "${SAMPLE_ID}#1#${SAMPLE_ID}.1#0" "${SAMPLE_ID}#2#${SAMPLE_ID}.2#0" "$OWN_HAP_GRAPH"
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    "$DIPLOID_READS" ${OWN_HAP_ALN}.real.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh "$CHM13_GRAPH" \
    "$DIPLOID_READS" ${CHM13_ALN}.real.minimap2

# align to nearest neighbor
./get_inputs/subset_ref_by_path_name.sh ${BIG_GRAPH}.pg \
    "$nearest_path_1" "$nearest_path_2" "$NEIGHBOR_GRAPH"
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    "$DIPLOID_READS" ${NEIGHBOR_ALN}.real.giraffe

# ---- run typing ----

kmc -k29 -m128 -okff -t16 -hp "$DIPLOID_READS" \
    $SAMPLE_TMP/${PREFIX}.real "$SAMPLE_TMP"

# Sample 10 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $SAMPLE_TMP/${PREFIX}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 10 -d ${BIG_GRAPH}.dist --absent-score "$ABSENT_SCORE" \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> ${GUESS_LOG}.real.log

# Use logfile to guess
./guess_n_and_cenhap.py --cenhap-table "$CENHAP_TABLE" --ploidy 2 --fall-threshold 2000 \
        --dist-matrix "$DISTS" ${GUESS_LOG}.real.log &>> ${GUESS_LOG}.real.log
cat ${GUESS_LOG}.real.log

# ---- align against sampled graph ----

n_to_sample=`fgrep Best ${GUESS_LOG}.real.log | cut -d " " -f4`

vg haplotypes -k $SAMPLE_TMP/${PREFIX}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" -d ${BIG_GRAPH}.dist --absent-score "$ABSENT_SCORE" \
    -g ${SAMPLED_GRAPH}.real.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix ${SAMPLED_GRAPH}.real --no-guessing \
    --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.real.gbz 2> /dev/null

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.real.gbz "$DIPLOID_READS" ${SAMPLED_ALN}.real.giraffe

# ---- get stats! ----

./helper_scripts/calculate_alignment_stats.py -c "$CHROM" -n "$SAMPLE_ID" \
    -g ${BIG_GRAPH}.gbz.gfa -r $PROJ_DIR/to_align -a "$ALN_DIR" > $ALN_DIR/${PREFIX}.stats.log

cat $ALN_DIR/${PREFIX}.stats.log

# Clean up behind for space reasons
rm $ALN_DIR/${PREFIX}.*.bam $ALN_DIR/${PREFIX}.*.sam $ALN_DIR/${PREFIX}.*.gam $ALN_DIR/${PREFIX}.*.tsv
rm $GRAPH_DIR/${PREFIX}.own_hap.* $GRAPH_DIR/${PREFIX}.neighbor.* $GRAPH_DIR/${PREFIX}.sampled.*
rm -rf "$KMER_DIR"