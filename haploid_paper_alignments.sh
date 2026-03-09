#!/bin/bash
# Run haplotype sampling with a leave-one-out graph on chr12
# Usage: haploid_paper_alignments.sh <original path name>
# Example: haploid_paper_alignments.sh HG00099.1

set -e

# ---- process arguments ----

echo "Top of haploid_paper_alignments.sh with ${1} input"

ORIG_PATH_NAME=$1
SAMPLE_ID=`echo "$ORIG_PATH_NAME" | cut -f1 -d "." `
HAPLO_NUM=`echo "$ORIG_PATH_NAME" | cut -f2 -d "." `
PATH_NAME="${SAMPLE_ID}#${HAPLO_NUM}#${ORIG_PATH_NAME}#0"

BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files

ls ${BED_DIR}/${SAMPLE_ID}_hap1_* &>/dev/null
if [ $? -eq 0 ]; then
    if [ "$HAPLO_NUM" -eq "1" ]; then
        SAMPLE_NAME=${SAMPLE_ID}_hap1
    else
        SAMPLE_NAME=${SAMPLE_ID}_hap2
    fi
else
    if [ "$HAPLO_NUM" -eq "1" ]; then
        SAMPLE_NAME=${SAMPLE_ID}_pat
    else
        SAMPLE_NAME=${SAMPLE_ID}_mat
    fi
fi

echo "Processing sample: $SAMPLE_NAME"

# ---- set up variables ----

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr12
DISTS=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=/private/groups/migalab/juklucas/centrolign/notes/correct_cenhaps_chr12/chr12_cenhap_assignments_final.csv

REAL_READS=$PROJ_DIR/to_align/real_${ORIG_PATH_NAME}.chr12.hifi

GRAPH_DIR=$PROJ_DIR/graph/haploid
ALN_DIR=$PROJ_DIR/alignments/haploid
KMER_DIR=$PROJ_DIR/to_align/kmers

CHM13_GRAPH=$GRAPH_DIR/chm13.chr12asat
OWN_HAP_GRAPH=$GRAPH_DIR/${ORIG_PATH_NAME}.own_hap
NEIGHBOR_GRAPH=$GRAPH_DIR/${ORIG_PATH_NAME}.neighbor
SAMPLED_GRAPH=$GRAPH_DIR/${ORIG_PATH_NAME}.sampled

OWN_HAP_ALN=$ALN_DIR/${ORIG_PATH_NAME}.own_hap
NEIGHBOR_ALN=$ALN_DIR/${ORIG_PATH_NAME}.neighbor
CHM13_ALN=$ALN_DIR/${ORIG_PATH_NAME}.chm13
SAMPLED_ALN=$ALN_DIR/${ORIG_PATH_NAME}.sampled

GUESS_LOG=$GRAPH_DIR/${ORIG_PATH_NAME}.guess.log

# ---- get reads to align ----

if [ ! -f ${REAL_READS}.fastq ]; then
    # Download reads
    echo "Downloading reads for $ORIG_PATH_NAME from AWS:"
    reads=`grep "^$SAMPLE_ID," $PROJ_DIR/to_align/aws_file_locations.csv | cut -f3 -d ","` 
    echo "$reads"
    aws s3 --no-sign-request cp $reads $PROJ_DIR/to_align/${ORIG_PATH_NAME}.bam &> /dev/null
    if [ ! -f "$PROJ_DIR/to_align/${ORIG_PATH_NAME}.bam" ]; then
        echo "ERROR: Could not find reads for $ORIG_PATH_NAME"
        exit 1
    fi
    grep chr12 ${BED_DIR}/${SAMPLE_NAME}_* > ${REAL_READS}.bed
    # Convert to FASTQ
    samtools view -@32 -L ${REAL_READS}.bed $PROJ_DIR/to_align/${ORIG_PATH_NAME}.bam | \
        awk '{print "@" $1 "\n" $10 "\n+\n" $11}' > ${REAL_READS}.fastq

    # Clean up memory
    rm -f $PROJ_DIR/to_align/${ORIG_PATH_NAME}*.bam
fi

if [ ! -f ${REAL_READS}.fastq ]; then
    echo "Failed to create a FASTQ file"
    exit 1
fi

if [ ! -s ${REAL_READS}.fastq ]; then
    echo "FASTQ file is empty"
    exit 1
fi

# ---- align to own haplotype ----

# Align to own haplotype
vg mod --keep-path "$PATH_NAME" ${BIG_GRAPH}.pg > ${OWN_HAP_GRAPH}.vg
# Convert to GBZ
vg gbwt --index-paths -x ${OWN_HAP_GRAPH}.vg -o ${OWN_HAP_GRAPH}.no_ref.gbwt
vg gbwt -o ${OWN_HAP_GRAPH}.with_ref.gbwt --set-reference "$SAMPLE_ID" ${OWN_HAP_GRAPH}.no_ref.gbwt
vg gbwt --gbz-format -x ${OWN_HAP_GRAPH}.vg ${OWN_HAP_GRAPH}.with_ref.gbwt -g ${OWN_HAP_GRAPH}.gbz
# Now that we have the GBZ we don't need these other things
rm ${OWN_HAP_GRAPH}.no_ref.gbwt ${OWN_HAP_GRAPH}.with_ref.gbwt ${OWN_HAP_GRAPH}.vg
# Index
vg autoindex --gbz ${OWN_HAP_GRAPH}.gbz -w lr-giraffe --prefix "$OWN_HAP_GRAPH" --no-guessing

# Minimap2 needs FASTA input
vg paths --extract-fasta -x ${OWN_HAP_GRAPH}.gbz > ${OWN_HAP_GRAPH}.fasta
# Avoid reusing an old index
rm -f "${OWN_HAP_GRAPH}.fasta.fai"
minimap2 -x lr:hqae -d ${OWN_HAP_GRAPH}.mmi ${OWN_HAP_GRAPH}.fasta

./helper_scripts/align_reads_minimap2.sh lr:hqae ${OWN_HAP_GRAPH}.mmi \
    ${OWN_HAP_GRAPH}.gbz ${REAL_READS}.fastq ${OWN_HAP_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${REAL_READS}.fastq ${OWN_HAP_ALN}.real.giraffe

# ---- align to CHM13 ----

./helper_scripts/align_reads_minimap2.sh map-hifi ${CHM13_GRAPH}.mmi \
    ${CHM13_GRAPH}.gbz ${REAL_READS}.fastq ${CHM13_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${REAL_READS}.fastq ${CHM13_ALN}.real.giraffe

# ---- align to nearest neighbor ----

# Get nearest neighbor
nearest_neighbor=`grep "$ORIG_PATH_NAME" "$DISTS" | sed 's/,/\t/g' | sort -k3 -n | head -1 \
    | cut -f1-2 | tr "\t" "\n" | grep -v "$ORIG_PATH_NAME"`
neighbor_sample_id=`echo $nearest_neighbor | cut -f1 -d "."`
neighbor_haplo_num=`echo $nearest_neighbor | cut -f2 -d "."`
neighbor_path_name="${neighbor_sample_id}#${neighbor_haplo_num}#${nearest_neighbor}#0"

echo "Aligning to nearest neighbor: $neighbor_path_name"

vg mod --keep-path "$neighbor_path_name" ${BIG_GRAPH}.pg > ${NEIGHBOR_GRAPH}.vg
# Convert to GBZ
vg gbwt --index-paths -x ${NEIGHBOR_GRAPH}.vg -o ${NEIGHBOR_GRAPH}.no_ref.gbwt
vg gbwt -o ${NEIGHBOR_GRAPH}.with_ref.gbwt --set-reference "$neighbor_sample_id" ${NEIGHBOR_GRAPH}.no_ref.gbwt
vg gbwt --gbz-format -x ${NEIGHBOR_GRAPH}.vg ${NEIGHBOR_GRAPH}.with_ref.gbwt -g ${NEIGHBOR_GRAPH}.gbz
# Now that we have the GBZ we don't need these other things
rm ${NEIGHBOR_GRAPH}.no_ref.gbwt ${NEIGHBOR_GRAPH}.with_ref.gbwt ${NEIGHBOR_GRAPH}.vg
# Index
vg autoindex --gbz ${NEIGHBOR_GRAPH}.gbz -w lr-giraffe --prefix "$NEIGHBOR_GRAPH" --no-guessing

# Minimap2 needs FASTA input
vg paths --extract-fasta -x ${NEIGHBOR_GRAPH}.gbz > ${NEIGHBOR_GRAPH}.fasta
# Avoid reusing an old index
rm -f "${NEIGHBOR_GRAPH}.fasta.fai"
minimap2 -x map-hifi -d ${NEIGHBOR_GRAPH}.mmi ${NEIGHBOR_GRAPH}.fasta

./helper_scripts/align_reads_minimap2.sh map-hifi ${NEIGHBOR_GRAPH}.mmi \
    ${NEIGHBOR_GRAPH}.gbz ${REAL_READS}.fastq ${NEIGHBOR_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${REAL_READS}.fastq ${NEIGHBOR_ALN}.real.giraffe

# ---- align to to haplotype-sampled graphs ----

echo "kmc -k29 -m128 -okff -t16 -hp ${REAL_READS}.fastq \
    $KMER_DIR/${ORIG_PATH_NAME}.real $KMER_DIR"
kmc -k29 -m128 -okff -t16 -hp ${REAL_READS}.fastq \
    $KMER_DIR/${ORIG_PATH_NAME}.real "$KMER_DIR"

echo "vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 10 --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g /dev/null --ban-sample $SAMPLE_ID ${BIG_GRAPH}.gbz 2> $GUESS_LOG"
# Sample 10 haps without alignment, so we can guess ideal # to sample
vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes 10 --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g /dev/null --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> "$GUESS_LOG"

# Use logfile to guess
./guess_n_and_cenhap.py --ploidy 1 --cenhap-table "$CENHAP_TABLE" "$GUESS_LOG" &>> "$GUESS_LOG"
cat "$GUESS_LOG"

n_to_sample=`fgrep Best "$GUESS_LOG" | cut -d " " -f4`

vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${SAMPLED_GRAPH}.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix "$SAMPLED_GRAPH" --no-guessing --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.gbz

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.gbz ${REAL_READS}.fastq $SAMPLED_ALN.real.giraffe