#!/bin/bash
# Run haplotype sampling with a leave-one-out graph on chr12
# Usage: haploid_paper_alignments.sh <original path name>
# Example: haploid_paper_alignments.sh HG00099.1

# ---- process arguments ----

echo "Top of haploid_paper_alignments.sh with ${1} input"

ORIG_PATH_NAME=$1
SAMPLE_ID=`echo "$ORIG_PATH_NAME" | cut -f1 -d "." `
HAPLO_NUM=`echo "$ORIG_PATH_NAME" | cut -f2 -d "." `
PATH_NAME="${SAMPLE_ID}#${HAPLO_NUM}#${ORIG_PATH_NAME}#0"

BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files

ls ${BED_DIR}/${SAMPLE_ID}_* &>/dev/null
if [ $? -ne 0 ]; then
    echo "No BED files available for $SAMPLE_ID"
    exit 1
fi

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

set -e

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr12
DISTS=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=/private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/chr12/chr12.cenhap_predictions.tsv

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

# ---- make extra references ----

# native ref
./helper_scripts/create_single_path_ref.sh ${BIG_GRAPH}.pg "$PATH_NAME" "$OWN_HAP_GRAPH"

# Nearest neighbor
nearest_neighbor=`grep "$ORIG_PATH_NAME" "$DISTS" | sed 's/,/\t/g' | sort -k3 -n | head -1 \
    | cut -f1-2 | tr "\t" "\n" | grep -v "$ORIG_PATH_NAME"`
neighbor_sample_id=`echo $nearest_neighbor | cut -f1 -d "."`
neighbor_haplo_num=`echo $nearest_neighbor | cut -f2 -d "."`
neighbor_path_name="${neighbor_sample_id}#${neighbor_haplo_num}#${nearest_neighbor}#0"
echo "Nearest neighbor: $neighbor_path_name"

./helper_scripts/create_single_path_ref.sh ${BIG_GRAPH}.pg "$neighbor_path_name" "$NEIGHBOR_GRAPH"

# ---- get reads to align ----

if [ ! -f ${REAL_READS}.fastq ]; then
    # Download reads
    echo "Downloading reads for $ORIG_PATH_NAME from AWS"
    reads=`grep "^$SAMPLE_ID," $PROJ_DIR/to_align/aws_file_locations.csv | cut -f3 -d ","` 
    full_bam=$PROJ_DIR/to_align/${ORIG_PATH_NAME}.bam
    aws s3 --no-sign-request cp "$reads" "$full_bam" &> /dev/null
    echo "Download complete"
    if [ ! -f "$full_bam" ]; then
        echo "ERROR: Could not find reads for $ORIG_PATH_NAME"
        exit 1
    fi
    # Subset BAM to only chr12 reads
    grep chr12 ${BED_DIR}/${SAMPLE_NAME}_* > ${REAL_READS}.bed
    samtools view -@32 -L ${REAL_READS}.bed -h "$full_bam" > ${REAL_READS}.sam
    # Get rid of giant BAM file for space
    rm $PROJ_DIR/to_align/${ORIG_PATH_NAME}.*bam*
    echo "Cleaned up BAM file"

    # Edit SAM to something compatible with the graph
    old_path_name=`cut -f1 ${REAL_READS}.bed`
    start=`cut -f2 ${REAL_READS}.bed`
    end=`cut -f3 ${REAL_READS}.bed`
    new_path_name=`echo $PATH_NAME | sed 's/#0//g'`
    samtools view ${REAL_READS}.sam > ${REAL_READS}.no_header.sam
    samtools view -H ${REAL_READS}.sam | sed "s/$old_path_name/$new_path_name/"> ${REAL_READS}.header
    # Filter for reads which appear within the BED file's boundaries
    # while also editing the coordinates to be graph-friendly
    awk -v contig="$new_path_name" -v start="$start" -v end="$end" \
        '{if (($4 - start + 1 >= 0) && ($4 + length($10) <= end)) 
            {print $1, $2, contig, $4 - start + 1, $5, $6, $7, $8, $9, $10, $11}}' \
        ${REAL_READS}.no_header.sam | tr " " "\t" > ${REAL_READS}.edited.sam
    cat ${REAL_READS}.header ${REAL_READS}.edited.sam > ${REAL_READS}.combined.sam
    
    # Get truth positions
    vg inject -x ${OWN_HAP_GRAPH}.gbz ${REAL_READS}.combined.sam > ${REAL_READS}.gam
    vg filter --tsv-out "name;nodes" ${REAL_READS}.gam > ${REAL_READS}.tsv
    # Convert to FASTQ
    vg view --fastq-out ${REAL_READS}.gam > ${REAL_READS}.fastq

    # Clean up memory; we only need FASTQ files & the nodes TSV
    rm -f ${REAL_READS}.bed ${REAL_READS}.sam ${REAL_READS}.no_header.sam ${REAL_READS}.header
    rm -f ${REAL_READS}.edited.sam ${REAL_READS}.combined.sam ${REAL_READS}.gam
fi

if [ ! -f ${REAL_READS}.fastq ]; then
    echo "Failed to create a FASTQ file"
    exit 1
fi

if [ ! -s ${REAL_READS}.fastq ]; then
    echo "FASTQ file is empty"
    exit 1
fi

# ---- basic alignments ----

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh ${OWN_HAP_GRAPH}.mmi \
    ${OWN_HAP_GRAPH}.gbz ${REAL_READS}.fastq ${OWN_HAP_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${REAL_READS}.fastq ${OWN_HAP_ALN}.real.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh ${CHM13_GRAPH}.mmi \
    ${CHM13_GRAPH}.gbz ${REAL_READS}.fastq ${CHM13_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${REAL_READS}.fastq ${CHM13_ALN}.real.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh ${NEIGHBOR_GRAPH}.mmi \
    ${NEIGHBOR_GRAPH}.gbz ${REAL_READS}.fastq ${NEIGHBOR_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${REAL_READS}.fastq ${NEIGHBOR_ALN}.real.giraffe

# ---- align to to haplotype-sampled graphs ----

kmc -k29 -m128 -okff -t16 -hp ${REAL_READS}.fastq \
    $KMER_DIR/${ORIG_PATH_NAME}.real "$KMER_DIR"

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