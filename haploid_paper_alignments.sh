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
SCRATCH=/scratch/fokamoto
BIG_GRAPH=$PROJ_DIR/graph/unsampled/$CHROM
BAM_LOCS=$PROJ_DIR/to_align/aws_file_locations.csv
MIRA_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2
BED_DIR=$MIRA_DIR/per_smp_asat_beds
DISTS=$MIRA_DIR/all_pairs/distance_matrices/${CHROM}_r2_QC_v2_centrolign_pairwise_distance.csv
CENHAP_TABLE=/private/groups/migalab/juklucas/centrolign/cenhap_assignment/cenhap_inference_out/${CHROM}/${CHROM}.cenhap_predictions.tsv

READS=$PROJ_DIR/to_align/${ORIG_PATH_NAME}.${CHROM}.hifi

GRAPH_DIR=$PROJ_DIR/graph/haploid
ALN_DIR=$PROJ_DIR/alignments/haploid
KMER_DIR=$PROJ_DIR/to_align/kmers/$ORIG_PATH_NAME

mkdir -p $KMER_DIR

CHM13_GRAPH=$GRAPH_DIR/chm13.${CHROM}asat
OWN_HAP_GRAPH=$GRAPH_DIR/${ORIG_PATH_NAME}.${CHROM}.own_hap
NEIGHBOR_GRAPH=$GRAPH_DIR/${ORIG_PATH_NAME}.${CHROM}.neighbor
SAMPLED_GRAPH=$GRAPH_DIR/${ORIG_PATH_NAME}.${CHROM}.sampled

OWN_HAP_ALN=$ALN_DIR/${ORIG_PATH_NAME}.${CHROM}.own_hap
NEIGHBOR_ALN=$ALN_DIR/${ORIG_PATH_NAME}.${CHROM}.neighbor
CHM13_ALN=$ALN_DIR/${ORIG_PATH_NAME}.${CHROM}.chm13
SAMPLED_ALN=$ALN_DIR/${ORIG_PATH_NAME}.${CHROM}.sampled

GUESS_LOG=$GRAPH_DIR/${ORIG_PATH_NAME}.${CHROM}.guess

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

if [ ! -f ${READS}.real.fastq ]; then
    # Download reads
    echo "Downloading reads for $ORIG_PATH_NAME from AWS"
    if [ `grep -L "^$SAMPLE_ID," "$BAM_LOCS" | wc -l` -eq 1 ]; then
        echo "ERROR: Could not find reads for $ORIG_PATH_NAME"
        exit 1
    fi
    reads=`grep "^$SAMPLE_ID," "$BAM_LOCS" | cut -f3 -d ","` 
    full_bam=$SCRATCH/${ORIG_PATH_NAME}.bam
    aws s3 --no-sign-request cp "$reads" "$full_bam" &> /dev/null
    echo "Download complete"
    # Subset BAM to only correct-chromosome reads
    grep ${CHROM} ${BED_DIR}/${ORIG_PATH_NAME}_asat_arrays.bed > ${READS}.real.bed
    samtools view -@32 -L ${READS}.real.bed -h "$full_bam" > ${READS}.real.sam
    # Get rid of giant BAM file for space
    rm $SCRATCH/${ORIG_PATH_NAME}.*bam*
    echo "Cleaned up BAM file"

    # Edit SAM to something compatible with the graph
    old_path_name=`cut -f1 ${READS}.real.bed`
    start=`cut -f2 ${READS}.real.bed`
    end=`cut -f3 ${READS}.real.bed`
    new_path_name=`echo $PATH_NAME | sed 's/#0//g'`
    samtools view ${READS}.real.sam | sed "s/$old_path_name/$new_path_name/" > ${READS}.real.no_header.sam
    samtools view -H ${READS}.real.sam | sed "s/$old_path_name/$new_path_name/" > ${READS}.real.header
    # Filter for reads which appear within the BED file's boundaries
    # while also editing the coordinates to be graph-friendly
    ./helper_scripts/edit_sam.py --start "$start" --end "$end" \
        ${READS}.real.no_header.sam > ${READS}.real.edited.sam
    cat ${READS}.real.header ${READS}.real.edited.sam > ${READS}.real.combined.sam
    
    # Get truth positions
    vg inject -x ${OWN_HAP_GRAPH}.gbz ${READS}.real.combined.sam > ${READS}.real.gam
    vg filter --tsv-out "name;nodes" ${READS}.real.gam > ${READS}.real.tsv
    # Convert to FASTQ
    vg view --fastq-out ${READS}.real.gam > ${READS}.real.fastq

    # Clean up memory; we only need FASTQ files & the nodes TSV
    rm -f ${READS}.real.bed ${READS}.real.sam ${READS}.real.no_header.sam ${READS}.real.header
    rm -f ${READS}.real.edited.sam ${READS}.real.combined.sam ${READS}.real.gam
fi

if [ ! -f ${READS}.real.fastq ]; then
    echo "Failed to create a FASTQ file"
    exit 1
fi

if [ ! -s ${READS}.real.fastq ]; then
    echo "FASTQ file is empty"
    exit 1
fi

# Generate similar simulated reads
num_reads=`grep -c "^@" ${READS}.real.fastq`
vg sim -x ${OWN_HAP_GRAPH}.gbz --num-reads "$num_reads" \
     --random-seed 42 --threads 20 --fastq ${READS}.real.fastq \
     --align-out --use-average-length > ${READS}.sim.gam
vg filter --tsv-out "name;nodes" ${READS}.sim.gam > ${READS}.sim.tsv
# Convert to FASTQ
vg view --fastq-out ${READS}.sim.gam > ${READS}.sim.fastq

# Clean up memory; we only need FASTQ files & the nodes TSV
rm -f ${READS}.sim.gam

# ---- basic alignments (real reads) ----

echo "====="
echo "Linear reference alignments for real reads"

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh ${OWN_HAP_GRAPH}.mmi \
    ${OWN_HAP_GRAPH}.gbz ${READS}.real.fastq ${OWN_HAP_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${READS}.real.fastq ${OWN_HAP_ALN}.real.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh ${CHM13_GRAPH}.mmi \
    ${CHM13_GRAPH}.gbz ${READS}.real.fastq ${CHM13_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${READS}.real.fastq ${CHM13_ALN}.real.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh ${NEIGHBOR_GRAPH}.mmi \
    ${NEIGHBOR_GRAPH}.gbz ${READS}.real.fastq ${NEIGHBOR_ALN}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${NEIGHBOR_GRAPH}.gbz \
    ${READS}.real.fastq ${NEIGHBOR_ALN}.real.giraffe

# ---- basic alignments (sim reads) ----

echo "====="
echo "Linear reference alignments for sim reads"

# align to own haplotype
./helper_scripts/align_reads_minimap2.sh ${OWN_HAP_GRAPH}.mmi \
    ${OWN_HAP_GRAPH}.gbz ${READS}.sim.fastq ${OWN_HAP_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.gbz \
    ${READS}.sim.fastq ${OWN_HAP_ALN}.sim.giraffe

# align to CHM13
./helper_scripts/align_reads_minimap2.sh ${CHM13_GRAPH}.mmi \
    ${CHM13_GRAPH}.gbz ${READS}.sim.fastq ${CHM13_ALN}.sim.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.gbz \
    ${READS}.sim.fastq ${CHM13_ALN}.sim.giraffe

# align to nearest neighbor
./helper_scripts/align_reads_minimap2.sh ${NEIGHBOR_GRAPH}.mmi \
    ${NEIGHBOR_GRAPH}.gbz ${READS}.sim.fastq ${NEIGHBOR_ALN}.sim.minimap2
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
        --dist-matrix $DISTS ${GUESS_LOG}.real.log &>> ${GUESS_LOG}.real.log
cat ${GUESS_LOG}.real.log

n_to_sample=`fgrep Best ${GUESS_LOG}.real.log | cut -d " " -f4`

vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.real.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${SAMPLED_GRAPH}.real.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix $SAMPLED_GRAPH.real --no-guessing \
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
        --dist-matrix $DISTS ${GUESS_LOG}.sim.log &>> ${GUESS_LOG}.sim.log
cat ${GUESS_LOG}.sim.log

n_to_sample=`fgrep Best ${GUESS_LOG}.sim.log | cut -d " " -f4`

vg haplotypes -k $KMER_DIR/${ORIG_PATH_NAME}.${CHROM}.sim.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes "$n_to_sample" --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${SAMPLED_GRAPH}.sim.gbz --ban-sample "$SAMPLE_ID" ${BIG_GRAPH}.gbz 2> /dev/null
vg autoindex --prefix $SAMPLED_GRAPH.sim --no-guessing \
    --workflow lr-giraffe --gbz ${SAMPLED_GRAPH}.sim.gbz 2> /dev/null

./helper_scripts/align_reads_giraffe.sh ${SAMPLED_GRAPH}.sim.gbz ${READS}.sim.fastq ${SAMPLED_ALN}.sim.giraffe