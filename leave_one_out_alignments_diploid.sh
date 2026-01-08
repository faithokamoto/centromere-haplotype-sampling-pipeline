#!/bin/bash
# Run haplotype sampling with a leave-one-out graph on chr12
# Usage: leave_one_out_alignments_diploid.sh <sample id>
# Example: leave_one_out_alignments.sh HG00099

# ---- set up variables ----

SAMPLE_ID=$1

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr12

CHM13_GRAPH=$PROJ_DIR/graph/linear_refs/chm13.chr12asat
OWN_HAP_GRAPH=$PROJ_DIR/graph/linear_refs/$SAMPLE_ID

REAL_READS=$PROJ_DIR/to_align/real_${SAMPLE_ID}_chr12_hor_array.hifi

KMER_DIR=$PROJ_DIR/to_align/kmers

OWN_HAP_PREFIX=$PROJ_DIR/alignments/linear_refs/${SAMPLE_ID}.own_hap
CHM13_PREFIX=$PROJ_DIR/alignments/linear_refs/${SAMPLE_ID}.chm13

BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files

echo "Processing sample: $SAMPLE_ID"

# ---- get reads to align ----

if [ ! -f ${REAL_READS}.fastq ]; then
    # Download reads
    echo "Downloading reads for $SAMPLE_ID from AWS:"
    reads=`grep "^$SAMPLE_ID," $PROJ_DIR/to_align/aws_file_locations.csv | cut -f3 -d ","` 
    echo "$reads"
    aws s3 --no-sign-request cp $reads $PROJ_DIR/to_align/${SAMPLE_ID}.bam &> /dev/null
    if [ ! -f "$PROJ_DIR/to_align/${SAMPLE_ID}.bam" ]; then
        echo "ERROR: Could not find reads for $SAMPLE_ID"
        exit 1
    fi
    grep chr12 ${BED_DIR}/${SAMPLE_ID}* | grep -e hap1 -e pat |cut -f2 -d ":" > ${REAL_READS}.1.bed
    grep chr12 ${BED_DIR}/${SAMPLE_ID}* | grep -e hap2 -e mat |cut -f2 -d ":" > ${REAL_READS}.2.bed
    # Convert to FASTQ
    samtools view -@32 -L ${REAL_READS}.1.bed $PROJ_DIR/to_align/${SAMPLE_ID}.bam | \
        awk '{print "@" $1 "\n" $10 "\n+\n" $11}' > ${REAL_READS}.1.fastq
    samtools view -@32 -L ${REAL_READS}.2.bed $PROJ_DIR/to_align/${SAMPLE_ID}.bam | \
        awk '{print "@" $1 "\n" $10 "\n+\n" $11}' > ${REAL_READS}.2.fastq
    cat ${REAL_READS}.1.fastq ${REAL_READS}.2.fastq > ${REAL_READS}.fastq
    
    # Clean up memory
    rm $PROJ_DIR/to_align/${SAMPLE_ID}.bam
fi

# ---- align to own haplotypes ----

# Align to own haplotypes
vg paths --paths-by $SAMPLE_ID --retain-paths -x $BIG_GRAPH.pg > ${OWN_HAP_GRAPH}.vg
vg paths --extract-fasta -x ${OWN_HAP_GRAPH}.vg > ${OWN_HAP_GRAPH}.fasta
# Avoid reusing an old index
rm -f ${OWN_HAP_GRAPH}.fasta.fai
# Convert to GBZ
vg gbwt --index-paths -x ${OWN_HAP_GRAPH}.vg -o ${OWN_HAP_GRAPH}.gbwt
vg gbwt --gbz-format -x ${OWN_HAP_GRAPH}.vg ${OWN_HAP_GRAPH}.gbwt \
    -g ${OWN_HAP_GRAPH}.giraffe.gbz
# Index
vg autoindex --gbz ${OWN_HAP_GRAPH}.giraffe.gbz -w lr-giraffe \
    --prefix $OWN_HAP_GRAPH --no-guessing
minimap2 -x map-hifi -d ${OWN_HAP_GRAPH}.mmi ${OWN_HAP_GRAPH}.fasta

./helper_scripts/align_reads_minimap2.sh ${OWN_HAP_GRAPH}.mmi ${OWN_HAP_GRAPH}.giraffe.gbz \
    ${REAL_READS}.fastq ${OWN_HAP_PREFIX}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${OWN_HAP_GRAPH}.giraffe.gbz \
    ${REAL_READS}.fastq ${OWN_HAP_PREFIX}.real.giraffe

# ---- align to CHM13 ----

./helper_scripts/align_reads_minimap2.sh ${CHM13_GRAPH}.mmi ${CHM13_GRAPH}.giraffe.gbz \
    ${REAL_READS}.fastq ${CHM13_PREFIX}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${CHM13_GRAPH}.giraffe.gbz \
    ${REAL_READS}.fastq ${CHM13_PREFIX}.real.giraffe

# ---- align to to haplotype-sampled graphs ----

kmc -k29 -m128 -okff -t16 -hp ${REAL_READS}.fastq \
    $KMER_DIR/real_${SAMPLE_ID} $KMER_DIR > $KMER_DIR/real_${SAMPLE_ID}.kff.log

for num_hap in {1..10}
do
    echo "Sampling $num_hap haplotypes"
    hap_prefix=${SAMPLE_ID}.${num_hap}haps
    real_graph=$PROJ_DIR/graph/leave_one_out/real_${hap_prefix}
    real_out=$PROJ_DIR/alignments/leave_one_out/real_${hap_prefix}

    # Haplotype sampling on leave-one-out graph
    vg haplotypes -t 1 -k $KMER_DIR/real_${SAMPLE_ID}.kff -i ${BIG_GRAPH}.hapl \
        --num-haplotypes $num_hap --haploid-scoring -d ${BIG_GRAPH}.dist \
        -g ${real_graph}.giraffe.gbz --ban-sample $SAMPLE_ID ${BIG_GRAPH}.gbz
    vg autoindex --prefix $real_graph --no-guessing \
        --workflow lr-giraffe --gbz ${real_graph}.giraffe.gbz 2> /dev/null

    # Align reads to sampled graph
    ./helper_scripts/align_reads_giraffe.sh ${real_graph}.giraffe.gbz ${REAL_READS}.fastq $real_out
done

./guess_optimal_num_sampled_haplo.py \
    --hap1-reads ${REAL_READS}.1.fastq \
    --hap2-reads ${REAL_READS}.2.fastq \
    $SAMPLE_ID