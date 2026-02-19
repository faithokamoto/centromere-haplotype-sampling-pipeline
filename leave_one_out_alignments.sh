#!/bin/bash
# Run haplotype sampling with a leave-one-out graph on chr12
# Usage: leave_one_out_alignments.sh <original path name> <version>
# Example: leave_one_out_alignments.sh HG00099.1 v1

# ---- process arguments ----

ORIG_PATH_NAME=$1
VERSION=$2
SAMPLE_ID=`echo "$ORIG_PATH_NAME" | cut -f1 -d "." `
HAPLO_NUM=`echo "$ORIG_PATH_NAME" | cut -f2 -d "." `
PATH_NAME="${SAMPLE_ID}#${HAPLO_NUM}#${ORIG_PATH_NAME}#0"

BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files
BED_SUFFIX=hprc_r2_${VERSION}_hor_arrays.bed

# Figure out what the sample name is
if [ ! -f "${BED_DIR}/${SAMPLE_ID}_hap1_${BED_SUFFIX}" ]; then
    if [ "$HAPLO_NUM" -eq "1" ]; then
        SAMPLE_NAME=${SAMPLE_ID}_pat
    else
        SAMPLE_NAME=${SAMPLE_ID}_mat
    fi
else
    if [ "$HAPLO_NUM" -eq "1" ]; then
        SAMPLE_NAME=${SAMPLE_ID}_hap1
    else
        SAMPLE_NAME=${SAMPLE_ID}_hap2
    fi
fi

echo "Processing sample: $SAMPLE_NAME"

# ---- set up variables ----

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr12
DISTS=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/all_pairs/distance_matrices/chr12_r2_QC_v2_centrolign_pairwise_distance.csv
TRUTH_CSV_DIR=/private/groups/migalab/juklucas/centrolign/variant_calling/rates/chr12/snv_calls/filtered/

CHM13_GRAPH=$PROJ_DIR/graph/linear_refs/chm13.chr12asat
OWN_HAP_GRAPH=$PROJ_DIR/graph/linear_refs/$SAMPLE_NAME

READ_SUFFIX=chr12_hor_array.hifi
REAL_READS=$PROJ_DIR/to_align/real_${SAMPLE_NAME}_$READ_SUFFIX

KMER_DIR=$PROJ_DIR/to_align/kmers

OWN_HAP_PREFIX=$PROJ_DIR/alignments/linear_refs/${SAMPLE_NAME}.own_hap
NEIGHBOR_PREFIX=$PROJ_DIR/alignments/linear_refs/${SAMPLE_NAME}.neighbor
CHM13_PREFIX=$PROJ_DIR/alignments/linear_refs/${SAMPLE_NAME}.chm13

# ---- get reads to align ----

if [ ! -f ${REAL_READS}.fastq ]; then
    # Download reads
    echo "Downloading reads for $SAMPLE_NAME from AWS:"
    reads=`grep "^$SAMPLE_ID," $PROJ_DIR/to_align/aws_file_locations.csv | cut -f3 -d ","` 
    echo "$reads"
    aws s3 --no-sign-request cp $reads $PROJ_DIR/to_align/${SAMPLE_NAME}.bam &> /dev/null
    if [ ! -f "$PROJ_DIR/to_align/${SAMPLE_NAME}.bam" ]; then
        echo "ERROR: Could not find reads for $SAMPLE_NAME"
        exit 1
    fi
    grep chr12 ${BED_DIR}/${SAMPLE_NAME}_${BED_SUFFIX} > ${REAL_READS}.bed
    # Convert to FASTQ
    samtools view -@32 -L ${REAL_READS}.bed $PROJ_DIR/to_align/${SAMPLE_NAME}.bam | \
        awk '{print "@" $1 "\n" $10 "\n+\n" $11}' > ${REAL_READS}.fastq

    # Clean up memory
    rm $PROJ_DIR/to_align/${SAMPLE_NAME}.bam
fi

# ---- align to own haplotype ----

# Align to own haplotype
vg paths --paths-by $PATH_NAME --extract-fasta -x $BIG_GRAPH.gbz > ${OWN_HAP_GRAPH}.fasta
# Avoid auto-conversion of name
sed "s/${PATH_NAME}/${ORIG_PATH_NAME}/" -i ${OWN_HAP_GRAPH}.fasta
# Avoid reusing an old index
rm ${OWN_HAP_GRAPH}.fasta.fai
# Convert to GBZ
vg construct --reference ${OWN_HAP_GRAPH}.fasta -m 1024 > ${OWN_HAP_GRAPH}.vg
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

# ---- align to nearest neighbor ----

# Get nearest neighbor
nearest_neighbor_line=`grep $ORIG_PATH_NAME $DISTS | sed 's/,/\t/g' | sort -k3 -n | head -1`
echo "Nearest neighbor line: $nearest_neighbor_line"
nearest_neighbor=`echo $nearest_neighbor_line | cut -f1-2 -d " " | tr " " "\n" | grep -v $ORIG_PATH_NAME`
neighbor_sample_id=`echo $nearest_neighbor | cut -f1 -d "."`
neighbor_haplo_num=`echo $nearest_neighbor | cut -f2 -d "."`
neighbor_path_name="${neighbor_sample_id}#${neighbor_haplo_num}#${nearest_neighbor}#0"
echo "Aligning to nearest neighbor: $neighbor_path_name"
neighbor_graph=$PROJ_DIR/graph/linear_refs/$nearest_neighbor

vg paths --paths-by $neighbor_path_name --extract-fasta -x $BIG_GRAPH.gbz > ${neighbor_graph}.fasta
# Avoid auto-conversion of name
sed "s/${neighbor_path_name}/${nearest_neighbor}/" -i ${neighbor_graph}.fasta
# Avoid reusing an old index
rm -f ${neighbor_graph}.fasta.fai
# Convert to GBZ
vg construct --reference ${neighbor_graph}.fasta -m 1024 > ${neighbor_graph}.vg
vg gbwt --index-paths -x ${neighbor_graph}.vg -o ${neighbor_graph}.gbwt
vg gbwt --gbz-format -x ${neighbor_graph}.vg ${neighbor_graph}.gbwt \
    -g ${neighbor_graph}.giraffe.gbz
# Index
vg autoindex --gbz ${neighbor_graph}.giraffe.gbz -w lr-giraffe \
    --prefix $neighbor_graph --no-guessing
minimap2 -x map-hifi -d ${neighbor_graph}.mmi ${neighbor_graph}.fasta

./helper_scripts/align_reads_minimap2.sh ${neighbor_graph}.mmi ${neighbor_graph}.giraffe.gbz \
    ${REAL_READS}.fastq ${NEIGHBOR_PREFIX}.real.minimap2
./helper_scripts/align_reads_giraffe.sh ${neighbor_graph}.giraffe.gbz \
    ${REAL_READS}.fastq ${NEIGHBOR_PREFIX}.real.giraffe

# ---- align to to haplotype-sampled graphs ----

kmc -k29 -m128 -okff -t16 -hp ${REAL_READS}.fastq \
    $KMER_DIR/real_${SAMPLE_NAME} $KMER_DIR > $KMER_DIR/real_${SAMPLE_NAME}.kff.log

for num_hap in {1..8}
do
    echo "Sampling $num_hap haplotypes"
    hap_prefix=${SAMPLE_NAME}.${num_hap}haps
    real_graph=$PROJ_DIR/graph/leave_one_out/real_${hap_prefix}
    real_out=$PROJ_DIR/alignments/leave_one_out/real_${hap_prefix}

    # Haplotype sampling on leave-one-out graph
    vg haplotypes -t 1 -k $KMER_DIR/real_${SAMPLE_NAME}.kff -i ${BIG_GRAPH}.hapl \
        --num-haplotypes $num_hap --haploid-scoring -d ${BIG_GRAPH}.dist \
        -g ${real_graph}.giraffe.gbz --ban-sample $SAMPLE_ID \
        --set-reference CHM13 --include-reference ${BIG_GRAPH}.gbz
    vg autoindex --prefix $real_graph --no-guessing \
        --workflow lr-giraffe --gbz ${real_graph}.giraffe.gbz

    # Align reads to sampled graph
    ./helper_scripts/align_reads_giraffe.sh ${real_graph}.giraffe.gbz ${REAL_READS}.fastq $real_out

    # Call variants relative to CHM13
    vg convert ${real_graph}.giraffe.gbz -p > ${real_graph}.pg
    vg augment ${real_graph}.pg ${real_out}.gam -m 3 -q 5 -Q 5 \
        -A ${real_out}.augment.gam > ${real_out}.augment.pg
    vg snarls ${real_out}.augment.pg > ${real_out}.augment.snarls
    vg pack -x ${real_out}.augment.pg -g ${real_out}.augment.gam -o ${real_out}.augment.pack
    vg call ${real_out}.augment.pg -r ${real_out}.augment.snarls --ploidy 1 \
        -k ${real_out}.augment.pack -s $SAMPLE_NAME --ref-sample augref_CHM13 > ${real_out}.vcf

    if [ -f $TRUTH_CSV_DIR/CHM13.0_${ORIG_PATH_NAME}.snvs.500bp_95pct.csv ]; then
        echo "Comparing SNVs to truth for $num_hap haplotypes"
        ./compare_snvs.py \
            --vcf ${real_out}.vcf \
            --truth-csv $TRUTH_CSV_DIR/CHM13.0_${ORIG_PATH_NAME}.snvs.500bp_95pct.csv \
            --relaxed-truth-csv $TRUTH_CSV_DIR/CHM13.0_${ORIG_PATH_NAME}.snvs.10bp_95pct.csv
    else
        echo "WARNING: Could not find truth CSV for $ORIG_PATH_NAME, skipping comparison"
    fi
done

./guess_optimal_num_sampled_haplo.py $SAMPLE_NAME