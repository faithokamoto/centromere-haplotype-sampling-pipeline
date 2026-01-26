#!/bin/bash
# Run haplotype sampling with a leave-one-out graph on chr12, only 10 haps
# Usage: sample_for_cenhap_guess.sh <path name>
# Example: sample_for_cenhap_guess.sh HG00099.1

# ---- set up variables ----

PATH_NAME=$1
SAMPLE_ID=`echo "$PATH_NAME" | cut -f1 -d "." `
HAPLO_NUM=`echo "$PATH_NAME" | cut -f2 -d "." `

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr12

READ_SUFFIX=chr12_hor_array.hifi
REAL_READS=$PROJ_DIR/to_align/real_${PATH_NAME}_$READ_SUFFIX

BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/extract_hors_HPRC/release2/contiguous_HORs_bed_files
KMER_DIR=$PROJ_DIR/to_align/kmers

echo "Processing sample: $PATH_NAME"

# ---- get reads to align ----

if [ ! -f ${REAL_READS}.fastq ]; then
    # Download reads
    echo "Downloading reads for $PATH_NAME from AWS:"
    reads=`grep "^$SAMPLE_ID," $PROJ_DIR/to_align/aws_file_locations.csv | cut -f3 -d ","` 
    echo "$reads"
    aws s3 --no-sign-request cp $reads $PROJ_DIR/to_align/${SAMPLE_ID}.bam &> /dev/null
    if [ ! -f "$PROJ_DIR/to_align/${SAMPLE_ID}.bam" ]; then
        echo "ERROR: Could not find reads for $PATH_NAME"
        exit 1
    fi
    if [ "$HAPLO_NUM" -eq "1" ]; then
        grep chr12 ${BED_DIR}/${SAMPLE_ID}* | grep -e hap1 -e pat |cut -f2 -d ":" > ${REAL_READS}.bed
    else
        grep chr12 ${BED_DIR}/${SAMPLE_ID}* | grep -e hap2 -e mat |cut -f2 -d ":" > ${REAL_READS}.bed
    fi
    # Convert to FASTQ
    samtools view -@32 -L ${REAL_READS}.bed $PROJ_DIR/to_align/${SAMPLE_ID}.bam | \
        awk '{print "@" $1 "\n" $10 "\n+\n" $11}' > ${REAL_READS}.fastq
    
    # Clean up memory
    rm $PROJ_DIR/to_align/${SAMPLE_ID}.bam
fi

# ---- align to to haplotype-sampled graphs ----

kmc -k29 -m128 -okff -t16 -hp ${REAL_READS}.fastq \
    $KMER_DIR/real_${PATH_NAME} $KMER_DIR > $KMER_DIR/real_${PATH_NAME}.kff.log

num_hap=10
hap_prefix=${PATH_NAME}.${num_hap}haps
real_graph=$PROJ_DIR/graph/leave_one_out/real_${hap_prefix}

echo "Sampling ${num_hap} haplotypes"
# Haplotype sampling on leave-one-out graph (only 10 haps)
vg haplotypes -t 1 -k $KMER_DIR/real_${PATH_NAME}.kff -i ${BIG_GRAPH}.hapl \
    --num-haplotypes ${num_hap} --haploid-scoring -d ${BIG_GRAPH}.dist \
    -g ${real_graph}.giraffe.gbz --ban-sample $SAMPLE_ID ${BIG_GRAPH}.gbz