#!/bin/bash
# Run typing with different --absent-score for chr4
# Usage: testing_absent_score.sh

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr4

READS_DIR=/private/groups/patenlab/fokamoto/centrolign/to_align
LOG_DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/typing_tests
for file in $READS_DIR/chr4.*.real.fastq.gz; do 
    hap_name=`echo "$file" | cut -f8 -d "/" | cut -f2-3 -d "."`
    sample_id=`echo "$hap_name" | cut -f1 -d "." `

    prefix=chr4.$hap_name
    reads=$PROJ_DIR/to_align/$prefix
    kmer_dir=$TMPDIR/kmers/$prefix

    rm -rf $kmer_dir
    mkdir -p $kmer_dir

    kmc -k29 -m128 -okff -t16 -hp ${reads}.real.fastq.gz \
        $kmer_dir/$prefix "$kmer_dir"

    for discount in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 30 40 50 60 70 80; do
        vg haplotypes -k $kmer_dir/${prefix}.kff -i ${BIG_GRAPH}.hapl \
            --num-haplotypes 1 --haploid-scoring -d ${BIG_GRAPH}.dist --absent-score 0.$discount \
            -g /dev/null --ban-sample "$sample_id" ${BIG_GRAPH}.gbz 2> $LOG_DIR/$prefix.$discount.log
    done
done