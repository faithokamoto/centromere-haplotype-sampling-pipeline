#!/bin/bash
# Run typing with different --absent-score for chr4

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BIG_GRAPH=$PROJ_DIR/graph/unsampled/chr4

for hap_name in log/chr4/*; do 
    hap_name=`echo "$file" | cut -f2-3 -d "."`
    sample_id=`echo "$hap_name" | cut -f1 -d "." `

    prefix=chr4.$hap_name
    reads=$PROJ_DIR/to_align/$prefix
    kmer_dir=$PROJ_DIR/kmers/$prefix

    rm -rf $kmer_dir
    mkdir -p $kmer_dir

    kmc -k29 -m128 -okff -t16 -hp ${reads}.real.fastq.gz \
        $kmer_dir/${prefix}.real "$kmer_dir"

    for discount in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 30 40 50 60 70 80; do
        vg haplotypes -k $kmer_dir/${prefix}.real.kff -i ${BIG_GRAPH}.hapl \
            --num-haplotypes 1 --haploid-scoring -d ${BIG_GRAPH}.dist --absent-score 0.$discount \
            -g /dev/null --ban-sample "$sample_id" ${BIG_GRAPH}.gbz 2> log/typing_tests/$prefix.$discount.real.log
    done
done