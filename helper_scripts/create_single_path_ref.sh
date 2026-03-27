#!/bin/bash
# Extract a single path from larger graph and index it for alignment
# Usage: create_single_path_ref.sh <big graph> <path name> <small graph prefix>
# Example: create_single_path_ref.sh $PROJ_DIR/graph/unsampled/chr12.pg CHM13#0#CHM13.0#0 $PROJ_DIR/graph/haploid/chr12.CHM13

BIG_GRAPH=$1
PATH_NAME=$2
OUT=$3

SAMPLE_ID=`echo "$PATH_NAME" | cut -f1 -d "#"`

# Only bother if needed
if [ ! -f ${OUT}.mmi ]; then
    vg mod --keep-path "$PATH_NAME" "$BIG_GRAPH" > ${OUT}.vg
    # Convert to GBZ
    vg gbwt --index-paths -x ${OUT}.vg -o ${OUT}.no_ref.gbwt
    vg gbwt -o ${OUT}.with_ref.gbwt --set-reference "$SAMPLE_ID" ${OUT}.no_ref.gbwt
    vg gbwt --gbz-format -x ${OUT}.vg ${OUT}.with_ref.gbwt -g ${OUT}.gbz
    # Now that we have the GBZ we don't need these other things
    rm ${OUT}.no_ref.gbwt ${OUT}.with_ref.gbwt ${OUT}.vg
    # Index
    vg autoindex --gbz ${OUT}.gbz -w lr-giraffe --prefix "$OUT" --no-guessing
    # Minimap2 needs FASTA input
    vg paths --extract-fasta -x ${OUT}.gbz > ${OUT}.fasta
    # Avoid reusing an old index
    rm -f "${OUT}.fasta.fai"
    minimap2 -x map-hifi -d ${OUT}.mmi ${OUT}.fasta
fi