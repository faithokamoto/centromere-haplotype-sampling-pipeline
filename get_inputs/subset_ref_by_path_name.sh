#!/bin/bash
# Extract one or two paths from larger graph and index it for alignment
# Usage: subset_ref_by_path_name.sh <big graph> <path name(s)> <small graph prefix>
# Example: subset_ref_by_path_name.sh $PROJ_DIR/graph/unsampled/chr12.pg CHM13#0#CHM13.0#0 $PROJ_DIR/graph/haploid/chr12.CHM13

set -e

BIG_GRAPH=$1
if [ $# == 3 ]; then
    path_arg="--keep-path $2"
    OUT=$3
    num_paths=haploid
    sample_id=`echo "$2" | cut -f1 -d "#"`
else
    path_arg="--keep-path $2 --keep-path $3"
    OUT=$4
    num_paths=diploid
fi

# Only bother if needed
vg mod $path_arg $BIG_GRAPH > ${OUT}.vg
# Convert to GBZ
vg gbwt --index-paths -x ${OUT}.vg -o ${OUT}.no_ref.gbwt

if [ "$num_paths" == "haploid" ]; then
    # Stick graph will act as a reference
    vg gbwt -o ${OUT}.with_ref.gbwt --set-reference "$sample_id" ${OUT}.no_ref.gbwt
    vg gbwt --gbz-format -x ${OUT}.vg ${OUT}.with_ref.gbwt -g ${OUT}.gbz
else
    # Diploid graph will have no reference
    vg gbwt --gbz-format -x ${OUT}.vg ${OUT}.no_ref.gbwt -g ${OUT}.gbz
fi

# Now that we have the GBZ we don't need these other things
rm ${OUT}.*.gbwt ${OUT}.vg
# Index
vg autoindex --gbz ${OUT}.gbz -w lr-giraffe --prefix "$OUT" --no-guessing

if [ "$num_paths" == "haploid" ]; then
    # Minimap2 needs FASTA input
    vg paths --extract-fasta -x ${OUT}.gbz > ${OUT}.fasta
    # Avoid reusing an old index
    rm -f "${OUT}.fasta.fai"
    minimap2 -x map-hifi -d ${OUT}.mmi ${OUT}.fasta
fi