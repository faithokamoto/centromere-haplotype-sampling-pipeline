#!/bin/bash
# Prepare a GFA for being used a haplotype sampling reference
# Usage: gfa_to_gbz_ref.sh <graph prefix>
# Example: gfa_to_gbz_ref.sh $PROJ_DIR/graph/unsampled/chr12

GRAPH=$1

# PG format needs nodes with length <= 1024bp for distance indexing
vg convert --gfa-in $GRAPH.gfa | vg mod --chop 1024 - > $GRAPH.pg
# Augmented reference coordinates to use for variant calling
vg paths -Q CHM13 -x $GRAPH.pg \
    --compute-augref --augref-sample augref_CHM13 --min-augref-len 50 \
    --augref-segs $GRAPH.augref.segs.tsv > $GRAPH.augref.vg
# Convert to GBZ format by smooshing GBWT and VG together
vg gbwt --index-paths -x $GRAPH.augref.vg -o $GRAPH.gbwt
vg gbwt --gbz-format -x $GRAPH.augref.vg $GRAPH.gbwt -g $GRAPH.gbz

vg gbwt -r $GRAPH.ri -Z $GRAPH.gbz
# Distanceless index is much quicker to make and works for vg haplotypes
vg index -t 64 -w 1 -w 2 --no-nested-distance -j $GRAPH.dist $GRAPH.gbz
vg haplotypes -v 3 -t 16 -H $GRAPH.hapl -d $GRAPH.dist -r $GRAPH.ri $GRAPH.gbz