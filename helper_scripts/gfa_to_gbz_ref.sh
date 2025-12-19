#!/bin/bash
# Prepare a GFA for being used a haplotype sampling reference
# Usage: gfa_to_gbz_ref.sh <graph prefix>
# Example: gfa_to_gbz_ref.sh $PROJ_DIR/graph/unsampled/chr12

GRAPH=$1

# PG format needs nodes with length <= 1024bp for distance indexing
vg convert --gfa-in $GRAPH.gfa | vg mod --chop 1024 - > $GRAPH.pg
# Convert to GBZ format by smooshing GBWT and PG together
vg gbwt --index-paths -x $GRAPH.pg -o $GRAPH.gbwt
vg gbwt --gbz-format -x $GRAPH.pg $GRAPH.gbwt -g $GRAPH.gbz

vg gbwt -r $GRAPH.ri -Z $GRAPH.gbz
# Distanceless index is much quicker to make and works for vg haplotypes
vg index -t 64 -w 1 -w 2 --no-nested-distance -j $GRAPH.dist $GRAPH.gbz
vg haplotypes -v 3 -t 16 -H $GRAPH.hapl -d $GRAPH.dist -r $GRAPH.ri $GRAPH.gbz