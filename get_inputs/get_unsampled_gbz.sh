#!/bin/bash
# For all Centrolign graphs, get a single .gbz
# Usage: get_unsampled_gbz.sh <chrom>
# Example: get_unsampled_gbz.sh chr12

CHROM=$1

GFA_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA
GBZ_DIR=/private/groups/patenlab/fokamoto/centrolign/graph/unsampled

CHM13_NAME="CHM13#0#CHM13.0#0"
CHM13_DIR=/private/groups/patenlab/fokamoto/centrolign/graph/haploid

# Create combined GFA
if [ -f $GFA_DIR/$CHROM/${CHROM}.centrolign.gfa ]; then
    ./get_inputs/add_dummy_caps.py -o $GBZ_DIR/${CHROM}.gfa $GFA_DIR/$CHROM/${CHROM}.centrolign.gfa
else
    ./get_inputs/add_dummy_caps.py -o $GBZ_DIR/${CHROM}.gfa $GFA_DIR/$CHROM/*/*.gfa
fi

# Convert to GBZ
./get_inputs/gfa_to_gbz_ref.sh $GBZ_DIR/$CHROM
# Extract CHM13 ref
./get_inputs/create_single_path_ref.sh $GBZ_DIR/${CHROM}.pg "$CHM13_NAME" $CHM13_DIR/${CHROM}.CHM13