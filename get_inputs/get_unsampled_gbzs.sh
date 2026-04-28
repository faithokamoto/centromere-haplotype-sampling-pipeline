#!/bin/bash
# For all Centrolign graphs, get a single .gbz
# Usage: get_unsampled_gbzs.sh

GFA_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA
GBZ_DIR=/private/groups/patenlab/fokamoto/centrolign/graph/unsampled

CHM13_NAME="CHM13#0#CHM13.0#0"
CHM13_DIR=/private/groups/patenlab/fokamoto/centrolign/graph/haploid

for chrom in {1..22} X; do
    chrom="chr${chrom}"
    
    # Create combined GFA
    if [ -f $GFA_DIR/$chrom/${chrom}.centrolign.gfa ]; then
        ./get_inputs/add_dummy_caps.py -o $GBZ_DIR/${chrom}.gfa $GFA_DIR/$chrom/${chrom}.centrolign.gfa
    else
        ./get_inputs/add_dummy_caps.py -o $GBZ_DIR/${chrom}.gfa $GFA_DIR/$chrom/*/*.gfa
    fi

    # Convert to GBZ
    ./get_inputs/gfa_to_gbz_ref.sh $GBZ_DIR/$chrom
    # Extract CHM13 ref
    ./get_inputs/create_single_path_ref.sh $GBZ_DIR/${chrom}.pg "$CHM13_NAME" $CHM13_DIR/${chrom}.CHM13
done