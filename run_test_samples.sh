#!/bin/bash
# Run the 21 test samples through the pipeline
# Usage: run_test_samples.sh

while IFS= read -r line; do
    path_name=$(echo "$line" | cut -f1 -d ",")
    version=$(echo "$line" | cut -f2 -d ",")
    haplotype=$(echo "$line" | cut -f3 -d ",")
    echo "Running path: $path_name, version: $version"
    ./leave_one_out_alignments.sh chr12 "$path_name" "$version"
    ./guess_optimal_num_sampled_haplo.py "$haplotype"
    echo "============="
done < test_samples.txt