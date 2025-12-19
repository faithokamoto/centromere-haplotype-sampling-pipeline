#!/bin/bash
# Run the 21 test samples through the pipeline
# Usage: run_test_samples.sh

while IFS= read -r line; do
    haplotype=$(echo "$line" | cut -f3 -d ",")
    ./guess_optimal_num_sampled_haplo.py "$haplotype"
    echo "============="
done < test_samples.txt