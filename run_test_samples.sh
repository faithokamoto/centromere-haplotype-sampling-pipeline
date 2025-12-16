#!/bin/bash
# Run the 21 test samples through the pipeline
# Usage: run_test_samples.sh

while IFS= read -r line; do
    ./guess_optimal_num_sampled_haplo.py "$line"
    echo "============="
done < test_samples.txt