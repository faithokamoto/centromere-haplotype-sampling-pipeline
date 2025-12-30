#!/bin/bash
# Job name:
#SBATCH --job-name=chr12_centrolign
#
# Partition - This is the queue it goes in:
#SBATCH --partition=medium
#
# Where to send email (optional)
#SBATCH --mail-user=fokamoto@ucsc.edu
#
# Number of nodes you need per job:
#SBATCH --nodes=1
#
# Memory needed for the jobs.  Try very hard to make this accurate.  DEFAULT = 4gb
#SBATCH --mem=60gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=20
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=6:00:00
#
# Array job specification:
#SBATCH --array=1-21

READS_DIR=/private/groups/patenlab/fokamoto/centrolign/to_align
CODE_DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline
SCRIPT=$CODE_DIR/code/sim_and_real_LOO.sh

line=`head -n $SLURM_ARRAY_TASK_ID $CODE_DIR/test_samples.txt | tail -n 1`
sample_id=$(echo "$line" | cut -f1 -d "," | cut -f1 -d ".")
echo "Running sample: $sample_id"
$CODE_DIR/leave_one_out_alignments_diploid.sh "$sample_id"
$CODE_DIR/guess_optimal_num_sampled_haplo.py \
 --hap1-reads $READS_DIR/real_${sample_id}_chr12_hor_array.hifi.1.fastq \
 --hap2-reads $READS_DIR/real_${sample_id}_chr12_hor_array.hifi.2.fastq \
 $sample_id