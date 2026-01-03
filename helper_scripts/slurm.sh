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
#SBATCH --time=4:00:00
#
# Array job specification:
#SBATCH --array=1-21

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

line=`head -n $SLURM_ARRAY_TASK_ID $DIR/input_data/test_samples.txt | tail -n 1`
path_name=$(echo "$line" | cut -f1 -d ",")
version=$(echo "$line" | cut -f2 -d ",")
haplotype=$(echo "$line" | cut -f3 -d ",")
echo "Running path: $path_name, version: $version"
log=$DIR/log/${path_name}.log
$DIR/leave_one_out_alignments.sh "$path_name" "$version" &> $log
$DIR/plot_scripts/plot_identity_and_accuracy.py \
    --name $haplotype --logfile $log --output_file $DIR/plot_outputs/${path_name}.png
