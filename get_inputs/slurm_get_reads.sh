#!/bin/bash
# Job name:
#SBATCH --job-name=get_reads
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
#SBATCH --ntasks=1
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=8:00:00
#
# Array job specification:
#SBATCH --array=2-232%50

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

READ_LOCS=$DIR/input_data/aws_file_locations.csv
sample_id=`head -n "$SLURM_ARRAY_TASK_ID" "$READ_LOCS" | tail -n 1 | cut -f1 -d ","`
echo "Running sample: $sample_id"

log=$DIR/log/get_reads/${sample_id}.get_reads.log
$DIR/get_inputs/get_reads.sh $sample_id &> $log
