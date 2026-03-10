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
#SBATCH --time=8:00:00
#
# Array job specification:
#SBATCH --array=1-100%10

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

GFA=/private/groups/patenlab/fokamoto/centrolign/graph/unsampled/chr12.gfa
path_name=`grep "^P" "$GFA" | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1 | cut -f2 | cut -d "#" -f3`
echo "Running path: $path_name"

log=$DIR/log/${path_name}.log
$DIR/haploid_paper_alignments.sh $path_name &> $log
