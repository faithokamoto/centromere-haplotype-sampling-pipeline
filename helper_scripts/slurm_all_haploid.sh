#!/bin/bash
# Job name:
#SBATCH --job-name=chr12_centrolign
#
# Partition - This is the queue it goes in:
#SBATCH --partition=short
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
#SBATCH --time=1:00:00
#
# Array job specification:
#SBATCH --array=1-373%25

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline
GFA=/private/groups/patenlab/fokamoto/centrolign/graph/unsampled/chr12.gfa

path_name=`fgrep P $GFA | cut -f2 | head -$SLURM_ARRAY_TASK_ID | tail -1 | cut -d"#" -f3`
echo "Running path: $path_name"

log=$DIR/log/${path_name}_guess.log
$DIR/sample_for_cenhap_guess.sh "$path_name" &> $log
$DIR/guess_cenhap.py --ploidy 1 --logfile $log --max-n 10 --force-n 3 &>> $log