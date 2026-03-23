#!/bin/bash
# Job name:
#SBATCH --job-name=diploid_type
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
#SBATCH --time=10:00
#
# Array job specification:
#SBATCH --array=1-159

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

GFA=/private/groups/patenlab/fokamoto/centrolign/graph/unsampled/chr12.gfa
sample_id=`grep "^P" "$GFA" | cut -f1 -d "#" | sort | uniq -c | fgrep -v "1 P" \
    | cut -f2 | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1`
echo "Running sample: $sample_id"

log=$DIR/log/${sample_id}.typing.log
$DIR/haploid_paper_alignments.sh $sample_id chr12 &> $log
