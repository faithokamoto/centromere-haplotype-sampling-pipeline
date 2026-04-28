#!/bin/bash
# Job name:
#SBATCH --job-name=get_unsampled_gbz
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
#SBATCH --mem=25gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=30:00
#
# Array job specification:
#SBATCH --array=1-23

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

if [ "$SLURM_ARRAY_TASK_ID" == "23" ]; then
    chrom="chrX"
else
    chrom="chr${SLURM_ARRAY_TASK_ID}"
fi

log=$DIR/log/get_unsampled_gbz/${chrom}.log
$DIR/get_inputs/get_unsampled_gbz.sh $chrom &> $log
