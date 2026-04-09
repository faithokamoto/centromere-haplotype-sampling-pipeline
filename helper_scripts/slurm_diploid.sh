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
#SBATCH --mem=10gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=10:00
#
# Array job specification:
#SBATCH --array=1-784

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

READS_DIR=/private/groups/patenlab/fokamoto/centrolign/to_align
combo=`ls $READS_DIR/*.real.fastq.gz | cut -f8 -d "/" | cut -f1-2 -d "." | uniq -c \
    | fgrep -v "1 chr" | head -n "$SLURM_ARRAY_TASK_ID" | tail -1 | sed "s/      2 //g"`
chrom=`echo "$combo" | cut -f1 -d "."`
sample_id=`echo "$combo" | cut -f2 -d "."`
echo "Running $sample_id on $chrom"

log=$DIR/log/diploid_typing/${chrom}.${sample_id}.typing.log
$DIR/diploid_paper_typing.sh "$sample_id" "$chrom" &> "$log"