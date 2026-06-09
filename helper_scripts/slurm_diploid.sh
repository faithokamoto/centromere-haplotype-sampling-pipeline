#!/bin/bash
# Job name:
#SBATCH --job-name=diploid_alignment
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
#SBATCH --mem=25gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=20
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=8:00:00
#
# Array job specification:
#SBATCH --array=1-2686%100

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

READS_DIR=/private/groups/patenlab/fokamoto/centrolign/to_align
combo=`ls $READS_DIR/*.real.fastq.gz | cut -f8 -d "/" | cut -f1-2 -d "." | uniq -c | fgrep -v "1 chr" \
    | head -n "$SLURM_ARRAY_TASK_ID" | tail -1 | sed "s/      2 //g"`
chrom=`echo "$combo" | cut -f1 -d "."`
sample_id=`echo "$combo" | cut -f2 -d "."`
echo "Running $sample_id on $chrom"

log=$DIR/log/diploid_${chrom}/${chrom}.${sample_id}.log
$DIR/diploid_paper_alignments.sh "$sample_id" "$chrom" &> "$log"