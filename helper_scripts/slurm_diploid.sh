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
#SBATCH --mem=4gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=20
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=5:00
#
# Array job specification:
#SBATCH --array=1-785

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

GFA=/private/groups/patenlab/fokamoto/centrolign/graph/unsampled/chr12.gfa
combo=`ls $HAP_GRAPH_DIR/*.*.*.gbz | cut -f9 -d "/" | cut -f1-2 -d "." | \
    uniq -c | fgrep -v "1 chr" | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1`
chrom=`echo "$combo" | cut -f1 -d "."`
sample_id=`echo "$combo" | cut -f2 -d "."`
echo "Running $sample_id on $chrom"

log=$DIR/log/${chrom}.${sample_id}.typing.log
$DIR/diploid_paper_typing.sh "$sample_id" "$chrom" &> "$log"
