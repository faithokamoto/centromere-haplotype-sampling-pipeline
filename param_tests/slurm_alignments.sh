#!/bin/bash
# Job name:
#SBATCH --job-name=centrolign
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
#SBATCH --mem=125gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=20
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=1:00:00
#
# Array job specification:
#SBATCH --array=1-122

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

READS_DIR=/private/groups/patenlab/fokamoto/centrolign/to_align
cur_file=`ls $READS_DIR/chr4.*.real.fastq.gz | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1 | cut -f8 -d "/"`
hap_name=`echo "$cur_file" | cut -f2,3 -d "."`
echo "Running $hap_name default param alignments"

LOCAL_DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline
log=$LOCAL_DIR/log/default/${hap_name}.log
$LOCAL_DIR/param_tests/default_param_alignments.sh "$hap_name" &> "$log"