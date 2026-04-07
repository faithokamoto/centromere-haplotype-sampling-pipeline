#!/bin/bash
# Job name:
#SBATCH --job-name=centrolign
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
#SBATCH --mem=100gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=20
#
# Standard output and error log
#SBATCH --output=/private/home/fokamoto/centromere-haplotype-sampling-pipeline/log/slurm%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=6:00:00
#
# Array job specification:
#SBATCH --array=1-1984

# Activate Conda environment
source /private/home/${USER}/.bashrc
source activate cenhap-sample

HAP_GRAPH_DIR=/private/groups/patenlab/fokamoto/centrolign/graph/haploid
cur_hap_graph=`ls $HAP_GRAPH_DIR/*.*.1.gbz $HAP_GRAPH_DIR/*.*.2.gbz | head -n "$SLURM_ARRAY_TASK_ID" | tail -n 1 | cut -f9 -d "/"`
chrom=`echo "$cur_hap_graph" | cut -f1 -d "."`
hap_name=`echo "$cur_hap_graph" | cut -f2,3 -d "."`
echo "Running $hap_name on $chrom"

DIR=/private/home/fokamoto/centromere-haplotype-sampling-pipeline

if [ `echo $hap_name | cut -f1 -d "."` == "HG002" ]
then
    echo "HG002; ignore"
else
    log=$DIR/log/${chrom}/${chrom}.${hap_name}.log
    $DIR/haploid_paper_alignments.sh "$hap_name" "$chrom" &> "$log"
fi


