# Centromere haplotype sampling pipeline

This pipeline aligns centromeric reads against a pangenome graph which has
undergone [haplotype sampling][HaplotypeSampling].

## Dependencies

- [`vg`][vg] currently version `v1.73.0-33-gbf85ab04b` with extra haplotype
    sampling information ([this print][ExtraPrint]) turned on
- [`ModDotPlot`][ModDotPlot] version `v0.9.9` installed from source and using
    its default environemnt

Most commands are run within this Conda environment:
```bash
conda create -n cenhap-sample -c conda-forge -c bioconda \
    kmc matplotlib minimap2 samtools
```

The ModDotPlot commands, at the end of `moddotplot/massage_inputs.sh`, are run
within its default virtual environment, and NOT within the Conda environment.
```bash
git clone https://github.com/marbl/ModDotPlot.git
cd ModDotPlot
python -m venv venv
source venv/bin/activate
```

## Background

Centromeres are composed of highly complex arrays of tandem repeats. Unlike most
of the genome, while some sequence similarity exists within centromere haplotype
groups ("cenhaps"), high divergence in between groups makes cross-cenhap
alignments fairly biologically meaningless. Thus any single linear reference, no
matter how high quality, will lead to poor read-to-reference mappings.

The widely used read-to-graph aligner `vg giraffe` was [updated][LRgiraffe] to
work with long reads. [Centrolign][Centrolign] creates high-quality pangenome
graphs of centromeres. This pipeline combines these developments to align
centromeric long reads to a pangenome reference graph.

However, this graph is very complex. The chr12 reference has 373 haplotypes.
Aligning directly to the reference graph, or even indexing it, is impractical.

Haplotype sampling entails selecting a personalized set of *n* haplotypes which
best match the input read set's *k*-mers. The full reference graph is then
subset to only those haplotype paths. This dramatically decreases the number of
hits and thus makes read alignment possible. Selection of *n* is non-trivial.

## Workflow

What I did:
1. Acquire Centrolign GFAs, e.g. 
    `/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/chr4.centrolign.gfa`
2. Connect dummy nodes to the source and sink of each haplotype path. This step
is necessary because all haplotypes must participate in the same top-level chain
for the haplotype sampling algorithm. This uses `add_dummy_caps.py`
3. Convert to a GBZ and index for haplotype sampling: `gfa_to_gbz_ref.sh`
4. Extract CHM13 reference from the larger graph: `create_single_path_ref.sh`
5. Acquire reads & linear references: `get_reads.sh` (via `slurm_get_reads.sh`)
6. Refine parameters (`param_test`):
    - `default_param_alignments.sh` (via `slurm_alignments.sh`)
    - `testing_absent_score.sh`
    - `absent_score_fig.py`
7. Perform alignment and typing experiments:
    - `haploid_paper_alignments.sh` (via `slurm_haploid.sh`)
    - `diploid_paper_typing.sh` (via `slurm_diploid.sh`)
8. Collect experiment data (`data_scripts`)
    - `haploid_data.py`
    - `diploid_data.py`

## Inputs

Note that I use "all chromosomes" here, but that just means the ones I used,
which are chr4, chr6, chr9, chr10, chr11, chr12, and chr17.

- **Centrolign GFAs**
- **Distance matrices** (headerless CSV hap1,hap2,dist)  
  `input_dir/chr*_r2_QC_v2_centrolign_pairwise_distance.csv`
- **Truth cenhap assignments** (header TSV haplotype & cenhap)  
  `input_dir/chr*.cenhap_predictions.tsv`
- **Table of HPRC BAM locations** (CSV, first col is sample name, third is BAM)  
  `input_dir/aws_file_locations.csv`
- **BED array annotations for all haplotypes** (external, in Mira's dir)
  `/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds/`
- **Graph-unique k-mer counts per chr4 haplotype** (for parameter search)  
  `input_data/chr4.txt`
- **CIGAR strings for some chr10 HG01106.1 assembly alignments**  
    - `input_data/pairwise_cigar_CHM13.0_HG01106.1.txt`
    - `input_data/pairwise_cigar_HG01106.1_HG01891.2.txt`
- **Arial** the font, `input_data/arial.ttf`

## Files

- **Data collection scripts** (`data_scripts`)
    - `haploid_data.py`: get haploid alignment & typing stats
    - `diploid_data.py`: get diploid typing stats
- **Preprocessing scripts** (`get_inputs`)
    - `add_dummy_caps.py`: modify Centrolign GFAs to input to haplotype sampling
    - `create_single_path_ref.sh`: extract & index a linear ref from graph
    - `edit_sam.py`: update a SAM file from genome-wide to graph-space
    - `get_reads.sh`: prepare reads & linear references for a haplotype
    - `gfa_to_gbz_ref.sh`: prepare Centrolign GFA as a GBZ reference graph
    - `slurm_get_reads.sh`: run `get_reads.sh` on haplo/chrom pairs
- **Main pipeline helpers** (`helper_scripts`)
    - `align_reads_giraffe.sh`: align reads with Giraffe & get per-read stats
    - `align_reads_minimap2.sh`: align reads with Minimap2 & get per-read stats
    - `calculate_alignment_stats.py`: get stats for one haplo/chrom pair
    - `slurm_diploid.sh`: run `diploid_paper_typing.sh` on sample/chrom pairs
    - `slurm_haploid.sh`: run `haploid_paper_alignments.sh` on haplo/chrom pairs
- **`--absent-score` parameter testing** (`param_tests`)
    - `absent_score_fig.py`: create `--absent-score` supplementary figure
    - `arial.ttf`: font to make the figure fit within *Nature* rules
    - `calculate_private_depths.py`: get depth on private nodes of top haplotype
    - `default_param_alignments.sh`: align to default-parameter sampled graph
    - `slurm_alignments.sh`: run `default_param_alignments.sh` on chr4 haplos
    - `testing_absent_score.sh`: run chr4 typing with varying `--absent-score`
- **Main pipeline/methods**:
    - `diploid_paper_typing.sh`: run typing experiment for one sample/chrom pair
    - `guess_n_and_cenhap.py`: guess optimal *n* value and cenhap of input
    - `haploid_paper_alignments.sh`: run alignments for one haplo/chrom pair
- **Metadata**
    - `.gitignore`: some files that I don't feel like putting on version control
    - `LICENSE`: the MIT license as it applies to this repository
    - `README.md`: this file, which explains the repository

[Centrolign]: https://github.com/jeizenga/centrolign
[ExtraPrint]: https://github.com/vgteam/vg/blob/bf85ab04b251e7a2bc308750d6c8c44afda213f5/src/recombinator.cpp#L12
[HaplotypeSampling]: https://github.com/vgteam/vg/wiki/Haplotype-Sampling
[ModDotPlot]: https://github.com/marbl/ModDotPlot
[LRgiraffe]: https://doi.org/10.1101/2025.09.29.678807
[SelectionCode]: https://github.com/vgteam/vg/blob/2e664f07e49caca29a208b3b0f2f25c7100df5e9/src/recombinator.cpp#L2056-L2058
[vg]: https://github.com/vgteam/vg