# Centromere haplotype sampling pipeline

This pipeline aligns centromeric reads against a pangenome graph which has
undergone [haplotype sampling][HaplotypeSampling].

## Dependencies

- [`vg`][vg] currently version `v1.73.0-33-gbf85ab04b` with extra haplotype
    sampling information ([this print][ExtraPrint]) turned on

Most commands are run within this Conda environment:
```bash
conda create -n cenhap-sample -c conda-forge -c bioconda \
    kmc matplotlib minimap2 samtools
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
1. Get graph & CHM13 references: `get_unsampled_gbzs.sh`
2. Acquire reads & linear references: `get_reads.sh` (via `slurm_get_reads.sh`)
3. Refine parameters (`param_test`):
    - `default_param_alignments.sh` (via `slurm_alignments.sh`)
    - `testing_absent_score.sh`
    - `absent_score_fig.py`
4. Perform alignment and typing experiments:
    - `haploid_paper_alignments.sh` (via `slurm_haploid.sh`)
    - `diploid_paper_typing.sh` (via `slurm_diploid.sh`)
5. Collect experiment data (`data_scripts`)
    - `haploid_data.py`
    - `diploid_data.py`
6. Analayze HG01106.1 on chr10 (`hap_focus`)
    - `retain_files_alignments.sh`
    - `massage_inputs.sh`
    - `alignment_details_fig.py`
    - Finish up with two outside operations:
        - Take the annotated GFA and visualize with [BandageNG][Bandage].
        - Add a curly brace to the `alignment_details_fig.py` output.

## Inputs

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
    - `get_unsampled_gbzs.sh`: prepare main references (unsampled graph & CHM13)
    - `gfa_to_gbz_ref.sh`: prepare Centrolign GFA as a GBZ reference graph
    - `slurm_get_reads.sh`: run `get_reads.sh` on haplo/chrom pairs
- **Focused analysis of one haplotype** (`hap_focus`)
    - `alignment_details_fig.py`: create chr10 HG01106.1 figure
    - `annotate_depths.py`: add `DP` tags to a GFA, given `vg pack` output
    - `massage_inputs.sh`: reshape default alignment outputs for easier plotting
    - `retain_files_alignments.sh`: simplified alignments with less cleanup
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

Various scripts assume that certain files exist in `input_data`, and that
certain directories exist in `log` (`chrX` for each chromosome, `default`,
`diploid_typing`, `get_reads`, `slurm`, and `typing_tests` for me).

[Bandage]: https://github.com/asl/BandageNG
[Centrolign]: https://github.com/jeizenga/centrolign
[ExtraPrint]: https://github.com/vgteam/vg/blob/bf85ab04b251e7a2bc308750d6c8c44afda213f5/src/recombinator.cpp#L12
[HaplotypeSampling]: https://github.com/vgteam/vg/wiki/Haplotype-Sampling
[ModDotPlot]: https://github.com/marbl/ModDotPlot
[LRgiraffe]: https://doi.org/10.1101/2025.09.29.678807
[SelectionCode]: https://github.com/vgteam/vg/blob/2e664f07e49caca29a208b3b0f2f25c7100df5e9/src/recombinator.cpp#L2056-L2058
[vg]: https://github.com/vgteam/vg