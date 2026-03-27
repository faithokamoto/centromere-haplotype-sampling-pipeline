# Centromere haplotype sampling pipeline

This pipeline aligns centromeric reads against a pangenome graph which has
undergone [haplotype sampling][HaplotypeSampling].

## Dependencies

- [`vg`][vg]: currently version v`1.72.0` with a few modifications

Then do
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
Thus, aligning directly to the reference graph is impractical. Many possibly
informative seeds are dropped due to the "high hit cap", i.e. appearing too many
times in the read and reference. We can't shorten the reads without surrending
the very longness that makes them valuable. However, we can simplify the graph.

Haplotype sampling entails selecting a personalized set of *n* haplotypes which
best match the input read set's *k*-mers. The full reference graph is then
subset to only those haplotype paths. This dramatically decreases the number of
hits and thus makes read alignment possible. Selection of *n* is non-trivial.

## Workflow

What I did:
1. Acquire a Centrolign GFA (or GFAs), e.g. 
    `/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr4/chr4.centrolign.gfa`
2. Connect dummy nodes to the source and sink of each haplotype path. This step
is necessary because all haplotypes must participate in the same top-level chain
for the haplotype sampling algorithm. This uses `add_dummy_caps.py`
3. Convert to a GBZ and index for haplotype sampling: `gfa_to_gbz_ref.sh`
4. Extract CHM13 reference from the larger graph: `create_single_path_ref.sh`
5. Acquire reads & linear references: `get_reads.sh` via `slurm_get_reads.sh`
6. Perform alignment and typing experiments:
    - `haploid_paper_alignments.sh` (called with `slurm_haploid.sh`)
    - `diploid_paper_typing.sh` (called with `slurm_diploid.sh`)

## Files

- **Workflow**
    - `add_dummy_caps.py`: modify Centrolign GFAs to input to haplotype sampling
    - `create_single_path_ref.sh`: extract & index a linear ref from graph
    - `gfa_to_gbz_ref.sh`: prepare Centrolign GFA as a GBZ reference graph
    - `guess_n_and_cenhap.py`: guess optimal *n* value and cenhap of input
    - `haploid_paper_alignments.sh`: run the alignments
- **Plotting**
    - `haploid_data.py`: calculate stats using all the outputs
    - `plot_identity_and_accuracy.py`: plot ID % and SNV calling across *n*s
    - `plot_variants_vs_dist.py`: four-panel plot for a sample with variant info
- **Metadata**
    - `.gitignore`: some files that I don't feel like putting on version control
    - `LICENSE`: the MIT license as it applies to this repository
    - `README.md`: this file, which explains the repository

[Centrolign]: https://github.com/jeizenga/centrolign
[HaplotypeSampling]: https://github.com/vgteam/vg/wiki/Haplotype-Sampling
[LRgiraffe]: https://doi.org/10.1101/2025.09.29.678807
[SelectionCode]: https://github.com/vgteam/vg/blob/2e664f07e49caca29a208b3b0f2f25c7100df5e9/src/recombinator.cpp#L2056-L2058
[vg]: https://github.com/vgteam/vg