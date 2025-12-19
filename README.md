# Centromere haplotype sampling pipeline

> [!WARNING]
> Under active development.

This is a pipeline which aligns centromeric reads against a pangenome graph
which has undergone [haplotype sampling][HaplotypeSampling].

## Dependencies

- [`vg`][vg]: currently commit `c4c189272dc02750774bac921c2c5996e486ad7d`
- [Python 3][Python]: currently using v3.14.1

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
seed hits and thus makes read alignment possible.

I currently have a manual method to select *n* from 1-8 by looking at some
output plots. However, this is time-consuming and relies on human judgement. In
this repository I will attempt to develop an heuristic algorithm which selects
an optimal *n* by leveraging the intuition I have developed from my ad hoc work.

## Workflow

Starting from an empty directory and installations of the dependencies:
1. Acquire a Centrolign GFA (or GFAs). The chr12 ones I'm using are in:
    * `/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_0/chr12.subgroup_0.centrolign.gfa`
    * `/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/MSA/chr12/subgroup_1/chr12.subgroup_1.centrolign.gfa`
2. Connect dummy nodes to the source and sink of each haplotype path. This step
is necessary because all haplotypes must participate in the same top-level chain
for the haplotype sampling algorithm. This uses `add_dummy_caps.py`
3. Convert to a GBZ and index for haplotype sampling: `gfa_to_gbz_ref.sh`
4. Perform alignments (currently done with `leave_one_out_alignments.sh`)
5. Run `guess_optimal_num_sampled_haplo.py` for each haplotype to get *n*.

`./run_test_samples.sh` should guess, for each input haplotype in the first
column, that the optimal *n* is as in the second column.

## Files

- **Workflow**
    - `add_dummy_caps.py`: modify Centrolign GFAs to input to haplotype sampling
    - `align_reads_giraffe.sh`: align reads against a graph via `vg giraffe`
    - `align_reads_minimap2.sh`: align reads against a linear ref via `minimap2`
    - `gfa_to_gbz_ref.sh`: prepare Centrolign GFA as a GBZ reference graph
    - `guess_optimal_num_sampled_haplo.py`: guess optimal *n* value
    - `leave_one_out_alignments.sh`: run the leave-one-out alignments
- **Testing**
    - `run_test_samples.sh`: script to run the guesser on all test haplotypes
    - `test_samples.txt`: CSV with 21 hand-selected test samples; columns are
    path name, BED file version number, haplotype name, and optimal *n*
- **Metadata**
    - `.gitignore`: some files that I don't feel like putting on version control
    - `LICENSE`: the MIT license as it applies to this repository
    - `README.md`: this file, which explains the repository

## TODO

- [ ] Make basic and probably very bad guessing script
    - [X] Decide which samples/haplotypes are hopeless
    - [X] Rank based on identity scores
    - [ ] Rank (tiebreak?) based on node usage stats
- [X] Can the haplotype sampler itself guess when to stop?
    **NO**: no obvious drop
    - [X] Plot score changes around optimal number
- [ ] Find "correctness" using distance matrix / nearby samples
    - [ ] Plot # of neighbors found within optimal subsampling
        - Also care about "fake" neighbors that were swept up?
        - % of total that could be found?
        - 0.15 threshold?
    - [X] Plot dist to nearest neighbor vs. dist to nearest subsampled
- [ ] Have pipeline run diploid samples instead of haploid
- [ ] Attempt to do typing

[Centrolign]: https://github.com/jeizenga/centrolign
[HaplotypeSampling]: https://github.com/vgteam/vg/wiki/Haplotype-Sampling
[LRgiraffe]: https://doi.org/10.1101/2025.09.29.678807
[matplotlib]: https://matplotlib.org/
[Python]: https://www.python.org/downloads
[vg]: https://github.com/vgteam/vg