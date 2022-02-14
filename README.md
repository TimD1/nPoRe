# nPoRe: n-Polymer Realigner for improved pileup variant calling

## Introduction
`npore` is a read realigner which recalculates each read's fine-grained alignment in order to more accurately align ``n-polymers'' such as homopolymers (n=1) and tandem repeats (2 &leq; n &leq; 6). In other words, given an input BAM, it adjusts each read's CIGAR string to more accurately model the most likely sequencing errors and actual variants. Traditional affine gap penalties are context-agnostic, and do not model the higher likelihood of INDELs in low-complexity regions (particularly n-polymers), leading to poor or inconsistent alignments. We find that `npore` improves pileup concordance across reads and results in slightly better variant calling performance.
<div align="center">
<img src="img/npore_pileup.png" width="720p" alt="npore vs minimap2 pileup comparison">
</div>


## Citation
Please cite the following pre-print if you use `npore`:

<details>
<summary>
<a href=""><b>[bioRxiv]</b> nPoRe: n-Polymer Realigner for improved pileup variant calling</a>
</summary>

Authors: Tim Dunn, David Blaauw, Reetuparna Das, Satish Narayansamy
</details>

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Project Structure](#project-structure)
* [Data Sources](#data-sources)
* [Acknowledgements](#acknowledgements)

## Installation
Installation requires a working version of `virtualenv` and `python3`. Python packages in `requirements.txt` will be downloaded as part of the installation. Simply run:

```bash
git clone https://github.com/TimD1/npore && cd npore
make npore
```

This makefile will create a virtual environment, download the required packages, and build `npore`. After setup, run:

```bash
. ./venv3/bin/activate
python3 ./src/realign.py --help
```

to verify the setup has succeeded.

## Usage
Here's an example usage of the main `realign.py` program, which will store results in `realigned.sam`.

```bash
python3 ~/npore/src/realign.py \
    --bam reads.bam \
    --ref ref.fasta \
    --out_prefix realigned \
    --contigs chr20,chr21,chr22 \
    --stats_dir ~/npore/guppy5_stats
```
For additional options, run `python3 realign.py --help`.

## Project Structure
`src/` | nPoRe source code
---: | ---
`realign.py` | Module for realigning a BAM file.
`standardize_vcf.py` | Module for standardizing a ground-truth VCF file to report variants in the same manner that nPoRe would align the reads.
`bed.py` | Module for computing n-polymer BED regions.
`purity.py` | Module for computing a BAM pileup's Gini purity, for measuring read concordance.
`filter.py` | Simple module for filtering overlapping variants.
`aln.pyx` | Contains alignment-related functions.
`bam.pyx` | Contains BAM-related functions.
`cig.pyx` | Contains CIGAR-related functions.
`cfg.py` | Contains global variables and configuration info.
`vcf.py` | Contains VCF-related functions.
`util.py` | Contains helper functions.

`scripts/` | Helper scripts used during evaluation, example usage
---: | ---
`realign_pipeline.sh` | Main Clair3 retraining pipeline.
`happy.sh` | Runs `hap.py` evaluation of all configurations/regions.
`plot_results.py` | Plots final precision/recall graphs.
`plot_sankey.py` | Generates Sankey plot of actual/error INDELs by n-polymer BED region.
`calc_beds.sh` | Calculates n-polymer BEDS, running `bed.py`.
`sankey.py` | Custom Sankey plot library, extended from <a href="https://github.com/anazalea/pySankey">`pySankey`</a>.
`purity.sh` | Calculates Gini purity.
`align.sh` | Aligns reads to a reference, allowing multiple input formats.
`tag_unphased.py` | Tags unphased reads with `HP:i:0`.

`test/` | Testing directory
---: | ---
`align.py` | Tests `align()` kernel.
`get_np_info.py` | Tests n-polymer info generation.
`realign.sh` | Tests full read realignment.
`test_std_vcf.sh` | Tests VCF standardization.
`profile_alignment.ipynb` | Line-by-line profiling of `align()` kernel.

`*stats/` | Directory storing cached confusion matrices
---: | ---
`dels_cm.npy` | DELetions confusion matrix.
`inss_cm.npy` | INSertions confusion matrix.
`subs_cm.npy` | SUBsitutions confusion matrix.
`nps_cm.npy` | N-Polymers confusion matrix.

## Data Sources
The GRCh38 human reference sequence was downloaded from: s3://ont-open-data/gm24385_2020.09/config/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
The Genome In A Bottle GRCh38 v4.1 ground truth VCF and benchmarking regions were downloaded from: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh38/
R9.4.1 reads basecalled with the Guppy 5.0.6 super-accuracy model were downloaded from ONT Open Datasets: s3://ont-open-data/gm24385_2020.11/analysis/r9.4.1/20201026_1644_2-E5-H5_PAG07162_d7f262d5/guppy_v4.0.11_r9.4.1_hac_prom/align_unfiltered/chrN/guppy_v5.0.6_r9.4.1_sup_prom/basecalls.fastq.gz

## Acknowledgements
We would like to thank the developers of `samtools`, `minimap2`, `pysam`, `clair3`, `pepper-deepvariant`, `pysankey`, `igv`, and `swalign`. We would also like to thank GIAB and ONT for making their data available publicly.
