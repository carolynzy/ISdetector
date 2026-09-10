# ISdetector
**ISdetector** is a Python package designed to detect Insertion Sequences (IS) and associated Structural Variations (SVs) using paired-end or single-end Whole Genome Sequencing (WGS) data. It employs a hybrid strategy utilizing both discordant read pairs and split-read signals to identify novel insertion sites, known (reference) sites, and complex deletions.
## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Dataset](#dataset)
- [Arguments](#arguments)
- [Outputs](#outputs)

## Prerequisites

For the recommended installation, install **Conda or Mamba** and make it available in `PATH`. Run `install.sh` from Bash. The installer downloads Python and the tools below into a dedicated environment; you do not need to install them separately.

### External Tools

The pipeline uses:
* **BWA** (Burrows-Wheeler Aligner)
* **Samtools >=1.3** (ISdetector checks the version on startup)
* **BLASTn** (NCBI BLAST+)
### Python Dependencies

ISdetector requires Python >=3.10. The installer provisions a compatible Python and installs these dependencies:
* `pysam` (>=0.19.0)
* `biopython` (>=1.80)
* `pandas`
* `numpy`
## Installation
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/carolynzy/isdetector.git
    cd isdetector
    ```
2.  **Create a dedicated environment and install:**
    ```bash
    ./install.sh
    ```
    Conda or Mamba must already be installed. The installer always creates a new environment containing Python, BWA, BLAST+, Samtools >=1.3, and ISdetector. Existing system tools and environments are preserved. It validates the environment's tools and stores their absolute paths privately in `<environment>/etc/isdetector/tools.json`, outside the repository. Runtime commands use these saved paths and report an error if a tool disappears or becomes incompatible.

    Alternatively, for development, choose an editable installation instead of the command above so local Python source changes take effect on the next run:
    ```bash
    ./install.sh --editable
    ```
    Keep the repository at the same location for editable installations. Editable mode does not automatically download code updates or upgrade dependencies.
    The default environment is `~/.local/share/isdetector/env`. Use `--prefix /path/to/new/environment` to choose another location. Existing locations are rejected. To reinstall, choose a new prefix; the existing managed launcher is updated only after successful installation. Environments and launchers must be outside this repository.

    The launcher is installed at `~/.local/bin/isdetector` (customizable with `--launcher-dir`). With that directory first in `PATH`, simply run `isdetector` to show help, or `isdetector <arguments>` to analyze data. The launcher runs inside the dedicated environment automatically; manual `conda activate` is unnecessary. If needed, add `export PATH="$HOME/.local/bin:$PATH"` to your shell startup file. The launcher does not change the calling shell's active environment.

    Direct `python -m pip install .` remains supported for users who manage dependencies themselves. Without a saved configuration, ISdetector searches `PATH` for compatible tools at startup. Direct pip installation does not create the environment or launcher.

3.  **Make the launcher available and check installation:**

    ```bash
    export PATH="$HOME/.local/bin:$PATH"
    isdetector
    ```

    Add the `export` line to your shell startup file (for example, `~/.bashrc`) to keep it for future terminals. Running the launcher without arguments displays help. You can also invoke it directly as `~/.local/bin/isdetector`.

To reinstall into a new environment while preserving the old one:

```bash
./install.sh --prefix "$HOME/.local/share/isdetector/env-new"
```

After a successful installation, the managed launcher points to the new environment. The installer refuses to overwrite an unrelated command at the launcher location. Use `./install.sh --help` to see all options, including `--launcher-dir`.

## Usage

Run `isdetector` from any working directory once its launcher directory is in `PATH`. The launcher automatically runs the command in the dedicated environment, using the saved tool paths; no manual activation is needed.

**For Paired-End Reads:**
```bash
isdetector \
    -1 data/reads_R1.fastq.gz \
    -2 data/reads_R2.fastq.gz \
    -i data/is_elements.fasta \
    -r data/reference_genome.gb \
    -s MySample \
    -o ./results \
    -t 16
```
## Workflow

<p align="center">
    <img src="flowchart_v2.png" alt="ISdetector Workflow" width="300">
</p>

## Dataset 
The test data for this pipeline is archived on Zenodo:
**DOI:** [https://doi.org/10.5281/zenodo.18996276](https://doi.org/10.5281/zenodo.18996276)

It could also be downloaded by running: 
```bash
./download_data.sh
```
which will create a folder "data" and download the test data automatically.

## Arguments

```bash
  -h, --help            show this help message and exit
  -i IS_DB, --is_db IS_DB
                        FASTA file of IS sequences (Bait).
  -r REFERENCE, --reference REFERENCE
                        Reference Genome (.fasta or .gb/.gbk).
  -s SAMPLE, --sample SAMPLE
                        Prefix used for output files and log files.
  -o OUTDIR, --outdir OUTDIR
                        Directory for all outputs.
  -t THREADS, --threads THREADS
                        CPU threads (default: 16).
  --debug-fastq         Save Stage 1 extracted reads to file.
  --debug-signal        Save Stage 3 insertion signals to file.
```

Input Reads (Choose One Strategy):
```bash
  -1 FASTQ1, --fastq1 FASTQ1
                        Path to Raw Read 1 (FASTQ/GZ).
  -2 FASTQ2, --fastq2 FASTQ2
                        Path to Raw Read 2 (FASTQ/GZ).
  -f FASTQ, --fastq FASTQ
                        Path to Interleaved Paired-end FASTQ/GZ.
  -u UNPAIRED, --unpaired UNPAIRED
                        Path to Single-End FASTQ/GZ.
```
## Outputs

Each IS analysis writes `SAMPLE_ISNAME_report.tsv`. When no peaks are found, this file contains `No peaks detected.` instead of a table. Nonempty reports contain these columns, in order:

| Column | Description |
| :--- | :--- |
| Chromosome | Chromosome or contig in the original reference. |
| Position | Representative junction coordinate on the original reference, using zero-based boundary coordinates. |
| Gap | Absolute distance between paired peaks for non-SV calls; `NA` for SV or singular calls. |
| Orientation | IS orientation relative to the reference: `F` (forward) or `R` (reverse). |
| IS_Length | Absolute difference between the paired peaks' IS coordinates; `NA` for singular calls. |
| Category | `Known`, `Known(SV)`, `Novel`, or `Novel(SV)`. |
| Start_clipped | Number of right-clipped signals supporting the event. |
| End_clipped | Number of left-clipped signals supporting the event. |
| Discordant | Number of non-proper paired alignment records in the representative peak's cluster. Single-end records are also counted here. |
| Proper | Number of proper-pair alignment records in the representative peak's cluster. |
| SV | Deletion interval, such as `Deletion_150_300`, when a deletion pair is identified; otherwise `None`. |

A position ending in an asterisk (for example, `12345*`) represents a singular peak with no pairable counterpart. Support counts refer to alignment records or signals, not necessarily distinct read pairs. An `(SV)` category reflects a peak's depth-based status and can occur without an explicit deletion interval in `SV`.

When a GenBank reference is provided and annotations are generated, `SAMPLE_ISNAME_annotations.tsv` contains:

| Column | Description |
| :--- | :--- |
| Group_ID | Identifier of the detected event grouping one or more peaks. |
| Position | Representative event coordinate on the original reference. |
| SV_Type | Deletion interval if available; otherwise `None`. |
| Annotation_Type | `Intragenic`, `Intergenic`, `Exonic (Deletion)`, or `Intergenic (Deletion)`. Annotation examines CDS features. |
| Gene | Gene name or locus tag for a CDS hit; affected genes for a deletion; or nearest upstream/downstream CDS labels and distances for an intergenic insertion. |
| Protein | Duplicates the `Gene` field; it is not a separate protein-name annotation. |
| Product | CDS product qualifier(s), or `NA` when unavailable. |

With `--debug-fastq`, the combined extracted reads are retained in `SAMPLE_extracted_combined.fq`, and an additional copy is written to `SAMPLE_debug_extracted.fq`. Paired-end output is interleaved; single-end output contains individual reads.

With `--debug-signal`, when peaks and signals are present, `SAMPLE_ISNAME_signals.tsv` contains:

| Column | Description |
| :--- | :--- |
| Chromosome | Contig name in the IS-clean reference, including its `_clean` suffix. |
| Cluster_id | Cluster identifier in the form `{chromosome}_{start}_{end}`, using IS-clean reference coordinates. |
| Read_id | Read name from the alignment. |
| is_read1 | `True` for read 1, `False` for read 2; single-end reads are treated as read 1. |
| Position | Zero-based junction boundary on the IS-clean reference. Debug signals are not lifted back to original-reference coordinates. |
| IS_coordinate | Zero-based alignment boundary on the IS sequence. |
| Clipped_sequence | Soft-clipped sequence that matched the IS. |
| IS_strand | `1` for a forward match or `-1` for a reverse-complement match to the IS. |
| Orientation | `F` (forward) or `R` (reverse). |
| Clip_side | `LEFT` or `RIGHT`, identifying the soft-clipped end of the alignment. |

Pipeline progress and errors are written to `SAMPLE.log` and the screen. Intermediate alignment BAM files are removed after successful per-IS processing.
