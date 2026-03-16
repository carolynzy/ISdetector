# ISdetector
**ISdetector** is a Python package designed to detect Insertion Sequences (IS) and associated Structural Variations (SVs) using paired-end or single-end Whole Genome Sequencing (WGS) data. It employs a hybrid strategy utilizing both discordant read pairs and split-read signals to identify novel insertion sites, known (reference) sites, and complex deletions.
## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Dataset](#dataset)
- [Arguments](#arguments)
- [Output](#output)

## Prerequisites
### External Tools
The following tools must be installed and available in your system's `PATH`:
* **BWA** (Burrows-Wheeler Aligner)
* **Samtools**
* **BLASTn** (NCBI BLAST+)
### Python Dependencies
* Python 3.x
* `pysam` (>=0.19.0)
* `biopython` (>=1.80)
* `pandas`
* `numpy` (>=1.18.0)
## Installation
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/carolynzy/isdetector.git
    cd isdetector
    ```
2.  **Install the package and dependencies globally:**
    ```bash
    pip install .
    ```
    *(Note: Using `pip install .` automatically installs the dependencies listed in `setup.py` and registers `isdetector` as a command-line tool.)*
## Usage
Because the tool is installed via `setup.py`, you can run it directly using the `isdetector` command from anywhere.
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

The pipeline generates a final results file (typically final_results.txt or .csv) with the following columns:

| Column |  Description |
| :--- | :--- |
|Chromosome|The genomic scaffold/chromosome where the event was detected.|
| Position | The specific genomic coordinate of the insertion junction. | 
| Category| Classification of the hit: Known, Known(SV), Novel, or Novel(SV).| 
| Orientation| The direction of the IS element relative to the genome (+ or -).| 
| Gap| The distance between the paired peaks (represents deletion size if applicable).| 
| IS_Length| The total length of the IS element detected.| 
| Start/End_Clipped| The count of supporting soft-clipped reads at the start and end junctions.| 
| Discordant_Count| Number of paired-end reads where mates map to different locations (supporting evidence).| 
| SV_Type| Description of associated structural variants (e.g., Deletion_150_300).| 
