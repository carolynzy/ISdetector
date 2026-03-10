# ISdetector
**ISdetector** is a Python package designed to detect Insertion Sequences (IS) and associated Structural Variations (SVs) using paired-end or single-end Whole Genome Sequencing (WGS) data. It employs a hybrid strategy utilizing both discordant read pairs and split-read signals to identify novel insertion sites, known (reference) sites, and complex deletions.
## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Arguments](#arguments)
- [Output](#output)
- [Workflow](#workflow)
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
    git clone [https://github.com/carolynzy/isdetector.git](https://github.com/carolynzy/isdetector.git)
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
