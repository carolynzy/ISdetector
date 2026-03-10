**1. Software Introduction**

ISdetector is a command-line-based bioinformatics software tool designed to detect Insertion Sequence (IS) related events from high-throughput sequencing data. It outputs insertion sites, orientations, supporting evidence, and (optional) annotation results. The software uses a pipeline-style processing workflow:

•	Stage 1: Read extraction 

•	Stage 2: Alignment and sorting 

•	Stage 3: Clustering and peak detection 

•	Stage 4: Insertion event identification/classification 

•	Stage 5: Annotation and report output 


**2. System Requirements**

•	Supported Operating Systems: 

Linux (Ubuntu, Debian, or CentOS series recommended). It can also run on Windows via WSL2 (Ubuntu) provided external dependencies are available and the PATH is configured correctly.

•	Hardware Recommendations:

    o	CPU: 4 cores or more (8–32 cores recommended for batch samples).

    o	Memory: ≥16 GB (≥32 GB recommended for batch or large datasets).

    o	Storage: SSD, with 50–200 GB+ reserved based on the scale of FASTQ/BAM/output files.


**3. System Requirements**

•	Supported Operating Systems: Linux (Ubuntu, Debian, or CentOS series recommended). It can also run on Windows via WSL2 (Ubuntu) provided external dependencies are available and the PATH is configured correctly.

•	Hardware Recommendations:

    o	CPU: 4 cores or more (8–32 cores recommended for batch samples).

    o	Memory: ≥16 GB (≥32 GB recommended for batch or large datasets).

    o	Storage: SSD, with 50–200 GB+ reserved based on the scale of FASTQ/BAM/output files.

**3. Software Dependencies and Installation**

**3.1 External Dependencies (Required)**
ISdetector checks for the following commands in your PATH before running:

•	bwa: Used for aligning the IS library and the reference genome.

•	samtools: Used for sorting.

•	blastn (NCBI BLAST+): Used to locate and "clean" segments highly similar to the IS on the reference genome to generate a "clean reference".

Installation on Ubuntu/Debian:
```Bash
sudo apt-get update
sudo apt-get install -y bwa samtools ncbi-blast+
```
**3.2 Python Environment (Required)**

•	Recommended Python Version: 3.8+ (3.9, 3.10, or 3.11 are all acceptable).

•	Required Python Packages: pysam, pandas, numpy, and biopython.

Example Installation using Conda:
```Bash
conda create -n isdetector python=3.10 -y
conda activate isdetector
pip install pysam pandas numpy biopython
```
**3.3 Software Installation**
Run the following in the project root directory:
```Bash
cd ISdetector
pip install -e .
```

**4. Input File Preparation**

ISdetector supports two types of read inputs (choose one):

1.	Paired-end FASTQ: -1 R1.fq(.gz) and -2 R2.fq(.gz).
  
2.	Interleaved FASTQ: -f interleaved.fq(.gz).

Other Required Files:

•	IS Sequence Library (FASTA): -i is_db.fasta.

•	Reference Genome (FASTA or GenBank): -r ref.fasta or -r ref.gbk/.gb.

•	Sample ID: -s sample_id.

•	Output Directory: -o outdir.

**5. Running the Program**

**5.1 View Help**
```Bash
isdetector --help
```
**5.2 Typical Execution Examples**

•	Paired FASTQ:
```
python -m src.main -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -i IS6110.fasta -r H37Rv.gbk -s sample001 -o ./out_sample001 -t 16
```

•	Interleaved FASTQ:
```
python -m src.main -f sample_interleaved.fastq.gz -i IS6110.fasta -r H37Rv.fasta -s sample001 -o ./out_sample001 -t 16
```

**6. Parameter Settings and Tuning**

Parameter	Required	Descriptio

```
-1 / --fastq1	One of two	Read1 (FASTQ/GZ)
-2 / --fastq2	One of two	Read2 (FASTQ/GZ)
-f / --fastq	One of two	Interleaved FASTQ/GZ
-i / --is_db	Yes	IS sequence library (FASTA)
-r / --reference	Yes	Reference genome (FASTA or GenBank)
-s / --sample	Yes	Sample prefix (used for logs and output naming)
-o / --outdir	Yes	Output directory
-t / --threads	No	Number of threads (default: 16)
--debug-fastq	No	Save candidate reads extracted in Stage 1
--debug-signal	No	Output Stage 3 signal details (TSV)
```
