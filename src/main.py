"""
Module: src/main.py
Project: ISdetector
Version: 1.1.0
Author: Yang Zhou
Purpose: Parse arguments and validate external tools, extract candidate reads once,
         then prepare a clean reference, align, detect, annotate, and report per IS.
         Skip existing reports and continue after individual IS failures; return a
         nonzero exit status if any IS analysis fails.
Usage:
    isdetector -1 R1.fq -2 R2.fq -i is.fa -r ref.gb -s sample_id -o ./output
"""	
import os
import sys
import argparse
import subprocess
import logging
from src.dependencies import check_samtools_version, configured_executable
from src.processes import run_checked, streaming_process
import pysam
import pandas as pd
from Bio import SeqIO
from src import __version__, utils, extract_fastq, cluster, detector

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=f"ISdetector v{__version__}: Detect Insertion Sequences and SVs."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program version and exit.",
    )
    # Input Files
    input_group = parser.add_argument_group('Input Reads (Choose One Strategy)')
    input_group.add_argument("-1", "--fastq1", help="Path to Raw Read 1 (FASTQ/GZ).")
    input_group.add_argument("-2", "--fastq2", help="Path to Raw Read 2 (FASTQ/GZ).")
    input_group.add_argument("-f", "--fastq", help="Path to Interleaved Paired-end FASTQ/GZ.")
    input_group.add_argument("-u", "--unpaired", help="Path to Single-End FASTQ/GZ.")
    
    parser.add_argument("-i", "--is_db", required=True, help="FASTA file of IS sequences (Bait).")
    parser.add_argument("-r", "--reference", required=True, help="Reference Genome (.fasta or .gb/.gbk).")
    
    # Output and naming
    parser.add_argument("-s", "--sample", required=True, help="Prefix used for output files and log files.")
    parser.add_argument("-o", "--outdir", required=True, help="Directory for all outputs.")
    
    # Optional / Config
    parser.add_argument("-t", "--threads", type=int, default=16, help="CPU threads (default: 16).")
    
    # Debug flags
    parser.add_argument("--debug-fastq", action="store_true", help="Save Stage 1 extracted reads to file.")
    parser.add_argument("--debug-signal", action="store_true", help="Save Stage 3 insertion signals to file.")

    args = parser.parse_args()

    modes = []
    if args.fastq1 and args.fastq2: modes.append("paired")
    if args.fastq: modes.append("interleaved")
    if args.unpaired: modes.append("single")

    if len(modes) > 1:
        parser.error("Conflicting inputs: Choose ONE of paired (-1/-2), interleaved (-f), or unpaired (-u).")
    if len(modes) == 0:
        parser.error("Missing inputs: No read files provided.")
    
    # Store the mode in args for easy access later
    args.mode = modes[0]

    return args

def setup_logging(outdir, sample):
    log_file = os.path.join(outdir, f"{sample}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info(f"Logging initialized. Output directory: {outdir}")

def check_dependencies():
    try:
        for dep in ("bwa", "blastn"):
            logging.info("Using %s at %s", dep, configured_executable(dep))
        samtools_path = configured_executable("samtools")
        version = check_samtools_version(samtools_path)
    except RuntimeError as exc:
        logging.error("%s", exc)
        sys.exit(1)
    logging.info("Using %s at %s", version, samtools_path)
    return samtools_path

def stream_and_save(input_bam, output_bam_handle):
    """Yields reads from stream and saves to disk simultaneously."""
    for read in input_bam:
        output_bam_handle.write(read)
        yield read

def main():
    args = parse_args()
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        
    setup_logging(args.outdir, args.sample)
    samtools_path = check_dependencies()
    
    # -------------------------------------------------------------------------
    # 0. Load Data & Prepare Global Reference
    # -------------------------------------------------------------------------
    logging.info(f"Loading IS Database: {args.is_db}")
    try:
        is_records = list(SeqIO.parse(args.is_db, "fasta"))
    except Exception as e:
        logging.error(f"Failed to parse IS database: {e}")
        sys.exit(1)

    if not is_records:
        logging.error("No records found in IS database.")
        sys.exit(1)
    
    is_records_toprocess = []
        
    for record in is_records:
        is_name = record.id
        is_seq = str(record.seq)
        is_len = len(record.seq)
        safe_is_name = is_name.replace("|", "_").replace("/", "_")
        
        # Define output filename to check existence
        report_path = os.path.join(args.outdir, f"{args.sample}_{safe_is_name}_report.tsv")
        
        # Check if output exists
        if os.path.exists(report_path):
            logging.warning(f"Output for {is_name} already exists. Skipping.")
            continue
        is_records_toprocess.append(record)
    
    if len(is_records_toprocess) == 0:
        logging.info(f"{'=' * 60}")
        logging.info("All ISs are already processed.")
        logging.info(f"{'=' * 60}")
        sys.exit(1)    
    else:
        logging.info(f"Starting analysis for {len(is_records_toprocess)} IS records...")    
    
    # We write the loaded records to a clean temporary file in the output dir.
    # This ensures BWA can read it (fixing any format issues in the original file).    
    global_is_path = os.path.join(args.outdir, "is_database_combined.fasta")
    SeqIO.write(is_records_toprocess, global_is_path, "fasta")
    
    # Index the IS Database (Required for Stage 1 Extraction)
    if not os.path.exists(global_is_path + ".bwt"):
        logging.info(f"Indexing IS Database: {global_is_path}")
        try:
            run_checked(["bwa", "index", global_is_path])
        except (subprocess.CalledProcessError, OSError):
            logging.error("Failed to index IS Database. Check 'bwa' installation.")
            sys.exit(1)
    
    # -------------------------------------------------------------------------
    # 0b. Prepare Reference Genome
    # -------------------------------------------------------------------------       
    # Handle Reference Format & Cleaning
    ref_fasta_original = args.reference
    gb = None

    if args.reference.endswith(('.gb', '.gbk')):
        logging.info("Processing GenBank file...")
        gb = utils.Genebank(args.reference)
        # Extract FASTA from GenBank
        ref_temp_path = os.path.join(args.outdir, "ref_original.fasta")
        ref_fasta_original = gb.convert_to_fasta(ref_temp_path)
    elif args.reference.endswith(('.fasta', '.fa')):
        gb = None
    else:
        logging.error("Invalid reference format. Must be .fasta or .gb/.gbk")
        sys.exit(1)
        
    # -------------------------------------------------------------------------
    # Stage 1: Global Extraction (Run ONCE for ALL IS)
    # -------------------------------------------------------------------------
    # We extract reads that match ANY sequence in the IS database.
    # We MUST save these to a file to reuse them for each IS loop.
    extracted_fq_path = os.path.join(args.outdir, f"{args.sample}_extracted_combined.fq")

    logging.info("Stage 1: Global Extraction (Aligning to all IS sequences)...")
    
    # Determine input files (List format)
    if args.mode == "interleaved":
        input_reads = [args.fastq]
    elif args.mode == "single":
        input_reads = [args.unpaired]
    else: # paired
        input_reads = [args.fastq1, args.fastq2]
        
    # Stage 1: Extraction
    # Run Extraction Module directly
    try:
        with open(extracted_fq_path, "w") as f_out:
            extract_fastq.run_extraction(
                is_db=global_is_path,     # Use filtered DB
                read_files=input_reads,    # Pass list of files
                out_handle=f_out,          # Pass file handle
                threads=args.threads,
                debug_path=os.path.join(args.outdir, f"{args.sample}_debug_extracted.fq") if args.debug_fastq else None,
                mode=args.mode
            )
    except Exception as e:
        logging.error(f"Stage 1 Extraction failed: {e}")
        sys.exit(1)
        
    logging.info(f"Stage 1 Complete. Extracted reads saved to: {extracted_fq_path}")
        
    # -------------------------------------------------------------------------
    # Loop: Process Each IS Individually
    # -------------------------------------------------------------------------    
    failed_is = []
    for record in is_records_toprocess:
        is_name = record.id
        is_seq = str(record.seq)
        is_len = len(record.seq)
        safe_is_name = is_name.replace("|", "_").replace("/", "_")
        
        # Define output filename to check existence
        report_path = os.path.join(args.outdir, f"{args.sample}_{safe_is_name}_report.tsv")
            
        logging.info(f"{'=' * 60}")
        logging.info(f"Processing IS: {is_name} ({is_len} bp)")
        logging.info(f"{'=' * 60}")
        
        try:
            # 2. Generate IS-Specific Clean Reference
            # Remove matches to this IS from the reference and retain coordinate mapping.
            ref_clean_path, crossmap = utils.generate_clean_reference(
                ref_fasta_original, is_seq, is_len, args.outdir, suffix=f"_{safe_is_name}"
            )
            
            # Index the specific reference
            utils.indexing(ref_clean_path) 
            
            # 3. Stage 2: Specific Alignment (Combined Reads -> Specific Clean Ref)
            bam_path = os.path.join(args.outdir, f"{args.sample}_{safe_is_name}.bam")
            
            cmd_bwa = ["bwa", "mem", "-t", str(args.threads), "-v", "1"]
            if args.mode != "single":
                cmd_bwa.append("-p") # Interleaved input flag for BWA (Extracted reads are interleaved if paired)
            cmd_bwa.extend([ref_clean_path, extracted_fq_path])
            cmd_sort = [samtools_path, "sort", "-@", str(args.threads), "-"]
            
            logging.info(f"Stage 2: Aligning to {os.path.basename(ref_clean_path)}...")
            
            with streaming_process(cmd_bwa) as p2_bwa:
                with streaming_process(cmd_sort, stdin=p2_bwa.stdout) as p2_sort:
                    p2_bwa.stdout.close()
                    with pysam.AlignmentFile(p2_sort.stdout, "rb") as bam_stream:
                        with pysam.AlignmentFile(bam_path, "wb", template=bam_stream) as outfile_bam:
                            cluster_dict = cluster.create_microclusters(
                                stream_and_save(bam_stream, outfile_bam)
                            )

            # Index the specific BAM
            pysam.index(bam_path)
            
            # 5. Peak Detection
            logging.info(f"Stage 3: Detecting peaks for {is_name}...")
            signal_list, peak_list = cluster.scan_clusters_for_peaks(
                cluster_dict, is_seq, bam_path, crossmap
            )

            logging.info(f"Stage 3: Detected {len(peak_list)} peaks for {is_name}.")
            
            if not peak_list:
                logging.info(f"No peaks found for {is_name}.")
                # Create empty report
                with open(report_path, "w") as f:
                    f.write("No peaks detected.\n")
            else:
                # 6. Stage 4: Structure & Reporting
                logging.info(f"Stage 4: Detecting SVs for {is_name}...")
                final_results, groups, group_cores = detector.detect_insertions(cluster_dict, peak_list, is_len)
                
                logging.info(f"Stage 4: Detected {len(final_results)} insertions for {is_name}.")
                
                # Write Structural Report
                if len(final_results) > 0: 
                    utils.write_insertion_report(final_results, report_path)
                    logging.info(f"Report saved: {report_path}")

                # 7. Annotation
                if gb:
                    ann = utils.generate_annotation_report(final_results, gb)
                    
                    if ann:
                        df = pd.DataFrame(ann)
                        # Sort by Group_ID and Position
                        df.sort_values(by=["Group_ID", "Position"], inplace=True)
                        annot_path = f"{args.outdir}/{args.sample}_{safe_is_name}_annotations.tsv"
                        df.to_csv(annot_path, sep="\t", index=False)
                        
                # Debug Signals
                if args.debug_signal and signal_list:
                     utils.write_signal_report(signal_list, f"{args.outdir}/{args.sample}_{safe_is_name}_signals.tsv")

            # Cleanup for this IS and specific reference files
            if os.path.exists(bam_path): 
                os.remove(bam_path)
            if os.path.exists(bam_path + ".bai"): 
                os.remove(bam_path + ".bai")

            if os.path.exists(ref_clean_path): 
                os.remove(ref_clean_path)
            if os.path.exists(global_is_path): 
                os.remove(global_is_path)
            for ext in [".bwt", ".pac", ".ann", ".amb", ".sa", ".bai"]:
                if os.path.exists(ref_clean_path + ext): 
                    os.remove(ref_clean_path + ext)
                if os.path.exists(global_is_path + ext):
                    os.remove(global_is_path + ext)

        except Exception as e:
            failed_is.append(is_name)
            logging.error(f"Pipeline failed for IS {is_name}: {e}")
            logging.error("Proceeding to next IS...")
            continue

    # Final Cleanup
    if os.path.exists(extracted_fq_path) and not args.debug_fastq:
        os.remove(extracted_fq_path)
    if args.reference.endswith(('.gb', '.gbk')) and os.path.exists(ref_fasta_original):
        os.remove(ref_fasta_original)

    if failed_is:
        logging.error("Analysis failed for %d IS sequences: %s", len(failed_is), ", ".join(failed_is))
        sys.exit(1)

    logging.info(f"All IS sequences processed.\n{'='*60}\n{'='*60}")

if __name__ == "__main__":
    main()
