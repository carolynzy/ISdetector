"""
Module: src/main.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Main entry point and pipeline controller.
         1. Parses command-line arguments.
         2. Pre-processes reference files (GenBank to FASTA, Clean Reference, Indexing).
         3. Loop per IS:
            a. Checks if output exists (skips if yes).
            b. Creates specific IS-clean reference.
            c. Aligns "Bait" reads to specific reference.
            d. Detects, Annotates, Reports.
            e. Catches errors and proceeds to next IS.
            f. Generates final insertion and signal reports.
Usage:
    isdetector -1 R1.fq -2 R2.fq -i is.fa -r ref.gb -s sample_id -o ./output
"""	
import os
import sys
import argparse
import subprocess
import logging
import gzip
import threading
import pysam
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from src import utils, extract_fastq, cluster, detector, config

# Specific statistics tracking dictionary
stats = {
    "processed_reads": 0,
    "extracted_pairs": 0
}

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="ISdetector: Detect Insertion Sequences and SVs.")
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
    dependencies = ["bwa", "samtools", "blastn"]
    for dep in dependencies:
        if subprocess.call(f"which {dep}", shell=True, stdout=subprocess.DEVNULL) != 0:
            logging.error(f"{dep} not found in PATH.")
            sys.exit(1)

def log_stream(stream, log_level, parse_extraction_count=False):
    """
    Reads from a stream line by line and logs it.
    """
    # Define messages to suppress
    suppress_keywords = ["[M::mem_pestat]", "[M::process]", "[M::mem_process_seqs]"]
    
    for line in stream:
        try:
            line_str = line.decode().strip()
                
            if not line_str:
                continue
            
            # --- Feature: BWA Suppression ---
            if any(keyword in line_str for keyword in suppress_keywords):
                continue  # Skip this line
            
            # --- Feature: Real-time Progress Display ---
            # If the subprocess sends a progress update, print it with \r to overwrite 
            # the current line on the console. Do NOT send to log_level (file).
            if line_str.startswith("PROGRESS:"):
                sys.stderr.write(f"\r{line_str}")
                sys.stderr.flush()
                continue
            
            # --- Standard Logging ---
            # If we were printing progress bars, clear the line first to prevent ghost text
            if "\r" in sys.stderr.getvalue() if hasattr(sys.stderr, 'getvalue') else True:
                 sys.stderr.write("\r" + " " * 80 + "\r") # Clear line
                         
            log_level(line_str)
                
            if parse_extraction_count and ("Extracted" in line_str or "Processed" in line_str) and ("read pairs" in line_str or "reads" in line_str):
                try:
                    parts = line_str.split()
                    if "Extracted" in parts:
                        idx = parts.index("Extracted")
                        stats["extracted_pairs"] = int(parts[idx + 1])
                    if "Processed" in parts:
                        idx = parts.index("Processed")
                        stats["processed_reads"] = int(parts[idx + 1])
                except Exception: 
                    pass
        except Exception: 
            pass

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
    check_dependencies()
    
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
            subprocess.check_call(["bwa", "index", global_is_path], 
                                  stdout=subprocess.DEVNULL, 
                                  stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
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
            # We pass the specific IS seq to mask it in the reference
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
            cmd_sort = ["samtools", "sort", "-@", str(args.threads), "-"]
            
            logging.info(f"Stage 2: Aligning to {os.path.basename(ref_clean_path)}...")
            
            with open(os.devnull, "w") as devnull:
                p2_bwa = subprocess.Popen(cmd_bwa, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p2_sort = subprocess.Popen(cmd_sort, stdin=p2_bwa.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # Start logging threads
                t_bwa = threading.Thread(target=log_stream, args=(p2_bwa.stderr, logging.debug)) # Debug level to reduce noise
                t_sort = threading.Thread(target=log_stream, args=(p2_sort.stderr, logging.debug))
                t_bwa.start(); t_sort.start()
                
                p2_bwa.stdout.close() # Allow p2_bwa to receive SIGPIPE if sort closes
                
                # 4. Stage 3: Clustering
                # Read from Sort stream, write to disk BAM
                bam_stream = pysam.AlignmentFile(p2_sort.stdout, "rb")
                outfile_bam = pysam.AlignmentFile(bam_path, "wb", template=bam_stream)
                
                cluster_dict = cluster.create_microclusters(
                    stream_and_save(bam_stream, outfile_bam)
                )
                
                outfile_bam.close()
                bam_stream.close()
                p2_sort.wait()
                p2_bwa.wait()
                t_bwa.join(); t_sort.join()

            if p2_sort.returncode != 0:
                raise Exception("Alignment/Sorting failed.")
                
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
            logging.error(f"Pipeline failed for IS {is_name}: {e}")
            logging.error("Proceeding to next IS...")
            continue

    # Final Cleanup
    if os.path.exists(extracted_fq_path) and not args.debug_fastq:
        os.remove(extracted_fq_path)
    if args.reference.endswith(('.gb', '.gbk')) and os.path.exists(ref_fasta_original):
        os.remove(ref_fasta_original)

    logging.info(f"All IS sequences processed.\n{'='*60}\n{'='*60}")

if __name__ == "__main__":
    main()
