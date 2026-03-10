"""
Module: src/extract_fastq.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Stage 1 - Map raw reads to the IS database and extract qualified pairs 
         (soft-clipped or unmapped-mate-mapped) for further analysis.
"""
import sys
import os
import subprocess
import logging
import threading
import pysam
from Bio.Seq import Seq
from src import config

def log_stderr(process):
    """
    Dedicated thread to read stderr from the BWA process to catch errors immediately.
    """
    suppress_keywords = [
        "[M::mem_pestat]", 
        "[M::process]", 
        "[M::mem_process_seqs]", 
        "[M::mem_sort_dedup_and_mate]", 
        "[M::worker]",
        "[M::main_mem]"
    ]
    
    for line in process.stderr:
        line_str = line.decode()
        
        # If the line contains any of the suppression keywords, skip it
        if any(keyword in line_str for keyword in suppress_keywords):
            continue
            
        # Otherwise, write it to stderr (keeps errors visible)
        sys.stderr.write(f"[BWA] {line_str}")

def get_original_sequence(read):
    """
    Recover the original read sequence and quality from the BAM record.
    If the read is reverse complemented (mapped to reverse strand), flip it back
    to restore the original FASTQ sequence orientation.
    """
    seq = read.query_sequence
    qual = read.qual
    
    # If read is mapped to reverse strand, SAM format stores the reverse complement.
    # We need the original sequence for the next stage (FASTQ format).
    if read.is_reverse:
        if seq:
            seq = str(Seq(seq).reverse_complement())
        if qual:
            qual = qual[::-1] # Reverse quality string
        
    return seq, qual

def is_high_n_content(sequence):
    """
    Reject sequence with >10% 'N'.
    """
    if not sequence:
        return True
    return (sequence.count('N') / len(sequence)) > 0.1

def check_soft_clip(read):
    """
    Check if read is soft-clipped with length between MIN_CLIP and MAX_CLIP.
    Checks both left (cigartuples[0]) and right (cigartuples[-1]) clips.
    """
    if not read.cigartuples:
        return False
    # Cigar Op 4 is Soft Clip (BAM specification)
    
    seq_len = read.query_length
    # Sanity check: Read must be long enough to support anchors on both sides
    if seq_len <= (config.MIN_ANCHOR_LEN * 2):
        return False
        
    # Define Dynamic Range
    # Example: 150bp read, Min=20. Clip must be between 20 and 130.
    min_clip = config.MIN_ANCHOR_LEN
    max_clip = seq_len - config.MIN_ANCHOR_LEN
    
    # Check Left Clip (Index 0)
    op_start, len_start = read.cigartuples[0]
    if op_start == 4 and min_clip < len_start < max_clip:
        return True
    # Check Right Clip (Index -1)
    # If both ends are soft-clipped, treat each separately
    op_end, len_end = read.cigartuples[-1]
    if op_end == 4 and min_clip < len_end < max_clip:
        return True
    return False

def align_to_is_db(is_fasta, read_files, threads=16, mode='paired'):
    """
    Execute bwa mem to map reads to IS database.
    Stream processing: Do not write intermediate BAM to disk.
    read_files: List containing either [interleaved_path] or [r1_path, r2_path].
    """
    # Command constructs: bwa mem -t {threads} {is_fasta} {fastq1} {fastq2}
    # Note: Although spec mentions -p, that is for interleaved. 
    # Since we have separate R1/R2, we pass both to standard bwa mem.
    cmd = ["bwa", "mem", "-v", "1", "-t", str(threads), is_fasta] + read_files
    
    if mode == 'interleaved':
        cmd.append("-p")
        logging.info(f"Starting BWA alignment (Interleaved): {' '.join(cmd)}")
    elif mode == 'single':
        # No -p for single end
        logging.info(f"Starting BWA alignment (Single-End): {' '.join(cmd)}")
    else:
        logging.info(f"Starting BWA alignment (Paired): {' '.join(cmd)}")
        
    # Use subprocess.Popen with bufsize management to stream out BWA
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=10**6
    )
    
    # Ensure a dedicated thread reads stderr
    thread = threading.Thread(target=log_stderr, args=(process,))
    thread.daemon = True
    thread.start()
    
    return process

def format_fastq(name, seq, qual, suffix):
    """
    Format entry as standard 4-line FASTQ.
    Suffix (e.g., '/1' or '/2') is appended to the name for clarity in interleaved output.
    """
    return f"@{name}/{suffix}\n{seq}\n+\n{qual}\n"

def process_buffer(reads, out_handle, debug_handle=None, mode='paired'):
    """
    Evaluate a list of reads (pair) and write to output if they meet criteria.
    Logic: Keep pair if soft-clipped OR (unmapped but mate mapped).
    Constraint: Must extract both reads of the pair.
    """
    # Expecting primary alignments for R1 and R2.
    # Filter out secondary/supplementary if strictly focusing on primary.
    
    # --- SINGLE END LOGIC ---
    if mode == 'single':
        count = 0
        for r in reads:
            if r.is_secondary or r.is_supplementary: continue
            
            seq, qual = get_original_sequence(r)
            if is_high_n_content(seq): continue
            
            # Logic: Keep if Soft-Clipped
            # (Unmapped-Mate logic does not exist for Single End)
            if check_soft_clip(r):
                # We append /1 just to keep format consistent, though not strictly required
                out_str = format_fastq(r.query_name, seq, qual, "1")
                out_handle.write(out_str)
                if debug_handle: debug_handle.write(out_str)
                count += 1
        return count > 0
    
    # --- PAIRED END LOGIC ---
    r1 = None
    r2 = None
    
    for r in reads:
        if r.is_secondary or r.is_supplementary:
            continue
        if r.is_read1:
            r1 = r
        elif r.is_read2:
            r2 = r
            
    if not r1 or not r2:
        return # Incomplete pair info
    # 1. Reject sequence with >10% "N"
    seq1, qual1 = get_original_sequence(r1)
    seq2, qual2 = get_original_sequence(r2)
    
    if is_high_n_content(seq1) or is_high_n_content(seq2):
        return
    keep_pair = False
    # 2. Check Logic: Keep pair if...
    # Condition A: Read is soft-clipped (min_clip < len < max_clip)
    if check_soft_clip(r1) or check_soft_clip(r2):
        keep_pair = True
    
    # Condition B: Read is unmapped but mate is mapped
    # Note: If R1 is unmapped, R1.mate_is_unmapped check relies on flag.
    # But since we have both objects, we can check directly.
    if not keep_pair:
        if (r1.is_unmapped and not r2.is_unmapped) or \
           (r2.is_unmapped and not r1.is_unmapped):
            keep_pair = True
    # 3. Output
    if keep_pair:
        # Output qualified pairs (Read 1 + Read 2) for stage 2 alignment.
        # Stage 2 expects FASTQ input (interleaved or handled by stream).
        out_str1 = format_fastq(r1.query_name, seq1, qual1, "1")
        out_str2 = format_fastq(r2.query_name, seq2, qual2, "2")
        
        # Write to the provided file handle instead of sys.stdout
        out_handle.write(out_str1)
        out_handle.write(out_str2)
        
        # If --debug-fastq is turned on
        if debug_handle:
            debug_handle.write(out_str1)
            debug_handle.write(out_str2)
        return True
    return False

def run_extraction(is_db, read_files, out_handle, threads=16, debug_path=None, mode='paired'):
    """
    Main Module Function:
    1. Runs BWA MEM via subprocess.
    2. Streams BAM output.
    3. Filters reads and writes them to 'out_handle'.
    """
       
    # Setup logging
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # Start BWA stream
    bwa_process = align_to_is_db(is_fasta=is_db, read_files=read_files, threads=threads, mode=mode)
    
    # Open pysam AlignmentFile on stdin
    try:
        # "r" mode reads SAM/BAM. Since BWA outputs SAM by default to stdout, this works.
        bam_stream = pysam.AlignmentFile(bwa_process.stdout, "r")
    except ValueError as e:
        logging.error(f"Failed to open BAM stream: {e}")
        # Ensure we kill the process if stream fails
        bwa_process.kill()
        raise e

    # Optional debug file
    debug_fh = None
    if debug_path:
       debug_fh = open(debug_path, "w")
       logging.info(f"Debug FASTQ will be saved to: {debug_path}")
    
    current_qname = None
    read_buffer = []

    # Counter for statistics
    extracted_count = 0
    processed_count = 0

    for read in bam_stream:
        processed_count += 1
        
        # Update screen every 100,000 reads
        if processed_count % 1000000 == 0:
        # The newline \n is critical so log_stream can read it as a complete line
            sys.stderr.write(f"PROGRESS: Processed {processed_count:,} reads...\n")
            sys.stderr.flush()	
        
        # Grouping logic
        if current_qname is None:
            current_qname = read.query_name
            read_buffer.append(read)
        elif read.query_name == current_qname:
            read_buffer.append(read)
        else:
            # New query name encountered; process previous buffer
            if process_buffer(read_buffer, out_handle, debug_fh, mode=mode):
                extracted_count += 1	
            # Reset buffer
            current_qname = read.query_name
            read_buffer = [read]
    # Process the last buffer
    if read_buffer:
        if process_buffer(read_buffer, out_handle, debug_fh, mode=mode):
            extracted_count += 1
    
    # Log the summary
    logging.info(f"Stage 1 Complete: Processed {processed_count} reads.Extracted {extracted_count} candidate read pairs. ")    

    # Cleanup
    bam_stream.close()
    if debug_fh:
        debug_fh.close()
    
    # Check BWA return code
    bwa_process.stdout.close()
    return_code = bwa_process.wait()
    if return_code != 0:
        logging.error(f"Align to ISs failed with return code {return_code}")
        sys.exit(return_code)

