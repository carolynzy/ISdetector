"""
Module: src/extract_fastq.py
Project: ISdetector
Version: 1.1.0
Author: Yang Zhou
Purpose: Stage 1 - Map raw reads to the IS database and extract qualified pairs 
         (soft-clipped or unmapped-mate-mapped) for further analysis.
"""
import sys
import logging
from contextlib import ExitStack
import pysam
from Bio.Seq import Seq
from src import config
from src.processes import streaming_process

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
    Check for a soft clip strictly between MIN_ANCHOR_LEN and
    the read length minus MIN_ANCHOR_LEN.
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
    Return a managed BWA process that streams SAM alignments to the IS database.
    read_files contains one single-end/interleaved file or two paired-end files.
    """
    # BWA uses -p for interleaved input; separate paired files are positional inputs.
    cmd = ["bwa", "mem", "-a", "-Y", "-v", "1", "-t", str(threads)]
    
    if mode == 'interleaved':
        cmd.append("-p")
    cmd.extend([is_fasta] + read_files)
    logging.info("Starting BWA alignment (%s): %s", mode, ' '.join(cmd))
    return streaming_process(cmd)

def format_fastq(name, seq, qual, suffix):
    """
    Format entry as standard 4-line FASTQ.
    Suffix (e.g., '1' or '2') is appended after '/'  to the name for clarity in interleaved output.
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

def process_is_buffer(reads, out_handles, mode='paired'):
    """Route each query to every IS for which the original criteria hold.

    Alternative alignments establish IS membership; primary records supply the
    original FASTQ sequences. Mate mapping is evaluated relative to each IS.
    """
    primary = {}
    by_is = {}
    for read in reads:
        mate = 1 if mode == 'single' or read.is_read1 else 2
        if not read.is_secondary and not read.is_supplementary:
            primary[mate] = read
        if not read.is_unmapped and read.reference_name in out_handles:
            by_is.setdefault(read.reference_name, {}).setdefault(mate, []).append(read)

    required = (1,) if mode == 'single' else (1, 2)
    if any(mate not in primary for mate in required):
        return set()
    originals = {mate: get_original_sequence(primary[mate]) for mate in required}
    if any(is_high_n_content(seq) for seq, qual in originals.values()):
        return set()

    def qualifying_clip(read, mate):
        length = len(originals[mate][0])
        if not read.cigartuples or length <= 2 * config.MIN_ANCHOR_LEN:
            return False
        # BWA may omit sequence or hard-clip alternative records. Use the
        # primary sequence length and count terminal hard clips as clipped bases.
        return any(op in (4, 5) and config.MIN_ANCHOR_LEN < size < length - config.MIN_ANCHOR_LEN
                   for op, size in (read.cigartuples[0], read.cigartuples[-1]))

    written = set()
    for is_name, alignments in by_is.items():
        clipped = any(qualifying_clip(read, mate)
                      for mate, records in alignments.items() for read in records)
        one_mapped = mode != 'single' and len(alignments) == 1
        if clipped or one_mapped:
            for mate in required:
                seq, qual = originals[mate]
                out_handles[is_name].write(format_fastq(primary[mate].query_name, seq, qual, str(mate)))
            written.add(is_name)
    return written


def run_extraction(is_db, read_files, output_paths, threads=16, mode='paired'):
    """Align once with alternative hits enabled and write per-IS FASTQs.

    output_paths maps original IS record IDs to output paths. Each query is
    emitted at most once per IS, with both mates together for paired input.
    Returns per-IS counts of extracted reads (single) or pairs (paired).
    """
    counts = dict.fromkeys(output_paths, 0)
    with ExitStack() as stack:
        handles = {name: stack.enter_context(open(path, 'w'))
                   for name, path in output_paths.items()}
        process = stack.enter_context(align_to_is_db(is_db, read_files, threads, mode))
        stream = stack.enter_context(pysam.AlignmentFile(process.stdout, 'r'))
        current_name = None
        buffer = []
        processed = 0

        def flush():
            for name in process_is_buffer(buffer, handles, mode):
                counts[name] += 1

        for read in stream:
            processed += 1
            if processed % 1000000 == 0:
                logging.info("Stage 1: Processed %s alignment records.", processed)
            if current_name is not None and read.query_name != current_name:
                flush()
                buffer = []
            current_name = read.query_name
            buffer.append(read)
        if buffer:
            flush()
    unit = 'reads' if mode == 'single' else 'read pairs'
    for name, count in counts.items():
        logging.info('Stage 1: Extracted %d candidate %s for %s.', count, unit, name)
    logging.info('Stage 1 complete: Processed %d alignment records.', processed)
    return counts
