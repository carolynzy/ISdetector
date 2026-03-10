"""
Module: src/cluster.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Stage 3 - Group reads into clusters, identify insertion signals, 
         determine orientation, detect peaks, and map coordinates back 
         to the original reference genome.
"""
import statistics
import pysam
from Bio import Align
from Bio.Seq import Seq
from src import utils
from src import config
import numpy as np
from collections import Counter
from collections import defaultdict

class Cluster:
    """
    Data structure to store a group of reads and associated analysis data.
    """
    def __init__(self, chromosome, start_position, end_position, reads, discordant_count, proper_count):
        self.chromosome = chromosome
        self.start_position = start_position
        self.end_position = end_position
        self.reads = reads
        self.discordant_count = discordant_count
        self.proper_count = proper_count
    def __repr__(self):
        return f"Cluster(Chrom={self.chromosome}, Start={self.start_position}, End={self.end_position}, Reads={len(self.reads)}, Discordant={self.discordant_count}, Proper={self.proper_count})"

class InsertionSignal:
    """
    Intermediate object to store soft-clip data from a single read.
    """
    def __init__(self, chromosome, cluster_id, read_id, is_read1, position, clipped_sequence, is_coordinate, is_strand, orientation, clip_side):
        self.chromosome = chromosome
        self.cluster_id = cluster_id
        self.read_id = read_id
        self.is_read1 = is_read1
        self.position = position
        self.clipped_sequence = clipped_sequence
        self.is_coordinate = is_coordinate # 0-based index on IS
        self.is_strand = is_strand # IS Map Strand, 1 for Forward match, -1 for Reverse match
        self.orientation = orientation	
        self.clip_side = clip_side
    def __repr__(self):
        return f"InsertionSignal(Cluster={self.cluster_id}, Read={self.read_id}, Pos={self.position}, Clipped_seq={self.clipped_sequence}, IS_coordinate={self.is_coordinate}, IS_strand={self.is_strand}, Orientation={self.orientation}, Clip_side={self.clip_side})"

class Peak:
    """
    Object representing a merged consensus of insertion signals.
    """
    def __init__(self, chromosome, cluster_id, position, is_coordinate, orientation):
        self.chromosome = chromosome
        self.cluster_id = cluster_id
        self.position = position
        self.is_coordinate = is_coordinate
        self.orientation = orientation
        # SV Attributes
        self.sv_status = "Normal"  # Default
        self.sv_side = None        # "Left" or "Right"
        self.left_depth = 0.0
        self.right_depth = 0.0
        # Initialize support lists as empty (Best Practice)
        self.soft_clipped = [] # List of dictionary {'read_id': str, 'is_read1': bool}
        self.n_left_clips = 0  # Supports "End" of insertion
        self.n_right_clips = 0 # Supports "Start" of insertion
        self.flank_class = "NA"
        # Source Tracking
        self.source = "New"       # "New" (Novel Insertion) or "Reference" (Restored IS)
        self.known_ref = None     # Tag if it restores a known site
        self.group_id = None      # ID assigned during SV grouping/pairing steps
        
    def __repr__(self):
        return f"Peak(ID={self.cluster_id}, Pos={self.position}, IS_coordinate={self.is_coordinate}, Orientation={self.orientation}, Left_depth={self.left_depth}, Right_depth={self.right_depth}, Clipped={str(len(self.soft_clipped))}, Clipped_left={self.n_left_clips}, Clipped_right={self.n_right_clips}, Flank={self.flank_class}, Type={self.sv_status}, Side={self.sv_side}, Source={self.source}, Known_ref={self.known_ref}, Group_id={self.group_id})"
 
def filter_read(read):
    """
    Filter reads for peak detection
    Exclude MQ=0 and Secondary/Supplementary alignments.
    """
    if read.is_secondary or read.is_supplementary or read.is_unmapped:
        return False
    return True

def is_coordinate_strand(clipped_sequence, is_seq, clip_side):
    """
    Aligns clip to full IS. Returns (is_coordinate, is_strand).
    1. Score >= len(clip) * 1.5
    2. Aligned Length >= 20 bp
    3. Coverage >= 90% of clipped sequence
    """
    # Configure Aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    
    clip_len = len(clipped_sequence)
    
    def get_valid_score(aln):
        if aln is None: return -999
        # aln.aligned[0] refers to the first sequence passed to align() -> The Clip
        aligned_len = max(end - start for start, end in aln.aligned[0])
        coverage = aligned_len / clip_len
        if aligned_len < 20 or coverage < 0.90:
            return -999
        return aln.score
        
    # Check Forward
    try:
        aln_fwd = next(aligner.align(clipped_sequence, is_seq))
        score_fwd = get_valid_score(aln_fwd)
    except StopIteration:
        score_fwd = -999
        aln_fwd = None
    # Check Reverse Complement (of the clip)
    seq_rc = str(Seq(clipped_sequence).reverse_complement())
    try:
        aln_rev = next(aligner.align(seq_rc, is_seq))
        score_rev = get_valid_score(aln_rev)
    except StopIteration:
        score_rev = -999
        aln_rev = None
        
    # Calculate dynamic score threshold
    min_score = len(clipped_sequence) * 1.5
    if max(score_fwd, score_rev) < min_score:
        return None, 0
        
    # Determine best orientation
    if score_fwd >= score_rev:
        # Forward Match
        is_strand = 1
        target_intervals = aln_fwd.aligned[1]
        if clip_side == "RIGHT":
            # [Genome] -> [Clip (IS Start)] -> Junction at START of match
            is_coord = target_intervals[0][0]
        else: # LEFT
            # [Clip (IS End)] -> [Genome] -> Junction at END of match
            is_coord = target_intervals[0][1]
    else:
        is_strand = -1
        target_intervals = aln_rev.aligned[1]
        # Reverse Match (RC aligned)
        if clip_side == "RIGHT":
            # Orig: [Genome] -> [Clip (Start)]
            # RC:   [Clip_RC (End)] <- [Genome]
            is_coord = target_intervals[0][1]
        else: # LEFT
            # Orig: [Clip (End)] -> [Genome]
            # RC:   [Genome] <- [Clip_RC (Start)]
            is_coord = target_intervals[0][0]
            
    return is_coord, is_strand

def extract_signals_from_cluster(reads, cluster_id, chromosome, is_seq):
    """
    Generate InsertionSignal objects from reads in a cluster.
    """
    signals = []
    for read in reads:
        if not filter_read(read):
            continue
        cigartuples = read.cigartuples
        if not cigartuples:
            continue
        
        seq_len = read.query_length
        if seq_len <= (config.MIN_ANCHOR_LEN * 2):
            continue
            
        min_clip = config.MIN_ANCHOR_LEN
        max_clip = seq_len - config.MIN_ANCHOR_LEN
            
        # If paired, use actual flag. If unpaired (single), treat as 'Read 1' (Primary).
        is_read1_flag = read.is_read1 if read.is_paired else True
        
        # Left Clip (Index 0)
        op_start, len_start = cigartuples[0]
        if op_start == 4 and min_clip < len_start < max_clip:
            pos = read.reference_start
            clipped_seq = read.query_sequence[:len_start]
            clip_side = "LEFT"

            is_coord, is_strand = is_coordinate_strand(clipped_seq, is_seq,  clip_side)
            if is_strand != 0:
                orient = "F" if is_strand == 1 else "R"
                sig = InsertionSignal(chromosome, cluster_id, read.query_name, is_read1_flag, pos, clipped_seq, is_coord, is_strand, orient, clip_side)
                signals.append(sig)
        # Right Clip (Index -1)
        op_end, len_end = cigartuples[-1]
        if op_end == 4 and min_clip < len_end < max_clip:
            pos = read.reference_end
            clipped_seq = read.query_sequence[-len_end:]
            clip_side = "RIGHT"
            is_coord, is_strand = is_coordinate_strand(clipped_seq, is_seq, clip_side)
            if is_strand != 0:
                orient = "F" if is_strand == 1 else "R"
                sig = InsertionSignal(chromosome, cluster_id, read.query_name, is_read1_flag, pos, clipped_seq, is_coord, is_strand, orient, clip_side)
                signals.append(sig)
    return signals

def get_mode_or_median(values):
    if not values: 
        return 0        
    c = Counter(values)
    most_common = c.most_common() 
    top_val, top_count = most_common[0]    
    
    use_median = False    
    if top_count == 1:
        use_median = True       
    elif len(most_common) > 1 and most_common[1][1] == top_count:
        use_median = True        
    
    if use_median:
        return int(np.median(values))        
    return top_val

def _get_flanking_depths(bam_stream, chrom, pos, window=20):
    """
    Calculates average read depth in immediate left and right flanking regions.
    """
    l_start = max(0, pos - window)
    l_end = pos - 1       
    r_start = pos + 1
    r_end = pos + window
    
    if l_end > l_start:
        l_stats = bam_stream.count_coverage(chrom, start=l_start, stop=l_end, quality_threshold=1)
        l_depths = np.array(l_stats).sum(axis=0)
        l_mean = np.mean(l_depths) if l_depths.size > 0 else 0.0
    else:
        l_mean = 0.0
        
    r_stats = bam_stream.count_coverage(chrom, start=r_start, stop=r_end, quality_threshold=1)
    r_depths = np.array(r_stats).sum(axis=0)
    r_mean = np.mean(r_depths) if r_depths.size > 0 else 0.0
    
    return l_mean, r_mean
    
def _derive_sv_status(peak, bam_stream, window=20):
    """
    Calculates depth and assigns sv_status/sv_side.
    """
    l_d, r_d = _get_flanking_depths(bam_stream, peak.chromosome, peak.position, window)
    peak.left_depth = l_d
    peak.right_depth = r_d
    
    # Cliff Logic
    if l_d < config.ZERO_THRESH and r_d < config.ZERO_THRESH:
        peak.sv_status = "Low_Coverage"; peak.sv_side = None
    elif r_d < config.ZERO_THRESH:
        peak.sv_status = "Deletion_coexistence"; peak.sv_side = "Right"
    elif l_d < config.ZERO_THRESH:
        peak.sv_status = "Deletion_coexistence"; peak.sv_side = "Left"
    else:
        ratio = l_d / r_d if r_d > 0 else 999
        if ratio < config.DEPTH_DIFF_THRESHOLD:
            peak.sv_status = "Deletion_coexistence"; peak.sv_side = "Left"
        elif ratio > (1.0 / config.DEPTH_DIFF_THRESHOLD):
            peak.sv_status = "Deletion_coexistence"; peak.sv_side = "Right"
        else:
            peak.sv_status = "Normal"; peak.sv_side = None
    
def create_peak_from_signals(signal_list, bam_stream):
    """
    Helper to aggregate signal list into a Peak object.
    Removed known_sites logic; classification happens during lift-over.
    """
    first = signal_list[0]
    # 1. Genomic Position
    positions = [s.position for s in signal_list]
    peak_pos = get_mode_or_median(positions)

    # 2. IS Coordinate
    is_coords = [s.is_coordinate for s in signal_list]
    peak_is_coord = get_mode_or_median(is_coords)

    # 3. Orientation 
    peak_orient = first.orientation
    peak = Peak(first.chromosome, first.cluster_id, peak_pos, peak_is_coord, peak_orient)
    
    # Aggregate clip counts
    for s in signal_list:
        peak.soft_clipped.append({
            'read_id': s.read_id, 
            'is_read1': s.is_read1,
            'clip_side': s.clip_side 
        })
        if s.clip_side == "LEFT":
            peak.n_left_clips += 1
        elif s.clip_side == "RIGHT":
            peak.n_right_clips += 1
            
    # Assign Category
    if peak.n_left_clips > config.ZERO_THRESH:
        peak.flank_class = "End" 
    elif peak.n_right_clips > config.ZERO_THRESH:
        peak.flank_class = "Start" 

    
    _derive_sv_status(peak, bam_stream)

    return peak
    
def detect_all_peaks(reads, cluster_id, chromosome, is_seq, bam_stream):
    """
    Extracts signals and merges them into Peaks using a Centroid-based approach (Median).
    Returns signals and detected peaks.
    """
    # 1. Extract
    signals = extract_signals_from_cluster(reads, cluster_id, chromosome, is_seq)
    if not signals: return [], []
    
    # 2. Filter
    valid_signals = [s for s in signals if s.is_coordinate is not None]
    if not valid_signals:
        return [], []
        
    # 3. Bucket signals by (Orientation, Clip_dise)
    buckets = defaultdict(list)
    
    for sig in valid_signals:
        # Key example: ('F', 'LEFT') or ('R', 'RIGHT')
        key = (sig.orientation, sig.clip_side)
        buckets[key].append(sig)
        
    # 3. Merge Strategy
    peaks = []

    for (orient, clipside), sig_list in buckets.items():
        if not sig_list: continue
        sig_list.sort(key=lambda x: (x.chromosome, x.position, x.is_coordinate))
        
        # This list holds the current clusters being built
        # Format: [{'signals': [sig1, sig2...], 'median_pos': 100, 'median_is': 500}, ...]
        active_groups = []
        
        for sig in sig_list:
            # Get values for the candidate signal
            sig_pos = sig.position
            sig_is = sig.is_coordinate
            
            best_group_idx = -1
            best_dist_score = float('inf') # Closer is better
            
            for i, group_data in enumerate(active_groups):# Calculate Reference values (Median) of the CURRENT group
                ref_pos = group_data['median_pos']
                ref_is = group_data['median_is']
            
                # 4. COMPARISON: Candidate vs Group Median
                # Check Genomic Position Distance
                pos_diff = abs(sig_pos - ref_pos)
            
                # Check IS Coordinate Distance (Only if both have valid coordinates)
                is_diff = abs(sig_is - ref_is)
            
                if pos_diff <= config.PEAK_DISTANCE and is_diff <= config.PEAK_DISTANCE:
                    # It's a candidate! Is it the *best* candidate?
                    total_score = pos_diff + is_diff
                    if total_score < best_dist_score:
                        best_dist_score = total_score
                        best_group_idx = i
            
            if best_group_idx != -1:
                # Add to existing group
                target = active_groups[best_group_idx]
                target['signals'].append(sig)
                
                # Update Medians (Recalculate to keep centroid accurate)
                all_pos = [s.position for s in target['signals']]
                all_is = [s.is_coordinate for s in target['signals']]
                
                target['median_pos'] = statistics.median(all_pos)
                target['median_is'] = statistics.median(all_is)
            else:
                # Create NEW peak
                new_group = {
                    'signals': [sig],
                    'median_pos': sig_pos,
                    'median_is': sig_is
                }
                active_groups.append(new_group)
    
        for group_data in active_groups:
            if len(group_data['signals']) > config.MIN_PEAK_SUPPORT:
                peaks.append(create_peak_from_signals(group_data['signals'], bam_stream))
    return signals, peaks

def finalize_cluster(reads, current_chrom):
    """
    Creates a Cluster object if minimum support requirements are met.
    Returns: cluster_id, cluster_obj
    """
    if not reads:
        return None, None
        
    ref_starts = [r.reference_start for r in reads]
    ref_ends = [r.reference_end for r in reads if r.reference_end is not None]
    
    start_position = min(ref_starts)
    end_position = max(ref_ends) if ref_ends else max(ref_starts)
    cluster_id = f"{start_position}_{end_position}"

    discordant_count = 0
    proper_count = 0
    for r in reads:
        if r.is_paired and r.is_proper_pair:
            proper_count += 1
        elif r.is_paired and not r.is_proper_pair:
            discordant_count += 1
        else:
            # For Single-End, everything is technically "discordant" in the sense of 
            # not being a proper pair, or we can count them as a separate category.
            # Here we count them as discordant to ensure they are tracked as 'support'.
            discordant_count += 1

    if len(reads) < config.MIN_CLUSTER_SUPPORT:
        return None, None

    cluster_obj = Cluster(current_chrom, start_position, end_position, reads, discordant_count, proper_count)

    return cluster_id, cluster_obj

def create_microclusters(input_source):
    """
    Main entry point for Stage 3.
    Iterates through sorted reads, groups them by distance, and creates Cluster objects.
    """
    clusters = {}

    bam_file = None
    if isinstance(input_source, str):
        try:
            bam_file = pysam.AlignmentFile(input_source, "rb")
            read_stream = bam_file.fetch()
        except ValueError:
            bam_file = pysam.AlignmentFile(input_source, "r")
            read_stream = bam_file.fetch()
    else:
        read_stream = input_source

    current_reads = []
    prev_pos = None
    current_chrom = None
    
    try:
        for read in read_stream:
            if not filter_read(read):
                continue
                
            pos = read.reference_start
            chrom = read.reference_name
            
            # Check if we need to close the previous cluster
            if current_chrom is not None and prev_pos is not None:
                if (chrom != current_chrom or (pos - prev_pos) >= config.CLUSTER_DIST):
                    if len(current_reads) >= config.MIN_CLUSTER_SUPPORT: 
                        c_id, c_obj = finalize_cluster(current_reads, current_chrom)
                        if c_id:
                            clusters[c_id] = c_obj
                  
                    current_reads = []
            
            current_reads.append(read)
            prev_pos = read.reference_end
            current_chrom = chrom
            
        # Process the last cluster
        if current_reads and len(current_reads) >= config.MIN_CLUSTER_SUPPORT:
            c_id, c_obj = finalize_cluster(current_reads, current_chrom)
            if c_id: 
                clusters[c_id] = c_obj

    finally:
        if bam_file:
            bam_file.close()

    return clusters

def lift_over_peak(peak, crossmap):
    """
    Translates 'clean' coordinates to 'original' coordinates using the CrossMap.
    Classifies peaks as 'Reference' (Restored) or 'New' based on stitch point overlap.
    """
    clean_id = peak.chromosome
    if clean_id not in crossmap:
        return # No changes needed

    # Assume original name is clean_id without "_clean" suffix
    # or rely on upstream tracking. Here we strip suffix for display.
    original_id = clean_id.replace("_clean", "")
    
    deletions = crossmap[clean_id]
    cumulative_shift = 0
    is_restored = False
    clean_pos = peak.position
    
    # Deletions in crossmap are ordered by clean_pos (as they were built)
    for record in deletions:
        # Check for overlap with the stitch point (Junction of removal)
        # Tolerance of +/- 5bp for alignment jitter
        if abs(clean_pos - record['new_pos']) <= 5:
            is_restored = True
            
        # If the peak is AFTER this deletion, add the deleted length back
        if clean_pos > record['new_pos']:
            cumulative_shift += record['deleted_len']
            
    # Apply translation
    peak.position += cumulative_shift
    peak.chromosome = original_id
    
    if is_restored:
        peak.source = "Reference"
        peak.known_ref = "Restored_IS"
    else:
        peak.source = "New"

def scan_clusters_for_peaks(cluster_dict, is_seq, bam_path, crossmap=None):
    """
    Stage 3 (Part B): Signal Extraction and Peak Detection.
    1. Detect Peaks (on IS-Clean reference).
    2. Lift-Over coordinates to Original Reference.
    3. Classify (Known/Restored vs. Novel).
    """
    all_signals = []
    all_peaks = []
    
    with pysam.AlignmentFile(bam_path, "rb") as bam_stream:
        # 1. Detect Observed Peaks
        for c_id, cluster_obj in cluster_dict.items():
            sigs, pks = detect_all_peaks(
                cluster_obj.reads, 
                c_id, 
                cluster_obj.chromosome, 
                is_seq, 
                bam_stream
            )
            all_signals.extend(sigs)
            all_peaks.extend(pks)
            
    # 2. Lift-Over and Classify
    if crossmap:
        for peak in all_peaks:
            lift_over_peak(peak, crossmap)
    
    # 3. Sort (Now by original coordinates)
    all_peaks.sort(key=lambda x: (x.chromosome, x.position))
        
    return all_signals, all_peaks
