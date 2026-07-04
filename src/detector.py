"""
Module: src/detector.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Stage 4 - Detect insertions by pairing peaks, identifying known sites, 
         annotating features, and detecting complex SVs (e.g., deletions).
"""
import pysam
import logging
from src import config
from collections import defaultdict


class InsertionResult:
    """
    Final output object representing a detected insertion or SV.
    """
    def __init__(self, chromosome, position, gap, orientation, category,
                 start_clipped, end_clipped, discordant_count, proper_count, is_length,
                 sv_type="None", group_id=None, peaks=None, is_singular=False):
        self.chromosome = chromosome
        self.position = position
        self.gap = gap # Reference length (Known) or Peak Distance (Novel).
        self.orientation = orientation
        self.category = category
        self.start_clipped = start_clipped
        self.end_clipped = end_clipped
        self.discordant_count = discordant_count
        self.proper_count = proper_count
        self.sv_type = sv_type # String describing associated SVs (e.g., "Deletion_100_200").
        self.is_length = is_length
        self.group_id = group_id
        self.peaks = peaks if peaks else [] # Tuple (start_peak, end_peak) or single Peak
        self.is_singular = is_singular  # True if this peak had no pairable counterpart

    def __repr__(self):
        return f"InsertionResult(Pos={self.position}, Gap={self.gap}, Orientation={self.orientation}, Category={self.category}, SV={self.sv_type}, IS_length={self.is_length}, Group={self.group_id}, Peaks={self.peaks}, Singular={self.is_singular})"
        
class PeakPair:
    """
    Helper object to evaluate potential pairs of peaks.
    """
    def __init__(self, start_peak, end_peak):
        self.start_peak = start_peak
        self.end_peak = end_peak
        self.gap = abs(start_peak.position - end_peak.position)
        self.group_id = None 
        if start_peak.is_coordinate is not None and end_peak.is_coordinate is not None:
            self.is_length = abs(start_peak.is_coordinate - end_peak.is_coordinate)
        else:
            self.is_length = "NA"
    def __repr__(self):
        return f"Pair(Start={self.start_peak}, End={self.end_peak}, IS_len={self.is_length})"        
            
def detect_insertions(cluster_dict, peak_list, is_len):
    """
    Main entry point for Stage 5. 
    Takes the raw outputs from Stage 4 (Cluster) and produces final InsertionResults.
    
    Args:
        cluster_dict: Dict {cluster_id: ClusterObj}
        peak_list: List[Peak]
        is_len: Integer length of the IS element
        
    Returns:
        List[InsertionResult]: The final, merged, sorted list of results.
    """
    
    groups = defaultdict(list)     # map: group_id -> list[Peak] including peaks related to deletions
    group_cores = defaultdict(list)    # map: group_id -> PeakPair (The core insertion object)
    
    assigned_peaks = set()
    
    
    # =========================================================================
    # PASS 1: Identify Core Pairs (Insertions)
    # =========================================================================
    pairs = []
    
    for i in range(len(peak_list)):
        p1 = peak_list[i]
        if p1 in assigned_peaks: continue
        
        # Look ahead for a mate
        for j in range(i+1, len(peak_list)):
            p2 = peak_list[j]
            if p2 in assigned_peaks: continue
            
            # Dynamic Gap Threshold
            limit = config.PAIR_GAP
            
            # If both are Reference/Restored from the SAME site (Lift-Over logic)
            # The gap between peaks will be roughly the IS length.
            if p1.known_ref and p2.known_ref and p1.known_ref == p2.known_ref:
                limit = is_len + 50
                
            # Break if too far (list is sorted)
            if abs(p2.position - p1.position) > limit:
                break     
            # Constraints: Orientation, flank_class
            if p1.orientation != p2.orientation: continue    
            if p1.flank_class == p2.flank_class: continue
            
            if p1.position != p2.position:
                # Standard Case: Upstream genomic position is Start
                p1, p2 = sorted([p1, p2], key=lambda x: x.position)
            else:
                # Edge Case: Same Genomic Position
                # Tie-breaker based on Viral Coordinate (is_coordinate)
                if p1.orientation == "F":
                    # Forward: Smaller viral coord is Start (5' end)
                    p1, p2 = sorted([p1, p2], key=lambda x: x.is_coordinate)
                elif p1.orientation == "R":
                    # Reverse: Larger viral coord is Start
                    p1, p2 = sorted([p1, p2], key=lambda x: x.is_coordinate, reverse=True)
            pair = PeakPair(p1, p2)
            # Constraints: Conflict reads

            reads_p1 = {(s['read_id'], s['is_read1']) for s in p1.soft_clipped}
            reads_p2 = {(s['read_id'], s['is_read1']) for s in p2.soft_clipped}
            if len(reads_p1.intersection(reads_p2)) > 0:
                continue
            pairs.append(pair)

    # Select Best Pairs (Greedy approach: smallest gap, then largest IS coverage)
    logging.info(f"[detect_insertions] Pass 1 candidate pairs: {len(pairs)} from {len(peak_list)} total peaks.")
    pairs.sort(key=lambda p: (p.gap, -p.is_length if isinstance(p.is_length, (int, float)) else 0))
            
    for pair in pairs:
        # Check usage again after sort
        if pair.start_peak in assigned_peaks or pair.end_peak in assigned_peaks: continue
        
        gid = f"{pair.start_peak.chromosome}_{pair.start_peak.position}_{pair.end_peak.position}"
        
        pair.group_id = gid
        group_cores[gid] = pair
        
        pair.start_peak.group_id = gid
        pair.end_peak.group_id = gid
        
        groups[gid].extend([pair.start_peak, pair.end_peak])
        
        assigned_peaks.add(pair.start_peak)
        assigned_peaks.add(pair.end_peak)

    logging.info(f"[detect_insertions] Pass 1 complete: {len(groups)} paired groups, "
                 f"{len(assigned_peaks)} peaks assigned out of {len(peak_list)} total peaks.")

    # =========================================================================
    # PASS 2: Aggregate Deletion Peaks (Physical Pairing)
    # =========================================================================
    # Collect only unassigned peaks with Deletion status.
    # Both peaks in a deletion pair must be unassigned — assigned_peaks is the
    # single source of truth, consistent with Pass 1.
    deletion_peaks = [p for p in peak_list if "Deletion" in p.sv_status]
    deletion_peaks.sort(key=lambda x: (x.chromosome, x.position))
    group_svs = defaultdict(list)

    for i in range(len(deletion_peaks)):
        d1 = deletion_peaks[i]
        if d1 in assigned_peaks: continue       # already grouped — skip
        if d1.sv_side != "Right": continue      # must be the Right cliff (start of deletion)

        for j in range(i + 1, len(deletion_peaks)):
            d2 = deletion_peaks[j]
            if d2 in assigned_peaks: continue   # already grouped — skip

            if d2.chromosome != d1.chromosome: break
            if (d2.position - d1.position) > config.MAX_SV_SIZE: break

            # d2 must be the Left cliff (end of deletion), same orientation as d1
            if d2.sv_side == "Left" and d1.orientation == d2.orientation:
                # Both peaks are unassigned: create a new deletion-only group.
                target_gid = f"{d1.chromosome}_{d1.position}_{d2.position}"
                d1.group_id = target_gid
                d2.group_id = target_gid
                group_cores[target_gid] = PeakPair(d1, d2)
                groups[target_gid].extend([d1, d2])
                assigned_peaks.add(d1)
                assigned_peaks.add(d2)

                group_svs[target_gid].append(f"Deletion_{d1.position}_{d2.position}")

                break  # move to next d1

    logging.info(f"[detect_insertions] Pass 2 complete: {len(group_svs)} deletion-only groups added, "
                 f"{len(groups)} total groups, {len(assigned_peaks)} peaks assigned.")

    # =========================================================================
    # PASS 3: Collect Unassigned (Singular) Peaks
    # =========================================================================
    # Any peak that was never paired in Pass 1 or Pass 2 is reported as-is,
    # flagged with is_singular=True so the writer can mark its position with "*".

    pass3_start = len(groups)
    for idx, p in enumerate(peak_list):
        if p in assigned_peaks:
            continue

        gid = f"{p.chromosome}_{p.cluster_id}_{p.position}_solo"
        p.group_id = gid
        groups[gid].append(p)
        assigned_peaks.add(p)

    logging.info(f"[detect_insertions] Pass 3 complete: {len(groups) - pass3_start} singular peaks added.")

    # =========================================================================
    # REPORT GENERATION
    # =========================================================================
    output_results = []

    for gid, members in groups.items():
        # No re-scanning for SVs needed here!
        members.sort(key=lambda x: (x.chromosome,x.position))
        
        core_pair = group_cores.get(gid)

        # 1. Retrieve Pre-calculated SV String
        sv_string = ",".join(group_svs[gid]) if group_svs[gid] else None

        pA = members[0] # Default representative peak

        # 2. Set defaults (orient must always be assigned before use)
        orient = pA.orientation
        final_gap = "NA"
        cat = "Novel"
        is_len = "NA"
        is_singular = False

        if core_pair:
            # --- Paired insertion ---
            pA = core_pair.start_peak
            pB = core_pair.end_peak
            orient = pA.orientation

            if pA.known_ref or pB.known_ref:
                cat = "Known(SV)" if ("Deletion" in pA.sv_status or "Deletion" in pB.sv_status) else "Known"
            else:
                cat = "Novel(SV)" if ("Deletion" in pA.sv_status or "Deletion" in pB.sv_status) else "Novel"

            if pA.is_coordinate is not None and pB.is_coordinate is not None:
                if 'SV' in cat:
                    final_gap = 'NA'
                else:
                    final_gap = core_pair.gap
            is_len = core_pair.is_length

        else:
            # --- Singular (unpaired) peak ---
            is_singular = True
            if pA.known_ref:
                cat = "Known(SV)" if "Deletion" in pA.sv_status else "Known"
            else:
                cat = "Novel(SV)" if "Deletion" in pA.sv_status else "Novel"

        # 3. Metadata
        cid = pA.cluster_id
        cluster = cluster_dict.get(cid)

        res = InsertionResult(
            chromosome=pA.chromosome,
            position=pA.position,
            gap=final_gap,
            orientation=orient,
            category=cat,
            start_clipped=sum(p.n_right_clips for p in members),
            end_clipped=sum(p.n_left_clips for p in members),
            discordant_count=cluster.discordant_count if cluster else 0,
            proper_count=cluster.proper_count if cluster else 0,
            is_length=is_len,
            sv_type=sv_string,
            group_id=gid,
            peaks=members,
            is_singular=is_singular
        )
        output_results.append(res)
    
    output_results.sort(key=lambda x: (x.chromosome, x.position))
    return output_results, groups, group_cores
