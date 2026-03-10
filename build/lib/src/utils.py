"""
Module: src/utils.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Helper functions for I/O, filtering, and data management.
         Includes handling of FASTA indexing, GenBank parsing, and Reporting.
"""
import os
import sys
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pandas as pd
from Bio import SeqFeature

class Genebank:
    """
    Parser and manager for GenBank reference files.
    """
    def __init__(self, filepath):
        """
        Initialize and load GenBank records into memory.
        """
        self.filepath = filepath
        self.records = {} # Dictionary { "chromosome": SeqRecord }
        
        if not os.path.exists(filepath):
            logging.error(f"GenBank file not found: {filepath}")
            sys.exit(1)
            
        try:
            logging.info(f"Loading GenBank file: {filepath}")
            self.records = SeqIO.to_dict(SeqIO.parse(filepath, "genbank"))
        except Exception as e:
            logging.error(f"Failed to parse GenBank file: {e}")
            sys.exit(1)

    def convert_to_fasta(self, output_path):
        """
        Extract DNA sequences and write to a FASTA file.
        """
        try:
            if output_path is None:
                base, _ = os.path.splitext(self.filepath)
                output_path = f"{base}.fasta"
            
            logging.info(f"Converting GenBank to FASTA: {output_path}")
            SeqIO.write(self.records.values(), output_path, "fasta")
            return output_path
        except Exception as e:
            logging.error(f"Failed to convert to FASTA: {e}")
            sys.exit(1)

    def extract_tag_info(self, chromosome, position):
        """
        Finds the genomic feature at the given position.
        Prioritizes: CDS > tRNA/rRNA > gene > mobile_element > repeat_region.
        """
        record = self.records.get(chromosome)
        
        # Handle cases where record ID might differ slightly or file structure varies
        if record is None:
            # Fallback check if the dictionary keys are not matching directly
            for r_id, r_obj in self.records.items():
                if r_id == chromosome or r_obj.id == chromosome:
                    record = r_obj
                    break
        
        if record is None:
            return "NA", "NA"

        overlapping_features = []
        for feature in record.features:
            if feature.type == "source":
                continue
            if feature.location.start <= position <= feature.location.end:
                overlapping_features.append(feature)

        if not overlapping_features:
            return "intergenic", "NA"

        priority = {
            "CDS": 1, "tRNA": 2, "rRNA": 2, "gene": 3, "mobile_element": 4, "repeat_region": 5
        }
        best_feature = None
        best_score = 100

        for feat in overlapping_features:
            ftype = feat.type
            score = priority.get(ftype, 10)
            if score < best_score:
                best_score = score
                best_feature = feat

        if best_feature is None:
            return "intergenic", "NA"

        f_type = best_feature.type
        quals = best_feature.qualifiers
        name = "NA"

        if f_type == "CDS":
            name = quals.get("gene", quals.get("locus_tag", ["NA"]))[0]
        elif f_type in ["tRNA", "rRNA"]:
            name = quals.get("gene", quals.get("product", ["NA"]))[0]
        elif f_type == "gene":
            name = quals.get("gene", quals.get("locus_tag", ["NA"]))[0]
        elif f_type == "mobile_element":
            name = quals.get("mobile_element_type", ["NA"])[0]
        else:
            name = quals.get("locus_tag", ["NA"])[0]

        return f_type, name

def generate_clean_reference(ref_fasta, is_seq, is_len, output_dir, identity_cutoff=95.0, coverage_cutoff=0.90, suffix=""):
    """
    Maps IS to reference, removes it, and records coordinate mapping.
    
    Args:
        ref_fasta (str): Path to original reference FASTA.
        is_seq (str): The specific IS sequence to mask.
        is_len (int): Length of the IS.
        output_dir (str): Output directory.
        identity_cutoff (float): Minimum BLAST identity to detect IS.
        coverage_cutoff (float): Minimum BLAST coverage to detect IS.
        suffix (str): Optional suffix for the output filename (e.g., "_ISSd1").
        
    Returns:
        tuple: (path_to_clean_fasta, cross_mapping_dict)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Create a temporary IS Query File (Unique per suffix to avoid collisions)
    temp_is_path = os.path.join(output_dir, f"temp_is_query{suffix}.fasta")
    with open(temp_is_path, "w") as f:
        f.write(f">IS_Query\n{str(is_seq)}\n")

    print(f"[INFO] Running BLAST against reference genome (Suffix: {suffix})...")

    # 2. Run BLASTn ONCE against the entire reference file
    # We ask for 'sseqid' (Subject Sequence ID) to know which contig the hit belongs to.
    cmd = [
        "blastn", 
        "-query", temp_is_path, 
        "-subject", ref_fasta,  # <--- Pass the whole file here!
        "-outfmt", "6 sseqid sstart send length pident qlen", 
        "-task", "blastn" 
    ]
        
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate()
        
    # 3. Parse BLAST results and Group by Contig ID
    # Structure: hits_by_contig = { "chrom1": [(start, end), ...], "chrom2": [...] }
    hits_by_contig = defaultdict(list)
    
    if stdout:
        lines = stdout.strip().split('\n')
        for line in lines:
            if not line: continue
            cols = line.split('\t')
            contig_id = cols[0]        # sseqid
            s_start, s_end = int(cols[1]), int(cols[2])
            match_len = int(cols[3])
            pident = float(cols[4])
            
            # BLAST is 1-based, convert to 0-based
            start_0 = min(s_start, s_end) - 1
            end_0 = max(s_start, s_end)
            
            cov = match_len / is_len
            if pident >= identity_cutoff and cov >= coverage_cutoff:
                hits_by_contig[contig_id].append((start_0, end_0))

    # 4. Process Every Contig in the Reference
    ref_records = list(SeqIO.parse(ref_fasta, "fasta"))
    clean_records = []
    full_crossmap = {} 

    print(f"[INFO] Cleaning {len(ref_records)} contigs...")

    for record in ref_records:
        contig_id = record.id
        contig_seq = str(record.seq)
        
        # Get hits specific to THIS contig (if any)
        hits = hits_by_contig.get(contig_id, [])
        
        # --- Merge Overlapping Hits ---
        merged_hits = []
        if hits:
            hits.sort()
            curr_start, curr_end = hits[0]
            for next_start, next_end in hits[1:]:
                if next_start < curr_end:
                    curr_end = max(curr_end, next_end)
                else:
                    merged_hits.append((curr_start, curr_end))
                    curr_start, curr_end = next_start, next_end
            merged_hits.append((curr_start, curr_end))

        # --- Rebuild Sequence (Remove IS regions) ---
        new_seq_parts = []
        current_ref_pos = 0
        contig_crossmap = []

        for (del_start, del_end) in merged_hits:
            # Append sequence BEFORE the IS
            chunk = contig_seq[current_ref_pos:del_start]
            new_seq_parts.append(chunk)
            
            # Record the coordinate shift
            new_junction_pos = len("".join(new_seq_parts)) 
            deleted_len = del_end - del_start
            
            contig_crossmap.append({
                'new_pos': new_junction_pos,
                'deleted_len': deleted_len,
                'orig_start': del_start,
                'orig_end': del_end
            })
            
            current_ref_pos = del_end

        # Append remaining sequence after the last IS
        new_seq_parts.append(contig_seq[current_ref_pos:])
        final_clean_seq = "".join(new_seq_parts)

        # Create new SeqRecord
        new_id = f"{contig_id}_clean"
        clean_record = SeqRecord(
            Seq(final_clean_seq),
            id=new_id,
            description=f"IS-clean version of {contig_id}"
        )
        clean_records.append(clean_record)
        full_crossmap[new_id] = contig_crossmap

        if len(merged_hits) > 0:
            print(f"  -> {contig_id}: Removed {len(merged_hits)} IS regions.")

    # 5. Save Output
    output_fasta = os.path.join(output_dir, f"reference_clean{suffix}.fasta")
    SeqIO.write(clean_records, output_fasta, "fasta")
    
    if os.path.exists(temp_is_path):
        os.remove(temp_is_path)
    
    return output_fasta, full_crossmap

def indexing(fasta):
    """
    Check if indices exist. If not, run bwa index.
    """
    if not os.path.exists(f"{fasta}.bwt"):
        logging.info(f"Indexing {fasta}...")
        subprocess.check_call(["bwa", "index", fasta])

def bam_to_fastq(read):
    """
    Format the read object into a standard FASTQ string.
    """
    seq = read.query_sequence
    qual = read.qual
    return f"@{read.query_name}\n{seq}\n+\n{qual}\n"

def write_insertion_report(results, output_path):
    """
    Write the final report of insertions sample_insertion_reports.tsv.
    """
    header = [
        "Chromosome", "Position", "Gap", "Orientation", "IS_Length", "Category", 
        "Start_clipped", "End_clipped", "Discordant", "Proper", "SV"
    ]
    with open(output_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for res in results:
            line = [
                str(res.chromosome), 
                str(res.position), 
                str(res.gap), 
                str(res.orientation), 
                str(res.is_length), 
                str(res.category),
                str(res.start_clipped), 
                str(res.end_clipped), 
                str(res.discordant_count), 
                str(res.proper_count),
                str(res.sv_type)
            ]
            f.write("\t".join(line) + "\n")

def write_signal_report(signal_list, output_path):
    """
    Write the debug file sample_signal_report.tsv.
    """
    header = ["Chromosome", "Cluster_id", "Read_id", "is_read1", "Position", "IS_coordinate", "Clipped_sequence", "IS_strand", "Orientation", "Clip_side"]
    
    with open(output_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for sig in signal_list:
            line = [
                sig.chromosome,
                str(sig.cluster_id),
                sig.read_id,
                str(sig.is_read1),
                str(sig.position),
                str(sig.is_coordinate),
                sig.clipped_sequence,
                str(sig.is_strand),
                sig.orientation,
                sig.clip_side
            ]
            f.write("\t".join(line) + "\n")
            
def get_gene_name(feature):
    """
    Helper to extract gene name. 
    Priority: 'gene' > 'locus_tag' > 'Unknown'
    """
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    return "Unknown"            
            
def generate_annotation_report(insertion_list, gb_obj):
    """
    Annotates final insertion results.
    - If SV (Deletion): Reports genes overlapping the deleted region.
    - If Insertion: Reports gene at site or flanking genes.
    - 'Gene' and 'Protein' columns both use the gene name/locus_tag.
    """
    annotations = []
    
    # Access BioPython records from the Genebank object
    # Assuming gb_obj.records is a dict {chrom: SeqRecord} or list
    if isinstance(gb_obj.records, list):
        records_map = {r.id: r for r in gb_obj.records}
    else:
        records_map = gb_obj.records

    for ins in insertion_list:
        chrom = ins.chromosome
        pos = ins.position
        
        if chrom not in records_map:
            continue
            
        record = records_map[chrom]
        
        # --- Logic Branch: SV vs Point Insertion ---
        is_sv = False
        sv_start, sv_end = 0, 0
        
        # Check if it's a Deletion SV
        if ins.sv_type and "Deletion" in ins.sv_type:
            try:
                # Expected format: Deletion_Start_End
                parts = ins.sv_type.split('_')
                if len(parts) >= 3:
                    sv_start = int(parts[1])
                    sv_end = int(parts[2])
                    is_sv = True
            except ValueError:
                pass

        # 1. SV ANNOTATION (Range Overlap)
        if is_sv:
            affected_genes = []
            affected_products = []
            
            for feat in record.features:
                if feat.type != "CDS": continue
                
                # Check for Overlap: Feature Start < SV End AND Feature End > SV Start
                f_start = int(feat.location.start)
                f_end = int(feat.location.end)
                
                if f_start < sv_end and f_end > sv_start:
                    g_name = get_gene_name(feat)
                    prod = feat.qualifiers.get('product', ['NA'])[0]
                    affected_genes.append(g_name)
                    affected_products.append(prod)
            
            if affected_genes:
                gene_str = ";".join(affected_genes)
                prod_str = ";".join(affected_products)
                loc_type = "Exonic (Deletion)"
            else:
                gene_str = "Intergenic"
                prod_str = "NA"
                loc_type = "Intergenic (Deletion)"
                
            entry = {
                "Group_ID": ins.group_id,
                "Position": pos,
                "SV_Type": ins.sv_type,
                "Annotation_Type": loc_type,
                "Gene": gene_str,
                "Protein": gene_str,  # Requested: same as Gene
                "Product": prod_str
            }
            annotations.append(entry)

        # 2. POINT INSERTION ANNOTATION (Site/Nearest)
        else:
            found_hit = False
            # Check for direct insertion into a gene
            for feat in record.features:
                if feat.type != "CDS": continue
                if int(feat.location.start) <= pos <= int(feat.location.end):
                    g_name = get_gene_name(feat)
                    prod = feat.qualifiers.get('product', ['NA'])[0]
                    
                    entry = {
                        "Group_ID": ins.group_id,
                        "Position": pos,
                        "SV_Type": "None",
                        "Annotation_Type": "Intragenic",
                        "Gene": g_name,
                        "Protein": g_name, # Requested: same as Gene
                        "Product": prod
                    }
                    annotations.append(entry)
                    found_hit = True
                    break # Usually one insertion is in one gene (unless overlapping genes)
            
            # If Intergenic, find nearest Upstream/Downstream
            if not found_hit:
                upstream = None
                downstream = None
                min_up_dist = float('inf')
                min_down_dist = float('inf')
                
                for feat in record.features:
                    if feat.type != "CDS": continue
                    f_start = int(feat.location.start)
                    f_end = int(feat.location.end)
                    
                    # Gene is Upstream (Gene End < Pos)
                    if f_end < pos:
                        dist = pos - f_end
                        if dist < min_up_dist:
                            min_up_dist = dist
                            upstream = feat
                    
                    # Gene is Downstream (Gene Start > Pos)
                    elif f_start > pos:
                        dist = f_start - pos
                        if dist < min_down_dist:
                            min_down_dist = dist
                            downstream = feat
                
                # Format Result
                up_str = f"{get_gene_name(upstream)} (-{min_up_dist})" if upstream else "None"
                down_str = f"{get_gene_name(downstream)} (+{min_down_dist})" if downstream else "None"
                
                entry = {
                    "Group_ID": ins.group_id,
                    "Position": pos,
                    "SV_Type": "None",
                    "Annotation_Type": "Intergenic",
                    "Gene": f"{up_str} / {down_str}",
                    "Protein": f"{up_str} / {down_str}",
                    "Product": "NA"
                }
                annotations.append(entry)

    return annotations
