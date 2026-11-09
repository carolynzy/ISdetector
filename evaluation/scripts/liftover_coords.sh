#!/usr/bin/env bash
# =============================================================================
# liftover_coords.sh
#
# Build a reciprocal-best chain from an existing PAF alignment and optionally
# lift BED coordinates through it.
#
# Pipeline:
#   1. Compute chromosome sizes for sample and reference FASTAs
#   2. Convert the provided PAF to chain (inline Python)
#   3. Sort chain (chainSort)
#   4. Build reciprocal-best netted chain
#      (chainPreNet, chainNet, netSyntenic, netChainSubset,
#       chainStitchId, chainSwap, chainSort; no netFilter -syn)
#   5. If a BED was provided, perform chain-based direct liftover:
#        - Parse liftover chain once; key by source (sample) chromosome
#        - For each record, select target chromosome by source-bp overlap,
#          then chain score as tiebreaker
#        - Map each boundary via exact offset (if inside an aligned block)
#          or nearest block edge (if in a gap between blocks)
#        - Status flags:
#            exact              — both boundaries fall inside aligned blocks;
#                                 coordinate is a direct offset mapping
#            gap_snapped        — at least one boundary fell in an unaligned
#                                 gap; nearest block edge used; gap distance
#                                 reported
#            insertion_site     — record carries an encoded IS body; mapped
#                                 BED boundaries reported as the site
#            inversion_boundary — IS flanks map to opposite chain strands
#            unmapped           — no chain blocks found on this source chromosome
#
# Usage:
#   bash liftover_coords.sh -p <input_paf> -s <sample_fasta> \
#                           -r <reference_fasta> [-b <input_bed>] \
#                           [--chain [name]] [-o <output_dir>]
#
# PAF orientation:
#   PAF query names/coordinates must describe the sample/source assembly, and
#   PAF target names/coordinates must describe the reference/target assembly.
#
# Required tools:
#   faSize, chainSort, chainPreNet, chainNet, netSyntenic,
#   netChainSubset, chainStitchId, chainSwap, python3
#
# Example:
#   bash liftover_coords.sh \
#       -p sample_to_reference.paf \
#       -s sample.fasta \
#       -r reference.fasta \
#       -b sample_true.bed \
#       -o results/
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
OUTPUT_DIR="."
CHAINNET_MIN_SCORE=200
CHAINNET_MIN_SPACE=10
CHAINNET_MIN_FILL=10
SV_FLANK_TOL=10
PAF_TAG_MODE="all"
CHAIN_OUTPUT_NAME=""
INPUT_CHAIN=""

# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------
usage() {
    cat <<'EOF'
Usage:
  liftover_coords.sh -p PAF -s SAMPLE_FASTA -r REFERENCE_FASTA [options]
  liftover_coords.sh --input-chain CHAIN -b BED [options]

Build a reciprocal-best chain from an existing PAF alignment and optionally
lift BED intervals through it. Alternatively, --input-chain bypasses chain
generation and lifts BED intervals directly through a supplied chain.

Required for PAF/chain-generation mode:
  -p PAF              Sample-to-reference alignment in PAF format. Records
                      must contain a minimap2 cg:Z: CIGAR tag.
  -s SAMPLE_FASTA     FASTA for the source/sample assembly named by the BED.
  -r REFERENCE_FASTA  FASTA for the target/reference assembly.

Options:
  -b BED              BED file whose sample coordinates will be lifted. Omit
                      this option to generate only the reciprocal-best chain.
  -o OUTPUT_DIR       Directory for results and intermediate files
                      (default: current directory).
  --input-chain FILE  Lift BED records through an existing chain. This skips
                      PAF conversion and reciprocal-best chain generation;
                      -s, --chain, and chainNet options are ignored. Optional
                      -p supplies PAF CIGAR evidence for flanking SV details
                      (filtered by --tag).
                      The chain must have sample/source as t and target as q.
                      In this mode -b is required; -r is optional and is used
                      only to name BED/report outputs (default label: target).
  --chain [NAME]      Save the final reciprocal-best chain as NAME in the
                      output directory. The .chain suffix is added when
                      absent. If NAME is omitted, use <PAF-basename>.chain.
  --min-score N       Minimum chain score considered by chainNet
                      (default: 200).
  --min-space N       Minimum unoccupied gap size available for filling
                      (default: 10).
  --min-fill N        Minimum aligned fill size retained by chainNet
                      (default: 10).
  --tag MODE          Select PAF alignments admitted to reciprocal-best
                      chaining: "all" includes every alignment with a cg:Z:
                      CIGAR; "primary" includes only records tagged tp:A:P
                      (default: all).
  --sv-flank-tol N    Report gap SVs (deletions/insertions) lying within N bp
                      of the IS on both the source and the lifted reference
                      position, on either flank, and reverse alignments within
                      N bp of either flank as inversions (default: 10). Only
                      sv_type is affected; lifted positions are not.
  -h, --help          Show this help message and exit.

PAF orientation:
  PAF query names and coordinates must describe SAMPLE_FASTA; PAF target names
  and coordinates must describe REFERENCE_FASTA. In the generated chain, the
  sample is the source (t) assembly and the reference is the target (q).

Alignment filtering:
  Minimap2 PAF uses tp:A:P for primary records and other tp values for
  non-primary records. PAF does not carry the SAM supplementary-alignment flag,
  so primary split segments tagged tp:A:P remain in "primary" mode. Records
  without a tp tag are included in "all" mode and excluded in "primary" mode.

Main outputs:
  <bed>_to_<reference>.bed          Successfully lifted BED records.
  <bed>_to_<reference>.unmapped     Records with no usable chain mapping
                                    (written only when any are unmapped).
  <bed>_to_<reference>.lift_report.tsv
                                    Per-record mapping status, gap distances,
                                    insertion context, and structural-variant
                                    classification.
  <chain-name>.chain                Generated reciprocal-best chain (PAF mode).
  liftover_<date>_<paf>.log         Full execution log.

Mapping statuses:
  exact               Both BED boundaries map inside aligned chain blocks.
  gap_snapped         One or both boundaries are snapped to a nearby block.
  insertion_site      An encoded IS body with mapped BED boundary positions.
  inversion_boundary  Opposite-strand flanks indicate an inversion boundary.

Each input interval has one lift-report row. original_bed_col6 preserves input
BED column 6 (or '.' when absent). The final column, sv_type, reports chain-gap
evidence as deletion_START-END, insertion_LENGTH, or inversion_START-END.
Ranges use zero-based, half-open reference coordinates; insertion lengths are
source-gap lengths in bp after excluding the IS body itself. PAF CIGAR events
provide SV details when -p is supplied; --input-chain without -p uses chain-gap
evidence; gaps between split primary PAF records count as well. An SV on
either flank is reported when it touches the IS interval on the source, or lies
within --sv-flank-tol bp of the IS on the source and of the lifted position on
the reference. Inversion ranges are the aligned reference extent of the supporting
reverse PAF record (or chain), not inferred unaligned breakpoints.
Multiple details are separated by ';'; equal-size insertions on different
flanks remain separate entries. Translocation/rearrangement labels are
retained when present. '.' means no SV evidence reported.
For an encoded aligned IS body with an upstream reference deletion in
the same chain, both output endpoints use the upstream reference deletion
start (represented as a 1-bp site), matching the observable caller position.
Columns 5-6 (lifted_start/end) contain the zero-based mapped endpoint positions
in ascending order. If the mapped original start exceeds the mapped original
end, sv_type is set to 'inside inversion' and the positions in columns 5-6 are
swapped. Equal positions p retain the 1-bp representation p,p+1. Boundaries in
gaps use existing directional snapping. The BED output has one row per mapped input,
with positions sorted into a valid interval (or p,p+1 for an equal-position site).
A span between endpoints does not establish an IS body or a deletion.
  unmapped            No usable chain blocks were found for the source record.

Required programs:
  faSize, chainSort, chainPreNet, chainNet, netSyntenic, netChainSubset,
  chainStitchId, chainSwap, and python3. With --input-chain, only python3 is
  required.

Example:
  bash scripts/liftover_coords.sh \
    -p sample_to_reference.paf \
    -s sample.fasta \
    -r reference.fasta \
    -b sample_true.bed \
    --tag primary \
    --min-score 200 --min-space 10 --min-fill 10 \
    --chain sample_to_reference \
    -o results/

Chain-only example:
  bash scripts/liftover_coords.sh \
    -p sample_to_reference.paf \
    -s sample.fasta \
    -r reference.fasta \
    --chain \
    -o results/

Existing-chain example:
  bash scripts/liftover_coords.sh \
    --input-chain sample_to_reference.chain \
    -b sample_true.bed \
    -o results/
EOF
    exit "${1:-0}"
}

# -----------------------------------------------------------------------------
# Parse arguments
# -----------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --chain)
            if [[ $# -ge 2 && "$2" != -* ]]; then
                CHAIN_OUTPUT_NAME="$2"
                shift 2
            else
                CHAIN_OUTPUT_NAME=""
                shift
            fi
            ;;
        --chain=*)
            CHAIN_OUTPUT_NAME="${1#*=}"
            shift
            ;;
        -p|-s|-r|-b|-o|--min-score|--min-space|--min-fill|--tag|--input-chain|--sv-flank-tol)
            [[ $# -ge 2 ]] || {
                echo "ERROR: Option $1 requires an argument." >&2
                exit 1
            }
            case "$1" in
                -p) INPUT_PAF="$2" ;;
                -s) SAMPLE_FASTA="$2" ;;
                -r) REF_FASTA="$2" ;;
                -b) INPUT_BED="$2" ;;
                -o) OUTPUT_DIR="$2" ;;
                --min-score) CHAINNET_MIN_SCORE="$2" ;;
                --min-space) CHAINNET_MIN_SPACE="$2" ;;
                --min-fill)  CHAINNET_MIN_FILL="$2" ;;
                --tag)       PAF_TAG_MODE="$2" ;;
                --input-chain) INPUT_CHAIN="$2" ;;
                --sv-flank-tol) SV_FLANK_TOL="$2" ;;
            esac
            shift 2
            ;;
        -h|--help)
            usage 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "Run with -h or --help for usage." >&2
            exit 1
            ;;
    esac
done

# --tag also filters SV-evidence PAF records in --input-chain mode, so always
# validate it.
if [[ "$PAF_TAG_MODE" != "all" && "$PAF_TAG_MODE" != "primary" ]]; then
    echo "ERROR: --tag must be 'all' or 'primary'; got: $PAF_TAG_MODE" >&2
    exit 1
fi
if [[ -z "$INPUT_CHAIN" && ! "$CHAINNET_MIN_SCORE" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: --min-score must be a non-negative number; got: $CHAINNET_MIN_SCORE" >&2
    exit 1
fi
if [[ ! "$SV_FLANK_TOL" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --sv-flank-tol must be a non-negative integer; got: $SV_FLANK_TOL" >&2
    exit 1
fi
for option_and_value in \
    "--min-space:$CHAINNET_MIN_SPACE" \
    "--min-fill:$CHAINNET_MIN_FILL"
do
    option="${option_and_value%%:*}"
    value="${option_and_value#*:}"
    if [[ -z "$INPUT_CHAIN" && ! "$value" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $option must be a non-negative integer; got: $value" >&2
        exit 1
    fi
done

# -----------------------------------------------------------------------------
# Validate required arguments
# -----------------------------------------------------------------------------
INPUT_BED="${INPUT_BED:-}"
if [[ -n "$INPUT_CHAIN" ]]; then
    [[ -f "$INPUT_CHAIN" ]] || { echo "ERROR: File not found: $INPUT_CHAIN" >&2; exit 1; }
    [[ -n "$INPUT_BED" ]] || {
        echo "ERROR: -b <input_bed> is required with --input-chain." >&2
        exit 1
    }
    [[ -f "$INPUT_BED" ]] || { echo "ERROR: File not found: $INPUT_BED" >&2; exit 1; }
    if [[ -n "${INPUT_PAF:-}" && ! -f "$INPUT_PAF" ]]; then
        echo "ERROR: File not found: $INPUT_PAF" >&2
        exit 1
    fi
    ABS_INPUT_CHAIN=$(realpath "$INPUT_CHAIN")
    ABS_INPUT_BED=$(realpath "$INPUT_BED")
    ABS_INPUT_PAF=""
    [[ -n "${INPUT_PAF:-}" ]] && ABS_INPUT_PAF=$(realpath "$INPUT_PAF")
    ABS_SAMPLE="<not required>"
    ABS_REF="<not required>"
else
    if [[ -z "${INPUT_PAF:-}" || -z "${SAMPLE_FASTA:-}" || -z "${REF_FASTA:-}" ]]; then
        echo "ERROR: -p <input_paf>, -s <sample_fasta>, and -r <reference_fasta> are required unless --input-chain is used." >&2
        echo "Run with -h for usage." >&2
        exit 1
    fi
    for f in "$INPUT_PAF" "$SAMPLE_FASTA" "$REF_FASTA"; do
        [[ -f "$f" ]] || { echo "ERROR: File not found: $f" >&2; exit 1; }
    done
    if [[ -n "$INPUT_BED" && ! -f "$INPUT_BED" ]]; then
        echo "ERROR: File not found: $INPUT_BED" >&2
        exit 1
    fi
    ABS_INPUT_PAF=$(realpath "$INPUT_PAF")
    ABS_SAMPLE=$(realpath "$SAMPLE_FASTA")
    ABS_REF=$(realpath "$REF_FASTA")
    ABS_INPUT_BED=""
    [[ -n "$INPUT_BED" ]] && ABS_INPUT_BED=$(realpath "$INPUT_BED")
fi

# -----------------------------------------------------------------------------
# Derive base names
# -----------------------------------------------------------------------------
if [[ -n "$INPUT_CHAIN" ]]; then
    paf_label=$(basename "$INPUT_CHAIN" | sed 's/\.chain$//')
    sample="source"
    if [[ -n "${REF_FASTA:-}" ]]; then
        reference=$(basename "$REF_FASTA" | sed 's/\.\(fa\|fasta\|fa\.gz\|fasta\.gz\|fna\|fna\.gz\)$//')
    else
        reference="target"
    fi
else
    paf_label=$(basename "$INPUT_PAF" | sed 's/\.paf$//')
    sample=$(basename "$SAMPLE_FASTA" | sed 's/\.\(fa\|fasta\|fa\.gz\|fasta\.gz\|fna\|fna\.gz\)$//')
    reference=$(basename "$REF_FASTA"  | sed 's/\.\(fa\|fasta\|fa\.gz\|fasta\.gz\|fna\|fna\.gz\)$//')
fi

# -----------------------------------------------------------------------------
# Check required tools
# -----------------------------------------------------------------------------
# NB: transanno was previously used for PAF→chain conversion but produces an
# invalid chain for some PAFs (it splits long M runs that span IS-like
# regions, dropping ~1353 bp on the q side and tripping chainSort). The
# conversion is now done inline in python3 below.
if [[ -n "$INPUT_CHAIN" ]]; then
    REQUIRED_TOOLS=(python3)
else
    REQUIRED_TOOLS=(
        faSize
        chainSort
        chainPreNet
        chainNet
        netSyntenic
        netChainSubset
        chainStitchId
        chainSwap
        python3
    )
fi
missing=()
for tool in "${REQUIRED_TOOLS[@]}"; do
    command -v "$tool" &>/dev/null || missing+=("$tool")
done
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: Missing required tool(s): ${missing[*]}" >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

PREFIX="${OUTPUT_DIR}/${paf_label}"
if [[ -n "$INPUT_CHAIN" ]]; then
    LIFTOVER_CHAIN="$ABS_INPUT_CHAIN"
else
    SAMPLE_SIZES="${OUTPUT_DIR}/${sample}.sizes"
    REF_SIZES="${OUTPUT_DIR}/${reference}.sizes"
    CHAIN="${PREFIX}.input.chain"
    SORTED_CHAIN="${PREFIX}.sorted.chain"
    PRENET_CHAIN="${PREFIX}.prenet.chain"
    NET="${PREFIX}.net"
    SYN_NET="${PREFIX}.syn.net"
    SOURCE_BEST_RAW_CHAIN="${PREFIX}.source_best.raw.chain"
    SOURCE_BEST_CHAIN="${PREFIX}.source_best.chain"
    SWAPPED_UNSORTED_CHAIN="${PREFIX}.source_best.swapped.unsorted.chain"
    SWAPPED_CHAIN="${PREFIX}.source_best.swapped.sorted.chain"
    RBEST_PRENET_CHAIN="${PREFIX}.rbest.prenet.chain"
    RBEST_NET="${PREFIX}.rbest.net"
    RBEST_SYN_NET="${PREFIX}.rbest.syn.net"
    RBEST_RAW_CHAIN="${PREFIX}.rbest.raw.chain"
    RBEST_CHAIN="${PREFIX}.rbest.chain"
    FINAL_UNSORTED_CHAIN="${PREFIX}.rbest.final.unsorted.chain"
    if [[ -z "$CHAIN_OUTPUT_NAME" ]]; then
        CHAIN_OUTPUT_NAME="${paf_label}.chain"
    else
        CHAIN_OUTPUT_NAME=$(basename "$CHAIN_OUTPUT_NAME")
        [[ "$CHAIN_OUTPUT_NAME" == *.chain ]] || CHAIN_OUTPUT_NAME+=".chain"
    fi
    LIFTOVER_CHAIN="${OUTPUT_DIR}/${CHAIN_OUTPUT_NAME}"
fi

if [[ -n "$INPUT_BED" ]]; then
    input_bed_base=$(basename "$INPUT_BED" .bed)
    OUTPUT_BED="${OUTPUT_DIR}/${input_bed_base}_to_${reference}.bed"
    UNMAPPED_BED="${OUTPUT_DIR}/${input_bed_base}_to_${reference}.unmapped"
    LIFT_REPORT="${OUTPUT_DIR}/${input_bed_base}_to_${reference}.lift_report.tsv"
fi

LOG="${OUTPUT_DIR}/liftover_$(date +%Y%m%d)_${paf_label}.log"

echo "============================================================" | tee "$LOG"
if [[ -n "$INPUT_CHAIN" ]]; then
    echo "Liftover pipeline: existing chain"                       | tee -a "$LOG"
    echo "Input chain   : $LIFTOVER_CHAIN"                         | tee -a "$LOG"
else
    echo "Liftover pipeline: existing PAF"                         | tee -a "$LOG"
    echo "PAF file      : $ABS_INPUT_PAF"                          | tee -a "$LOG"
    echo "Sample FASTA  : $ABS_SAMPLE"                             | tee -a "$LOG"
    echo "Reference     : $ABS_REF"                                | tee -a "$LOG"
fi
echo "Input BED     : ${ABS_INPUT_BED:-<not provided; chain-only mode>}" | tee -a "$LOG"
echo "Chain file    : $LIFTOVER_CHAIN"                              | tee -a "$LOG"
echo "SV PAF        : ${ABS_INPUT_PAF:-<not supplied; chain-gap fallback>}" | tee -a "$LOG"
echo "Output dir    : $OUTPUT_DIR"                                  | tee -a "$LOG"
[[ -n "$INPUT_BED" ]] && \
    echo "SV flank tol  : $SV_FLANK_TOL bp"                         | tee -a "$LOG"
[[ -n "$ABS_INPUT_PAF" ]] && \
    echo "PAF tag mode  : $PAF_TAG_MODE"                            | tee -a "$LOG"
if [[ -z "$INPUT_CHAIN" ]]; then
    echo "chainNet opts : -minScore=$CHAINNET_MIN_SCORE -minSpace=$CHAINNET_MIN_SPACE -minFill=$CHAINNET_MIN_FILL" | tee -a "$LOG"
fi
echo "Log           : $LOG"                                         | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

# Helper: log + run a command whose output goes to the terminal/log
run() {
    echo "" | tee -a "$LOG"
    echo "[$(date +%H:%M:%S)] RUNNING: $*" | tee -a "$LOG"
    "$@" 2>&1 | tee -a "$LOG"
}

# Helper: log + run a command that writes its stdout to a file.
# Stderr is appended to the log only (not mixed into the output file).
#   run_to <output_file> <command> [args …]
run_to() {
    local outfile="$1"; shift
    echo "" | tee -a "$LOG"
    echo "[$(date +%H:%M:%S)] RUNNING: $* > $outfile" | tee -a "$LOG"
    "$@" > "$outfile" 2>> "$LOG"
}

if [[ -z "$INPUT_CHAIN" ]]; then
echo "" | tee -a "$LOG"
echo "--- Step 1: Chromosome sizes ---" | tee -a "$LOG"

run_to "$SAMPLE_SIZES" faSize -detailed "$ABS_SAMPLE"
run_to "$REF_SIZES"    faSize -detailed "$ABS_REF"

echo "" | tee -a "$LOG"
echo "--- Step 2: PAF to chain (inline python) ---" | tee -a "$LOG"

# Inline PAF→chain converter — replaces `transanno minimap2chain`. transanno
# was observed to corrupt long M runs spanning IS-like regions, emitting a
# chain whose cumulative q advance is short by the IS length (~1353 bp) and
# making chainSort reject the entry. This converter walks the CIGAR exactly
# once and emits the chain in the UCSC format the downstream tools expect:
#
#   chain.tName = PAF.qname (sample, source for liftover)
#   chain.qName = PAF.tname (reference, target for liftover)
#   Each <size dt dq> row = one M op plus the internal indels that immediately
#   follow it; the terminal row has just <size> (no trailing gaps).
#   Leading/trailing indels are trimmed into the chain header coordinates.
#   I in the CIGAR (sample-only advance) → dt; D (ref-only advance) → dq.
#
# strand handling: PAF's strand is between qname and tname; for '-' the q
# coordinates in the chain header are reverse-strand-mapped (qSize - qend,
# qSize - qstart) per UCSC chain convention.
[[ -f "$CHAIN" ]] && rm -f "$CHAIN"
ABS_INPUT_PAF="$ABS_INPUT_PAF" CHAIN_OUT="$CHAIN" PAF_TAG_MODE="$PAF_TAG_MODE" python3 <<'PYEOF' 2>&1 | tee -a "$LOG"
import os, re, sys
paf = os.environ['ABS_INPUT_PAF']
out = os.environ['CHAIN_OUT']
tag_mode = os.environ['PAF_TAG_MODE']
n_chains = 0
n_filtered_by_tag = 0
n_missing_cigar = 0
with open(out, 'w') as fh:
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12: continue
        tp_tag = next((tag[5:] for tag in f[12:] if tag.startswith('tp:A:')), None)
        if tag_mode == 'primary' and tp_tag != 'P':
            n_filtered_by_tag += 1
            continue
        qname, qlen, qst, qen, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
        tname, tlen, tst, ten         = f[5], int(f[6]), int(f[7]), int(f[8])
        score = f[9]                                  # number of matches
        cigar = ''
        for tag in f[12:]:
            if tag.startswith('cg:Z:'):
                cigar = tag[5:]; break
        if not cigar:
            n_missing_cigar += 1
            continue
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        # minimap2 writes the cg CIGAR for '-' PAF records in the
        # reverse-complement query frame. This chain uses PAF.qname as the
        # source/t sequence, whose coordinates must increase in the original
        # sample plus frame, so reverse the operation order before building
        # chain rows. I/D semantics are unchanged: I consumes sample, D/N
        # consumes reference.
        if strand == '-':
            ops = list(reversed(ops))
        first_m = next((i for i, (_, op) in enumerate(ops) if op in ('M', '=', 'X')), None)
        last_m = next((i for i in range(len(ops) - 1, -1, -1) if ops[i][1] in ('M', '=', 'X')), None)
        if first_m is None:
            continue
        n_chains += 1
        lead_t = lead_q = trail_t = trail_q = 0
        for n, op in ops[:first_m]:
            n = int(n)
            if op == 'I':
                lead_t += n
            elif op in ('D', 'N'):
                lead_q += n
        for n, op in ops[last_m + 1:]:
            n = int(n)
            if op == 'I':
                trail_t += n
            elif op in ('D', 'N'):
                trail_q += n
        # chain header — t = sample (PAF q), q = ref (PAF t). PAF qname has
        # an implicit '+' strand; PAF strand applies to the alignment of q
        # against t, so we map it onto chain's qStrand and reverse-complement
        # the q coordinates when '-' per UCSC convention.
        chain_t_st = qst + lead_t
        chain_t_en = qen - trail_t
        if strand == '+':
            q_lo, q_hi = tst + lead_q, ten - trail_q
            q_strand   = '+'
        else:
            q_lo, q_hi = tlen - ten + lead_q, tlen - tst - trail_q
            q_strand   = '-'
        fh.write(f"chain {score} {qname} {qlen} + {chain_t_st} {chain_t_en} "
                 f"{tname} {tlen} {q_strand} {q_lo} {q_hi} {n_chains}\n")
        # Ops are now trimmed to start and end with a match, so every row
        # starts with an M run and the last row carries no trailing gap.
        # Ops that consume neither sequence (S/H/P) are skipped. Consecutive
        # M runs (e.g. '=' followed by 'X') merge into one gap-free block.
        rows = []
        for n, op in ops[first_m:last_m + 1]:
            n = int(n)
            if op in ('M', '=', 'X'):
                if rows and rows[-1][1] == 0 and rows[-1][2] == 0:
                    rows[-1][0] += n
                else:
                    rows.append([n, 0, 0])
            elif op == 'I':
                rows[-1][1] += n
            elif op in ('D', 'N'):
                rows[-1][2] += n
        for r in rows[:-1]:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
        fh.write(f"{rows[-1][0]}\n\n")
print(f"Wrote {n_chains} chain entries → {out}")
print(f"PAF tag mode: {tag_mode}; filtered by tp tag: {n_filtered_by_tag}; "
      f"skipped without cg:Z: CIGAR: {n_missing_cigar}")
PYEOF

echo "" | tee -a "$LOG"
echo "--- Step 3: chainSort ---" | tee -a "$LOG"

run chainSort "$CHAIN" "$SORTED_CHAIN"

echo "" | tee -a "$LOG"
echo "--- Step 4: Reciprocal-best chain/net filtering ---" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Step 4a: source-side pre-net" | tee -a "$LOG"
run chainPreNet "$SORTED_CHAIN" "$SAMPLE_SIZES" "$REF_SIZES" "$PRENET_CHAIN"

echo "" | tee -a "$LOG"
echo "Step 4b: source-side net" | tee -a "$LOG"
run chainNet \
    "-minScore=$CHAINNET_MIN_SCORE" \
    "-minSpace=$CHAINNET_MIN_SPACE" \
    "-minFill=$CHAINNET_MIN_FILL" \
    "$PRENET_CHAIN" "$SAMPLE_SIZES" "$REF_SIZES" "$NET" /dev/null

echo "" | tee -a "$LOG"
echo "Step 4c: annotate source-side net synteny (no netFilter -syn)" | tee -a "$LOG"
run netSyntenic "$NET" "$SYN_NET"

echo "" | tee -a "$LOG"
echo "Step 4d: extract and stitch source-best chain" | tee -a "$LOG"
run netChainSubset "$SYN_NET" "$PRENET_CHAIN" "$SOURCE_BEST_RAW_CHAIN"
run chainStitchId "$SOURCE_BEST_RAW_CHAIN" "$SOURCE_BEST_CHAIN"

echo "" | tee -a "$LOG"
echo "Step 4e: swap source-best chain into reference orientation" | tee -a "$LOG"
run chainSwap "$SOURCE_BEST_CHAIN" "$SWAPPED_UNSORTED_CHAIN"
run chainSort "$SWAPPED_UNSORTED_CHAIN" "$SWAPPED_CHAIN"

echo "" | tee -a "$LOG"
echo "Step 4f: reference-side pre-net" | tee -a "$LOG"
run chainPreNet "$SWAPPED_CHAIN" "$REF_SIZES" "$SAMPLE_SIZES" "$RBEST_PRENET_CHAIN"

echo "" | tee -a "$LOG"
echo "Step 4g: reference-side net" | tee -a "$LOG"
run chainNet \
    "-minScore=$CHAINNET_MIN_SCORE" \
    "-minSpace=$CHAINNET_MIN_SPACE" \
    "-minFill=$CHAINNET_MIN_FILL" \
    "$RBEST_PRENET_CHAIN" "$REF_SIZES" "$SAMPLE_SIZES" "$RBEST_NET" /dev/null

echo "" | tee -a "$LOG"
echo "Step 4h: annotate reference-side net synteny (no netFilter -syn)" | tee -a "$LOG"
run netSyntenic "$RBEST_NET" "$RBEST_SYN_NET"

echo "" | tee -a "$LOG"
echo "Step 4i: extract, stitch, swap back final reciprocal-best chain" | tee -a "$LOG"
run netChainSubset "$RBEST_SYN_NET" "$RBEST_PRENET_CHAIN" "$RBEST_RAW_CHAIN"
run chainStitchId "$RBEST_RAW_CHAIN" "$RBEST_CHAIN"
run chainSwap "$RBEST_CHAIN" "$FINAL_UNSORTED_CHAIN"
run chainSort "$FINAL_UNSORTED_CHAIN" "$LIFTOVER_CHAIN"

if [[ -z "$INPUT_BED" ]]; then
    echo "" | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"
    echo "Done (chain-only mode)."                                    | tee -a "$LOG"
    echo "  Chain file: $LIFTOVER_CHAIN"                             | tee -a "$LOG"
    echo "  Full log  : $LOG"                                        | tee -a "$LOG"
    echo "============================================================" | tee -a "$LOG"
    exit 0
fi
fi

echo "" | tee -a "$LOG"
echo "--- Step 5: Chain-based direct liftover ---" | tee -a "$LOG"

export _LIFTOVER_CHAIN="$LIFTOVER_CHAIN"
export _INPUT_BED="$ABS_INPUT_BED"
export _OUTPUT_BED="$OUTPUT_BED"
export _UNMAPPED_BED="$UNMAPPED_BED"
export _LIFT_REPORT="$LIFT_REPORT"
export _OUTPUT_DIR="$OUTPUT_DIR"
export _LOG="$LOG"
export _SV_PAF="$ABS_INPUT_PAF"
export _PAF_TAG_MODE="$PAF_TAG_MODE"
export _SV_FLANK_TOL="$SV_FLANK_TOL"

python3 << 'PYEOF'
import os, sys, time, re
from collections import defaultdict

# ── Config from environment ────────────────────────────────────────────────────
liftover_chain = os.environ['_LIFTOVER_CHAIN']
input_bed_path = os.environ['_INPUT_BED']
output_bed     = os.environ['_OUTPUT_BED']
unmapped_bed   = os.environ['_UNMAPPED_BED']
lift_report    = os.environ['_LIFT_REPORT']
output_dir     = os.environ['_OUTPUT_DIR']
log_file       = os.environ['_LOG']

# ── Helpers ────────────────────────────────────────────────────────────────────
def log(msg):
    ts   = time.strftime('%H:%M:%S')
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(log_file, 'a') as fh:
        fh.write(line + '\n')

def read_bed(path):
    """Return list-of-lists for each non-header BED line."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#') \
               or line.startswith('track') or line.startswith('browser'):
                continue
            rows.append(line.split('\t'))
    return rows

def write_bed(rows, path):
    with open(path, 'w') as fh:
        for r in rows:
            fh.write('\t'.join(str(c) for c in r) + '\n')

def insertion_body_interval(extra_cols):
    if len(extra_cols) < 2:
        return None
    m = re.match(r'^(\d+)_(\d+)$', extra_cols[0])
    if m is None:
        return None
    # Encoded IS body intervals are 0-based half-open coordinates, matching BED
    # and PAF conventions.
    return int(m.group(1)), int(m.group(2)), extra_cols[1].upper()

def insertion_body_site_edge(extra_cols, orig_start, orig_end):
    iv = insertion_body_interval(extra_cols)
    if iv is None:
        return None
    q_start, q_end, ori = iv
    if q_end <= orig_start:
        return 'start'
    if q_start >= orig_end:
        return 'end'
    if ori == 'F':
        return 'start'
    if ori == 'R':
        return 'end'
    return None

def read_paf_sv_events(path):
    """Walk CIGAR in reference order, as plot_paf_region.cigar_events does.

    Also track source positions to associate each individual operation with
    either IS flank. Reverse PAF records consume the source from qEnd downward.

    Returns (gaps, inversions, junctions). junctions holds the gaps between
    consecutive primary records of a split alignment (same strand), in the same
    (source_lo, source_hi, ref_lo, ref_hi, strand, score) form as gaps; they are
    not CIGAR operations, so they are kept apart and used only for SV reporting.
    """
    if not path:
        return None, None, None
    gaps, inversions = defaultdict(list), defaultdict(list)
    primary = defaultdict(list)
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                raise ValueError(f'{path}:{number}: expected at least 12 PAF columns')
            tags = fields[12:]
            if os.environ['_PAF_TAG_MODE'] == 'primary' and 'tp:A:P' not in tags:
                continue
            cigar = next((tag[5:] for tag in tags if tag.startswith('cg:Z:')), None)
            if cigar is None:
                continue
            ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
            if ''.join(n + op for n, op in ops) != cigar:
                raise ValueError(f'{path}:{number}: malformed CIGAR')
            strand = fields[4]
            source = int(fields[2]) if strand == '+' else int(fields[3])
            target = int(fields[7])
            direction = 1 if strand == '+' else -1
            score = int(fields[9])
            key = (fields[0], fields[5])
            for n, op in ops:
                size = int(n)
                if op in 'M=X':
                    source += direction * size
                    target += size
                elif op in 'DN':
                    gaps[key].append((source, source, target, target + size, strand, score))
                    target += size
                elif op == 'I':
                    other = source + direction * size
                    gaps[key].append((min(source, other), max(source, other),
                                      target, target, strand, score))
                    source = other
            if strand == '-':
                inversions[key].append((int(fields[2]), int(fields[3]),
                                        int(fields[7]), int(fields[8])))
            if 'tp:A:P' in tags:
                primary[key].append((int(fields[2]), int(fields[3]), strand,
                                     int(fields[7]), int(fields[8]), score))
    junctions = defaultdict(list)
    for key, records in primary.items():
        records.sort()
        left = None
        for record in records:
            if left is not None and record[1] <= left[1]:
                continue                      # contained in the previous record
            if left is not None and record[2] == left[2]:
                qs1, qe1, strand, ts1, te1, s1 = left
                qs2, qe2, _, ts2, te2, s2 = record
                # Reference gap in alignment order; negative means the two
                # records overlap on the reference (e.g. a target-site
                # duplication), which is not a deletion.
                lo, hi = (te1, ts2) if strand == '+' else (te2, ts1)
                if hi < lo:
                    lo = hi = min(lo, hi)
                junctions[key].append((qe1, max(qe1, qs2), lo, hi, strand, min(s1, s2)))
            left = record
    return dict(gaps), dict(inversions), dict(junctions)


paf_sv_gaps, paf_sv_inversions, paf_sv_junctions = read_paf_sv_events(os.environ['_SV_PAF'])
sv_flank_tol = int(os.environ['_SV_FLANK_TOL'])


def flanking_chain_gaps(extra_cols, orig_chrom, orig_start, orig_end,
                        all_blocks, target_chrom):
    """Collect distinct chain gaps touching the BED interval or encoded body."""
    iv = insertion_body_interval(extra_cols)
    start = min(orig_start, iv[0]) if iv is not None else orig_start
    end = max(orig_end, iv[1]) if iv is not None else orig_end
    if paf_sv_gaps is not None:
        events = paf_sv_gaps.get((orig_chrom, target_chrom), [])
        # Include deletions adjacent in CIGAR order to an insertion spanning
        # either boundary, even when its source anchor is outside the IS body.
        changed = True
        while changed:
            old = start, end
            for a, b, lo, hi, strand, score in events:
                if a <= end and b >= start:
                    start, end = min(start, a), max(end, b)
            changed = old != (start, end)
        return [event for event in events
                if event[0] <= end and event[1] >= start]
    chains = defaultdict(list)
    for block in all_blocks:
        if block[2] == target_chrom:
            chains[(block[7], block[5])].append(block)
    gaps = {}
    for (_, strand), blocks in chains.items():
        blocks.sort(key=lambda b: b[0])
        for left, right in zip(blocks, blocks[1:]):
            lo, hi = (left[4], right[3]) if strand == '+' else (right[4], left[3])
            if right[0] < left[1] or hi < lo:
                continue
            if left[1] > end or right[0] < start:
                continue
            key = (left[1], right[0], lo, hi, strand)
            score = max(left[6], right[6])
            gaps[key] = max(gaps.get(key, 0), score)
    return [(*key, score) for key, score in sorted(gaps.items())]


def all_gap_events(orig_chrom, target_chrom, all_blocks):
    """Every gap event between this source and target chromosome.

    With a PAF: CIGAR indels plus gaps between consecutive primary records of a
    split alignment. Without: gaps between source-consecutive chain blocks on
    the same strand, whether or not they belong to the same chain ID.
    """
    if paf_sv_gaps is not None:
        key = (orig_chrom, target_chrom)
        return paf_sv_gaps.get(key, []) + paf_sv_junctions.get(key, [])
    events = []
    for strand in ('+', '-'):
        blocks = sorted(b for b in all_blocks
                        if b[2] == target_chrom and b[5] == strand)
        for left, right in zip(blocks, blocks[1:]):
            lo, hi = (left[4], right[3]) if strand == '+' else (right[4], left[3])
            if right[0] < left[1] or hi < lo:
                continue
            events.append((left[1], right[0], lo, hi, strand, max(left[6], right[6])))
    return events


def reference_adjacent_gaps(extra_cols, orig_chrom, orig_start, orig_end,
                            all_blocks, target_chrom, ref_lo, ref_hi):
    """Gap events adjacent to the IS on both the source and the reference.

    An event qualifies when its source span lies within sv_flank_tol bp of the
    BED interval/encoded body and its reference span lies within sv_flank_tol bp
    of the lifted reference footprint [ref_lo, ref_hi]. This catches SVs on
    either flank whose alignment anchor sits a few bp off the IS boundary, and
    SVs between split-alignment records, which the source-touching rule misses.
    """
    iv = insertion_body_interval(extra_cols)
    start = min(orig_start, iv[0]) if iv is not None else orig_start
    end = max(orig_end, iv[1]) if iv is not None else orig_end
    tol = sv_flank_tol
    return [event for event in all_gap_events(orig_chrom, target_chrom, all_blocks)
            if event[0] <= end + tol and event[1] >= start - tol
            and event[2] <= ref_hi + tol and event[3] >= ref_lo - tol]


def chain_gap_sv_types(extra_cols, orig_chrom, orig_start, orig_end,
                       all_blocks, target_chrom, ref_lo, ref_hi):
    """Classify source/reference gaps on either flank of the IS.

    Candidates are gaps touching the IS interval on the source plus every gap
    adjacent to it on both source and reference (reference_adjacent_gaps).
    Reference gaps indicate deletions. Source gaps touching the body may contain
    flanking insertions, but the IS body's own bases are excluded from their size.
    """
    iv = insertion_body_interval(extra_cols)
    start, end = iv[:2] if iv is not None else (orig_start, orig_end)
    events = set(flanking_chain_gaps(extra_cols, orig_chrom, orig_start, orig_end,
                                     all_blocks, target_chrom))
    events.update(reference_adjacent_gaps(extra_cols, orig_chrom, orig_start, orig_end,
                                          all_blocks, target_chrom, ref_lo, ref_hi))
    details = set()
    for source_lo, source_hi, lo, hi, strand, score in events:
        # Keep source coordinates in the identity: equal-length insertions on
        # opposite flanks are distinct events and must both appear in the TSV.
        for a, b in ((source_lo, min(source_hi, start)),
                     (max(source_lo, end), source_hi)):
            if b > a:
                details.add(('insertion', a, b, lo, hi))
        if hi > lo:
            details.add(('deletion', lo, hi, lo, hi))
    return [f'insertion_{hi - lo}' if kind == 'insertion' else f'deletion_{lo}-{hi}'
            for kind, lo, hi, _, _ in sorted(details)]


def inversion_sv_details(extra_cols, orig_chrom, orig_start, orig_end,
                         all_blocks, target_chrom):
    """Report the reference extent of reverse chains supporting this IS.

    Reverse records within sv_flank_tol bp of either flank count, so an
    inversion starting just past an IS boundary is reported too.
    """
    iv = insertion_body_interval(extra_cols)
    start = min(orig_start, iv[0]) if iv is not None else orig_start
    end = max(orig_end, iv[1]) if iv is not None else orig_end
    start, end = start - sv_flank_tol, end + sv_flank_tol
    if paf_sv_inversions is not None:
        ranges = {(lo, hi) for a, b, lo, hi in
                  paf_sv_inversions.get((orig_chrom, target_chrom), [])
                  if a <= end and b >= start}
        return [f'inversion_{lo}-{hi}' for lo, hi in sorted(ranges)]
    chains = defaultdict(list)
    for block in all_blocks:
        if block[2] == target_chrom and block[5] == '-':
            chains[block[7]].append(block)
    ranges = set()
    for blocks in chains.values():
        if any(b[0] <= end and b[1] >= start for b in blocks):
            ranges.add((min(b[3] for b in blocks), max(b[4] for b in blocks)))
    return [f'inversion_{lo}-{hi}' for lo, hi in sorted(ranges)]


def upstream_flanking_deletion_site(extra_cols, orig_chrom, orig_start, orig_end,
                                    all_blocks, target_chrom):
    """Return a deletion start on the lower-reference-coordinate IS flank.

    The IS boundary may lie in a source gap, or aligned flank sequence may
    separate it from a deletion inside the BED interval. Exact block/body edge
    equality and full alignment of the body are not required.
    """
    iv = insertion_body_interval(extra_cols)
    start, end = iv[:2] if iv is not None else (orig_start, orig_end)
    if end <= start:
        return None
    candidates = []
    gaps = flanking_chain_gaps(extra_cols, orig_chrom, orig_start, orig_end,
                               all_blocks, target_chrom)
    left_edge, right_edge = start, end
    for a, b, lo, hi, strand, score in gaps:
        if a <= start <= b:
            left_edge = max(left_edge, b)
        if a <= end <= b:
            right_edge = min(right_edge, a)
    for source_lo, source_hi, lo, hi, strand, score in gaps:
        upstream = source_lo <= left_edge if strand == '+' else source_hi >= right_edge
        if hi > lo and upstream:
            candidates.append((score, lo))
    return min(candidates, key=lambda item: (-item[0], item[1]))[1] if candidates else None


def lift_is_body_to_ref(extra_cols, blocks):
    """
    Project the sample IS body interval (q_start, q_end) through the chain
    to its corresponding reference position, then return a 1-bp insertion
    site. This is a post-rescue interpretation step: minimap2/flank rescue
    has already produced the PAF CIGAR, and the PAF has already been converted
    to chain blocks. Unifies three cases that previously needed separate code
    paths:

      (a) Same-site / column-aligned — sample IS body fits inside a single
          M block (e.g. minimap2 aligned the sample IS to the ref IS as M).
          block_map projects q_start to its ref column; that column is the
          site. This is the case that the old flank-based logic missed
          because the containing block neither ends before q_start nor
          starts after q_end.

      (b) Pure insertion — sample IS body sits in a query-only gap between
          two flanking blocks (no ref deletion). The gap collapses to a
          single ref edge: qe of the left flank block (= qs of the right
          flank block) on forward, qs of the left flank block on reverse.

      (c) Insertion + ref deletion — sample IS body sits in a query gap
          adjacent, in alignment order, to a reference-side gap. The sample
          IS is reported at the adjacent deletion edge, not at the local
          pure-insertion edge. The edge choice depends on the chain strand;
          this is independent of the biological F/R orientation of the IS.

    No reference-IS BED is consulted; the chain alone determines the site.

    Returns (site, orientation_letter, context) where context is
    'insertion', a list of inversion-boundary sites when the two adjacent
    flanks disagree in strand, or None when neither a containing block nor
    a flanking block can be found.
    """
    iv = insertion_body_interval(extra_cols)
    if iv is None:
        return None
    q_start, q_end, ori = iv

    # (a) Same-site: the sample IS body overlaps an M block. Use the mapped
    #     sample-interval boundary, not a ref-side annotation. On '+' the lower
    #     insertion edge is the mapped body start. On '-' the lower insertion
    #     edge is one base after the lowest mapped body base, i.e. the boundary
    #     after body_end - 1.
    for ts, te, qs, qe, strand in blocks:
        if ts <= q_start < te:
            body_start = max(q_start, ts)
            body_end = min(q_end, te)
            if body_end <= body_start:
                site = block_map(q_start, ts, te, qs, qe, strand)
                return site, ori, 'insertion'
            if strand == '+':
                site = block_map(body_start, ts, te, qs, qe, strand)
            else:
                site = block_map(body_end - 1, ts, te, qs, qe, strand) + 1
            return site, ori, 'insertion'

    # (b)/(c) IS body lies in a query-side gap. Find the left flank block
    # (largest te <= q_start) and right flank block (smallest ts >= q_end).
    left_candidates  = [b for b in blocks if b[1] <= q_start]
    right_candidates = [b for b in blocks if b[0] >= q_end]
    left  = max(left_candidates,  key=lambda b: b[1]) if left_candidates  else None
    right = min(right_candidates, key=lambda b: b[0]) if right_candidates else None

    if left is None and right is None:
        return None
    if left is not None and right is not None and left[4] != right[4]:
        sites = []
        # Opposite-strand flanks mark an inversion/rearrangement boundary.
        # Keep the two breakends separate; min/max would invent a continuous
        # reference interval spanning unrelated sides of the inversion.
        left_site = left[3] if left[4] == '+' else left[2]
        right_site = right[2]
        sites.append((left_site, ori, 'inversion_boundary_left'))
        sites.append((right_site, ori, 'inversion_boundary_right'))
        return sites

    if left is not None and right is not None:
        # The sample IS lies in an unaligned sample gap. The gap may span a
        # pure insertion, insertion plus deletion, or a rearrangement-adjacent
        # interval. The exact internal insertion coordinate is not observable,
        # so report the lower reference breakend edge consistently.
        left_edge = left[3] if left[4] == '+' else left[2]
        right_edge = right[2] if right[4] == '+' else right[3]
        site = min(left_edge, right_edge)
    else:
        anchor = left if left is not None else right
        if anchor[4] == '+':
            site = left[3] if left is not None else right[2]
        else:
            site = left[2] if left is not None else right[3]
    return site, ori, 'insertion'

# ── Parse liftover chain ───────────────────────────────────────────────────────
def parse_liftover_chain(chain_path):
    """
    Parse a UCSC chain file.

    Required orientation:
      t = SAMPLE  (source genome — the assembly being lifted from)
      q = REF     (target genome — the reference being lifted to)

    Returns dict: t_chrom → sorted list of
        (t_start, t_end, q_chrom, q_start, q_end, q_strand, chain_score, chain_id)
         ↑ source coords ↑         ↑ target coords ↑
    All coordinates are 0-based half-open in the PLUS-STRAND frame.
    chain_score is stored for use as a tiebreaker in target-chromosome selection.
    chain_id is the UCSC chain record ID; blocks from different chain IDs but
    anchoring the same source interval indicate a structural variant boundary.
    """
    blocks = defaultdict(list)

    with open(chain_path) as fh:
        chain_score = 0
        q_chrom = q_size = q_strand = None
        t_chrom = None
        q_pos = t_pos = 0

        for raw in fh:
            line = raw.strip()
            if not line:
                q_chrom = None        # blank line → chain ended
                continue

            if line.startswith('chain'):
                # chain <score> <tName> <tSize> <tStrand> <tStart> <tEnd>
                #               <qName> <qSize> <qStrand> <qStart> <qEnd> <id>
                p           = line.split()
                chain_score = int(p[1])
                t_chrom     = p[2]
                t_pos       = int(p[5])        # tStart (always + strand, 0-based)
                q_chrom     = p[7]
                q_size      = int(p[8])
                q_strand    = p[9]
                q_pos       = int(p[10])       # qStart in chain-strand coords
                chain_id    = p[12] if len(p) > 12 else p[1]   # chain record ID
                continue

            if q_chrom is None:
                continue

            p    = line.split()
            size = int(p[0])

            # Convert q-position to plus-strand genomic coordinates.
            # For '-' strand the chain accumulates q_pos from the rc-start:
            #   plus_start = q_size - (q_pos + size)
            #   plus_end   = q_size - q_pos
            if q_strand == '+':
                qs, qe = q_pos, q_pos + size
            else:
                qs = q_size - (q_pos + size)
                qe = q_size - q_pos

            # Key by t_chrom (source/sample); source coords first, target second.
            blocks[t_chrom].append(
                (t_pos, t_pos + size, q_chrom, qs, qe, q_strand, chain_score, chain_id)
            )

            if len(p) == 3:             # non-terminal block: advance both cursors
                q_pos += size + int(p[2])   # dq gap
                t_pos += size + int(p[1])   # dt gap
            # terminal block (len==1): q_chrom reset at next blank line

    for ch in blocks:
        blocks[ch].sort(key=lambda b: b[0])   # sort by source start (t_pos)

    return dict(blocks)


# ── Select target chromosome ───────────────────────────────────────────────────
def select_target_chrom(orig_start, orig_end, all_blocks):
    """
    Choose the best target chromosome for lifting [orig_start, orig_end].

    Priority (applied in order):
      1. Largest total source-bp overlap with [orig_start, orig_end]
      2. Highest sum of chain scores among overlapping blocks
      3. If no block overlaps the interval at all: smallest distance to interval

    Returns the chosen q_chrom string, or None if all_blocks is empty.
    """
    overlap_bp    = defaultdict(int)
    overlap_score = defaultdict(int)
    min_dist      = defaultdict(lambda: float('inf'))

    for (ts, te, q_chrom, qs, qe, strand, score, cid) in all_blocks:
        ovl = max(0, min(te, orig_end) - max(ts, orig_start))
        overlap_bp[q_chrom] += ovl
        if ovl > 0:
            overlap_score[q_chrom] += score
        dist = min(abs(ts - orig_start), abs(ts - orig_end),
                   abs(te - orig_start), abs(te - orig_end))
        if dist < min_dist[q_chrom]:
            min_dist[q_chrom] = dist

    if not overlap_bp:
        return None

    return min(overlap_bp,
               key=lambda ch: (-overlap_bp[ch], -overlap_score[ch], min_dist[ch]))


# ── Find the best anchor block for a single source position ───────────────────
def find_anchor_block(t_pos, all_blocks):
    """
    Return the single best block anchoring t_pos from all_blocks.

    Preference order:
      1. Blocks that contain t_pos (ts <= t_pos < te); highest chain_score wins.
      2. If t_pos is in a gap: nearest block by source distance; highest
         chain_score breaks ties.

    Used to identify which chain record anchors each boundary of a source
    interval, so that different chain IDs signal a structural variant boundary.

    Returns the full 8-tuple (ts, te, q_chrom, qs, qe, strand, score, chain_id),
    or None if all_blocks is empty.
    """
    if not all_blocks:
        return None

    inside = [b for b in all_blocks if b[0] <= t_pos < b[1]]
    if inside:
        return max(inside, key=lambda b: b[6])   # highest chain_score

    # Not inside any block: pick nearest by source distance, break ties by score
    return min(all_blocks,
               key=lambda b: (min(abs(b[0] - t_pos), abs(b[1] - t_pos)), -b[6]))


# ── Map a source coordinate to a target coordinate within one block ────────────
def block_map(t_pos, ts, te, qs, qe, strand):
    """
    Map source coordinate t_pos (inside [ts, te)) to a target coordinate.

    Blocks are internally gap-free (no indels), so mapping is a direct offset:
      '+' strand: source[ts + k]  →  target[qs + k]
      '-' strand: source[ts + k]  →  target[qe - 1 - k]

    Returns a 0-based inclusive integer target position.
    """
    k = t_pos - ts
    return qs + k if strand == '+' else qe - 1 - k


# ── Refine one boundary to the nearest chain block ────────────────────────────
def refine_boundary(t_pos, blocks, for_start):
    """
    Return the tightest real target coordinate that bounds t_pos.

    blocks   : list of (ts, te, qs, qe, strand) — source coords first,
               target coords second — sorted by ts, pre-filtered to one
               target chromosome.
    for_start: True  → refining the START boundary
               False → refining the END   boundary

    EXACT (t_pos is inside a block):
        Direct offset mapping via block_map; gap = 0.

    IN A GAP (t_pos falls between blocks):
        START → target end of the block immediately LEFT of the gap
                  '+': qe  (exclusive end = first pos past aligned region)
                  '-': qs
        END   → target start of the block immediately RIGHT of the gap
                  '+': qs
                  '-': qe
        gap = source-bp distance from t_pos to the nearest block edge.

        These two edges together give the tightest interval in the target
        guaranteed to contain the true (unaligned) position.

    Returns (target_coord, gap_bp), or None if blocks is empty.
    target_coord is a 0-based inclusive integer; caller adds +1 for the END
    to produce a half-open BED interval.
    """
    if not blocks:
        return None

    # ── Exact: t_pos is inside a block ────────────────────────────────────────
    for (ts, te, qs, qe, strand) in blocks:
        if ts <= t_pos < te:
            return block_map(t_pos, ts, te, qs, qe, strand), 0

    # ── t_pos is in a gap ─────────────────────────────────────────────────────
    # Blocks are sorted by ts only; with several chains on one target chrom
    # the block with the largest ts is not necessarily the one ending
    # closest to t_pos, so select by the relevant edge explicitly.
    left  = [b for b in blocks if b[1] <= t_pos]
    right = [b for b in blocks if b[0] >= t_pos]
    if for_start:
        # Use the target END of the block immediately to the LEFT of t_pos.
        if left:
            ts, te, qs, qe, strand = max(left, key=lambda b: b[1])
            return (qe if strand == '+' else qs), (t_pos - te)
        # t_pos is before all blocks: snap to the start of the first block.
        ts, te, qs, qe, strand = min(right, key=lambda b: b[0])
        return (qs if strand == '+' else qe), (ts - t_pos)

    else:  # for_start=False  →  END boundary
        # Use the target START of the block immediately to the RIGHT of t_pos.
        if right:
            ts, te, qs, qe, strand = min(right, key=lambda b: b[0])
            return (qs if strand == '+' else qe), (ts - t_pos)
        # t_pos is after all blocks: snap to the end of the last block.
        ts, te, qs, qe, strand = max(left, key=lambda b: b[1])
        return (qe if strand == '+' else qs), (t_pos - te)


# ── Read input BED ─────────────────────────────────────────────────────────────
raw_rows = read_bed(input_bed_path)
if not raw_rows:
    log("ERROR: Input BED is empty or missing.")
    sys.exit(1)

log(f"    {len(raw_rows)} input record(s).")

# ── Parse liftover chain ───────────────────────────────────────────────────────
log("    Parsing liftover chain …")
chain_index = parse_liftover_chain(liftover_chain)
total_blks  = sum(len(v) for v in chain_index.values())
log(f"    Loaded {total_blks} alignment blocks across "
    f"{len(chain_index)} source sequence(s).")
if total_blks == 0:
    log("    WARNING: chain_index is EMPTY — chain file may have failed to parse.")

# ── Step 6: Chain-based direct liftover ───────────────────────────────────────
log("")
log("=== Step 5: Chain-based direct liftover ===")

lifted        = {}   # row_idx → one BED row, stored in a list
report_positions = {}  # row_idx → mapped source-start/source-end positions
gap_info      = {}   # row_idx → (gap_start_bp, gap_end_bp)
n_blocks_used = {}   # row_idx → number of blocks on selected target chrom
status_map    = {}   # row_idx → exact | gap_snapped | insertion_site | inversion_boundary | unmapped
sv_map        = {}   # row_idx → None or ';'-joined SV details
insert_orientation_map = {}  # row_idx → forward | reverse | .
insert_context_map     = {}  # row_idx → insertion | inversion_boundary | .

status_counts = defaultdict(int)

for idx, row in enumerate(raw_rows):
    orig_chrom = row[0]
    orig_start = int(row[1])
    orig_end   = int(row[2])
    extra_cols = row[3:]
    insert_orientation = '.'
    insert_context = '.'

    all_blocks = chain_index.get(orig_chrom, [])

    if not all_blocks:
        status_map[idx] = 'unmapped'
        status_counts['unmapped'] += 1
        log(f"    [unmapped] {orig_chrom}:{orig_start}-{orig_end} — "
            f"no chain blocks on this chromosome.")
        continue

    # ── Detect structural variant by comparing anchor blocks for each boundary ──
    # If orig_start and orig_end are anchored by blocks from different chain IDs,
    # the source interval spans a structural variant breakpoint.
    # Chain ID is the discriminator: same chain → indel (normal); different → SV.
    anchor_s = find_anchor_block(orig_start, all_blocks)
    anchor_e = find_anchor_block(orig_end,   all_blocks)
    sv_type  = None
    if anchor_s is not None and anchor_e is not None:
        cid_s, qc_s, strand_s = anchor_s[7], anchor_s[2], anchor_s[5]
        cid_e, qc_e, strand_e = anchor_e[7], anchor_e[2], anchor_e[5]
        if cid_s != cid_e:
            if qc_s != qc_e:
                sv_type = 'translocation'    # ends map to different chromosomes
            elif strand_s != strand_e:
                sv_type = 'inversion'        # same chrom, opposite strand
            else:
                sv_type = 'rearrangement'    # same chrom, same strand, diff chain

    # Select the best target chromosome by overlap then chain score
    target_chrom = select_target_chrom(orig_start, orig_end, all_blocks)
    if target_chrom is None:
        status_map[idx] = 'unmapped'
        status_counts['unmapped'] += 1
        continue

    # Filter to blocks targeting the selected chromosome; strip chain_score/chain_id
    blks = [(ts, te, qs, qe, strand)
            for (ts, te, qc, qs, qe, strand, score, cid) in all_blocks
            if qc == target_chrom]
    # Project the sample IS body's sample-coord interval through the chain
    # to its ref-coord position. This is the unified replacement for the
    # old insertion_flank_site + deletion_gap_site_for_insertion pair —
    # handles same-site (IS body in an M block), pure insertion (IS body
    # in a query-only gap), and insertion+ref-deletion (IS body in a query
    # gap with a ref gap nearby) under one rule.
    is_site = lift_is_body_to_ref(extra_cols, blks)

    res_start = refine_boundary(orig_start, blks, for_start=True)
    res_end   = refine_boundary(orig_end,   blks, for_start=False)

    if res_start is None or res_end is None:
        status_map[idx] = 'unmapped'
        status_counts['unmapped'] += 1
        continue

    t_start, gap_start = res_start
    t_end,   gap_end   = res_end
    # Reference footprint of the mapped boundaries, before any collapse to a
    # deletion site; SV reporting looks for events adjacent to it.
    ref_lo, ref_hi = min(t_start, t_end), max(t_start, t_end)

    deletion_site = upstream_flanking_deletion_site(extra_cols, orig_chrom,
                                                    orig_start, orig_end,
                                                    all_blocks, target_chrom)
    if deletion_site is not None:
        # The upstream flank identifies a deletion breakpoint. Use that breakpoint
        # as the observable truth site on either chain strand.
        t_start = t_end = deletion_site

    stat = 'exact' if gap_start == 0 and gap_end == 0 else 'gap_snapped'
    insertion_edge = insertion_body_site_edge(extra_cols, orig_start, orig_end)
    if isinstance(is_site, list):
        stat = 'inversion_boundary'
        insert_orientation = 'reverse' if is_site[0][1] == 'R' else 'forward'
        insert_context = 'inversion_boundary'
        sv_type = 'inversion'
    elif is_site is not None and is_site[2] == 'insertion':
        stat = 'insertion_site'
        insert_orientation = 'reverse' if is_site[1] == 'R' else 'forward'
        insert_context = 'insertion'
    elif insertion_edge is not None:
        # Retain both mapped BED boundaries in the insertion fallback too.
        stat = 'insertion_site'
        iv = insertion_body_interval(extra_cols)
        if iv is not None:
            insert_orientation = 'reverse' if iv[2] == 'R' else 'forward'
            insert_context = 'insertion'
    # Keep the mapped boundaries (or the upstream deletion breakpoint) together.
    # Do not replace them with the IS-body representative or split them into
    # alternative truth locations. Equal positions retain the 1-bp convention.
    if t_start > t_end:
        sv_type = 'inside inversion'
        report_start, report_end = t_end, t_start
    else:
        report_start = t_start
        report_end = t_end if t_end != t_start else t_start + 1
    gap_sv_types = chain_gap_sv_types(extra_cols, orig_chrom, orig_start, orig_end,
                                     all_blocks, target_chrom,
                                     min(ref_lo, report_start), max(ref_hi, report_end))
    sv_details = [sv_type] if sv_type else []
    if paf_sv_inversions is not None or sv_type in {'inversion', 'inside inversion'}:
        inversion_details = inversion_sv_details(extra_cols, orig_chrom, orig_start, orig_end,
                                                 all_blocks, target_chrom)
        if inversion_details:
            sv_details = ([] if sv_type in {'inversion', 'inside inversion'} else sv_details) + inversion_details
    sv_type = ';'.join(sv_details + gap_sv_types) or None
    report_positions[idx] = (report_start, report_end)
    new_start = min(t_start, t_end)
    new_end = max(t_start, t_end) if t_start != t_end else t_start + 1
    lifted_rows = [[target_chrom, str(new_start), str(new_end)] + extra_cols]

    lifted[idx]        = lifted_rows
    gap_info[idx]      = (gap_start, gap_end)
    n_blocks_used[idx] = len(blks)
    status_map[idx]    = stat
    sv_map[idx]        = sv_type
    insert_orientation_map[idx] = insert_orientation
    insert_context_map[idx] = insert_context

    status_counts[stat] += 1

    sv_tag = f"  sv={sv_type}" if sv_type else ''
    log(f"    [{stat}] {orig_chrom}:{orig_start}-{orig_end} → "
        f"{target_chrom}:{report_start},{report_end}  "
        f"start_gap={gap_start}bp  end_gap={gap_end}bp  "
        f"blocks={len(blks)}{sv_tag}")

log('    ' + '  '.join(
    f"{name}={status_counts[name]}"
    for name in ('exact', 'gap_snapped', 'insertion_site',
                 'inversion_boundary', 'unmapped')))

# ── Write outputs ──────────────────────────────────────────────────────────────
final_rows    = [r for i in range(len(raw_rows)) if i in lifted for r in lifted[i]]
unmapped_rows = [raw_rows[i] for i in range(len(raw_rows)) if i not in lifted]

write_bed(final_rows, output_bed)
log(f"\n    Written {len(final_rows)} lifted record(s) → {output_bed}")

if unmapped_rows:
    write_bed(unmapped_rows, unmapped_bed)
    log(f"    Written {len(unmapped_rows)} unmapped record(s) → {unmapped_bed}")
elif os.path.exists(unmapped_bed):
    # Drop a stale file from an earlier run so the summary does not list it.
    os.remove(unmapped_bed)

# Lift report TSV
with open(lift_report, 'w') as fh:
    fh.write('\t'.join([
        'original_chrom', 'original_start', 'original_end',
        'lifted_chrom',   'lifted_start',   'lifted_end',
        'start_gap_to_nearest_block_bp',
        'end_gap_to_nearest_block_bp',
        'insert_orientation', 'insert_context',
        'n_chain_blocks', 'status', 'original_bed_col6', 'sv_type'
    ]) + '\n')
    for idx, row in enumerate(raw_rows):
        oc, os_, oe_ = row[0], row[1], row[2]
        original_bed_col6 = row[5] if len(row) > 5 else '.'
        if idx in lifted:
            lifted_rows = lifted[idx]
            gs, ge = gap_info[idx]
            nb     = str(n_blocks_used[idx])
            stat   = status_map[idx]
            sv     = sv_map.get(idx) or '.'
            ins_ori = insert_orientation_map.get(idx, '.')
            ins_ctx = insert_context_map.get(idx, '.')
            report_start, report_end = report_positions[idx]
            fh.write('\t'.join([oc, os_, oe_,
                                lifted_rows[0][0], str(report_start), str(report_end),
                                str(gs), str(ge), ins_ori, ins_ctx, nb, stat,
                                original_bed_col6, sv]) + '\n')
        else:
            fh.write('\t'.join([oc, os_, oe_,
                                '.', '.', '.',
                                '.', '.', '.', '.', '.', 'unmapped',
                                original_bed_col6, '.']) + '\n')

log(f"    Lift report → {lift_report}")
PYEOF

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "Done." | tee -a "$LOG"
echo "  Lifted BED  : $OUTPUT_BED"      | tee -a "$LOG"
echo "  Lift report : $LIFT_REPORT"     | tee -a "$LOG"
[[ -f "$UNMAPPED_BED" ]] && \
    echo "  Unmapped BED: $UNMAPPED_BED"    | tee -a "$LOG"
echo "  Chain file  : $LIFTOVER_CHAIN"  | tee -a "$LOG"
echo "  Full log    : $LOG"             | tee -a "$LOG"

# Count records the same way read_bed does: skip blank, '#', track, browser.
total=$(grep -cvE '^(#|track|browser|$)' "$ABS_INPUT_BED" 2>/dev/null || true)
lifted=$(wc -l < "$OUTPUT_BED" 2>/dev/null || echo "?")
echo "" | tee -a "$LOG"
echo "  Input regions : $total"  | tee -a "$LOG"
echo "  Lifted regions: $lifted" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
