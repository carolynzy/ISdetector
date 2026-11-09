#!/usr/bin/env python3
"""
Random genome variation simulator.

The simulator accepts command-line parameters and/or a JSON template. When a config is supplied, all other command-line
parameters are ignored. Coordinates in truth outputs are 1-based and
refer to the original input reference. Dry runs save novel sequence fragments and source recipes in gzip tables;
--simulate LOG uses their replay data to construct genomes without planning new events.
Independent placements keep their nearest occupied reference bases at least 1,000 bp apart.
Preferred BED placements and boundary SNPs are exempt; existing inversion/SV interactions
and linked components within a compound IS placement remain supported.
Placement searches try at most 100 candidates. Three consecutive failed events stop
their category for the replicate; a successful event resets the failure count.
"""
from __future__ import annotations

import argparse
import copy
import contextlib
import csv
import hashlib
import gzip
import io
import platform
import shlex
import time
import datetime as dt
import json
import math
import random
import sys
import textwrap
from bisect import bisect_right
from collections import UserDict
from dataclasses import dataclass
from pathlib import Path


BASES = "ACGT"
DEFAULT_FASTA_COMPRESSION_THRESHOLD = 100_000_000
DEFAULT_BATCH_START = 1
PLANNED_EVENT_MIN_DISTANCE = 1000
MAX_PLACEMENT_ATTEMPTS = 100
MAX_CONSECUTIVE_EVENT_FAILURES = 3
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
SV_ALIASES = {
    "del": "deletion",
    "deletion": "deletion",
    "ins": "insertion",
    "insertion": "insertion",
    "inv": "inversion",
    "inversion": "inversion",
    "invertion": "inversion",
    "dup": "duplication",
    "duplication": "duplication",
    "trans": "translocation",
    "translocation": "translocation",
}
TRANSITIONS = {"A": "G", "G": "A", "C": "T", "T": "C"}

DEFAULTS = {
    "input": None,
    "outdir": "./simulation",
    "prefix": "simulation",
    "batch": 1,
    "seed": 378,
    "modes": {"complex_mode": False, "is_enabled": True, "snp_enabled": False, "sv_enabled": False},
    "snp": {
        "frequency": 0.0001,
        "transition": 2.0,
        "transversion": 1.0,
        "boundary_snp_frequency": 2,
        "boundary_snp_window": 50,
    },
    "sv": {
        "specs": [
            {"subtype": "deletion", "size_range": "100-500", "frequency": 5},
            {"subtype": "insertion", "size_range": "100-500", "frequency": 5},
            {"subtype": "inversion", "size_range": "100-500", "frequency": 5},
            {"subtype": "duplication", "size_range": "100-500", "frequency": 5},
            {"subtype": "translocation", "size_range": "100-500", "frequency": 5, "cross_chromosome_count": 2},
        ],
        "long_insertion_min_size": 100,
        "long_insertion_candidates_fasta": None,
        "long_insertion_local_window": 100000,
        "long_insertion_k": 5,
    },
    "reference_is": {
        "known_is_bed": None,
        "reference_removal_frequency_range": 0.5,
    },
    "preferred_bed_targeting": {
        "preferred_is_bed": None,
        "preferred_bed_target_proportion": 0.3,
        "preferred_bed_window": 10,
    },
    "is_elements": [
        {
            "name": None,
            "is_seq": None,
            "simple_insertion_frequency_range": "100-500",
            "partial_insertion_frequency_range": "100-500",
            "partial_insertion_proportion_range": "0.1-0.9",
            "nested_insertion_frequency_range": "100-500",
            "nested_mixed_percentage": 0.5,
            "reference_adjacent_insertion_percentage": 0.2,
            "reference_adjacent_mixed_percentage": 0.5,
            "consecutive_insertion_frequency_range": "100-500",
            "consecutive_mixed_percentage": 0.5,
            "nearby_insertion_frequency_range": "100-500",
            "nearby_insertion_gap_range": "1-100",
            "nearby_mixed_percentage": 0.5,
            "flanking_deletion_size_ranges": ["10-100", "101-500", "501-1000", "1001-5000", "5001-10000"],
            "flanking_deletion_frequency_ranges": ["20-100", "20-100", "20-100", "20-100", "20-100"],
            "flanking_insertion_size_ranges": ["10-20", "21-50", "51-100", "101-200", "201-500"],
            "flanking_insertion_frequency_ranges": ["20-100", "20-100", "20-100", "20-100", "20-100"],
            "flanking_inversion_size_ranges": ["101-500", "501-5000", "5001-50000", "50001-100000", "100001-1000000"],
            "flanking_inversion_frequency_ranges": ["20-100", "20-100", "1-5", "1-5", "0-1"],
        }
    ],
}


@dataclass
class Record:
    name: str
    desc: str
    seq: str


def deep_update(base, update):
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def validate_config_completeness(config):
    """Reject incomplete JSON before defaults or CLI overrides can hide omissions."""
    errors = []

    def check(value, template, path):
        if isinstance(template, dict):
            if not isinstance(value, dict):
                errors.append(f"{path or 'config'}: expected an object")
                return
            for key, default in template.items():
                field = f"{path}.{key}" if path else key
                if key not in value:
                    errors.append(f"{field}: missing required field")
                else:
                    check(value[key], default, field)
        elif isinstance(template, list):
            if not isinstance(value, list):
                errors.append(f"{path}: expected a list")
                return
            for index, item in enumerate(value):
                field = f"{path}[{index}]"
                if path == "sv.specs" and isinstance(item, str):
                    try:
                        parse_sv_spec(item)
                    except (ValueError, TypeError, argparse.ArgumentTypeError) as exc:
                        errors.append(f"{field}: {exc}")
                else:
                    check(item, template[0], field)
                    if path == "sv.specs" and isinstance(item, dict) and all(
                        key in item for key in ("subtype", "size_range", "frequency")
                    ):
                        try:
                            parse_sv_spec(item)
                        except (ValueError, TypeError, argparse.ArgumentTypeError) as exc:
                            errors.append(f"{field}: {exc}")
        elif value is None:
            if template is not None or path == "input":
                errors.append(f"{path}: must not be null")
        elif isinstance(value, str) and not value.strip():
            errors.append(f"{path}: must not be empty; use null for optional inputs")

    check(config, DEFAULTS, "")
    if isinstance(config, dict) and "batch_start" in config:
        try:
            parse_batch_start(config["batch_start"])
        except argparse.ArgumentTypeError as exc:
            errors.append(str(exc))
    if isinstance(config, dict) and "chromosome_weight" in config:
        try:
            normalize_chromosome_weights(config["chromosome_weight"])
        except ValueError as exc:
            errors.append(str(exc))
    if isinstance(config, dict) and "dryrun" in config and not isinstance(config["dryrun"], bool):
        errors.append("dryrun: expected a boolean")
    if isinstance(config, dict) and "debug" in config and not isinstance(config["debug"], bool):
        errors.append("debug: expected a boolean")
    if isinstance(config, dict) and isinstance(config.get("modes"), dict):
        for flag in ("complex_mode", "is_enabled", "snp_enabled", "sv_enabled"):
            if flag in config["modes"] and not isinstance(config["modes"][flag], bool):
                errors.append(f"modes.{flag}: expected a boolean")
    if isinstance(config, dict) and isinstance(config.get("is_elements"), list):
        for index, element in enumerate(config["is_elements"]):
            if not isinstance(element, dict):
                continue
            if element.get("is_seq") is not None:
                name = element.get("name")
                if not isinstance(name, str) or not name.strip():
                    errors.append(f"is_elements[{index}].name: required nonempty string when is_seq is provided")
            for subtype in ("deletion", "insertion", "inversion"):
                sizes = element.get(f"flanking_{subtype}_size_ranges")
                frequencies = element.get(f"flanking_{subtype}_frequency_ranges")
                if isinstance(sizes, list) and isinstance(frequencies, list) and len(sizes) != len(frequencies):
                    errors.append(
                        f"is_elements[{index}].flanking_{subtype}_frequency_ranges: "
                        f"expected {len(sizes)} entries to match flanking_{subtype}_size_ranges, "
                        f"got {len(frequencies)}"
                    )
    if errors:
        raise ValueError("Incomplete or invalid config:\n  " + "\n  ".join(errors))


@contextlib.contextmanager
def open_sequence_text(path, mode="r"):
    """Read plain/gzip sequence files and write deterministic gzip streams."""
    if Path(path).suffix.lower() != ".gz":
        with open(path, mode, encoding="utf-8") as handle:
            yield handle
    elif mode == "r":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield handle
    elif mode == "w":
        # Omit filename and timestamp so identical runs have identical checksums.
        with open(path, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8") as handle:
                    yield handle
    else:
        raise ValueError(f"Unsupported sequence file mode: {mode}")


def parse_fasta(path):
    records, name, desc, chunks = [], None, None, []
    with open_sequence_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(Record(name, desc, "".join(chunks).upper()))
                desc = line[1:].strip()
                name = desc.split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append(Record(name, desc, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    return records


def write_fasta(records, path, width=80):
    with open_sequence_text(path, "w") as handle:
        for record in records:
            handle.write(f">{record.desc}\n")
            handle.write(textwrap.fill(record.seq, width))
            handle.write("\n")


def parse_bed(path):
    intervals = []
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"BED line {line_no} has fewer than 3 columns")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"BED line {line_no} must satisfy 0 <= start < end")
            intervals.append({
                "chromosome": fields[0],
                "start": start,
                "end": end,
                "name": fields[3] if len(fields) > 3 else f"bed_{line_no}",
                "line": line_no,
            })
    return intervals


def parse_range(value, cast=float):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return cast(value), cast(value)
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return cast(value[0]), cast(value[1])
    text = str(value).replace(":", "-")
    if "-" not in text:
        return cast(text), cast(text)
    lo, hi = text.split("-", 1)
    lo, hi = cast(lo), cast(hi)
    if hi < lo:
        raise ValueError(f"Invalid range: {value}")
    return lo, hi


def parse_batch_start(value):
    """Accept positive integers or decimal strings such as 001."""
    if isinstance(value, str):
        value = value.strip()
        if not value.isascii() or not value.isdecimal():
            raise argparse.ArgumentTypeError("batch_start must be a positive integer, for example 001")
        try:
            value = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError("batch_start is too large") from exc
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise argparse.ArgumentTypeError("batch_start must be a positive integer, for example 001")
    return value


def parse_count_range(value):
    lo, hi = parse_range(value, float)
    if not all(math.isfinite(bound) and bound >= 0 and bound.is_integer() for bound in (lo, hi)):
        raise ValueError(f"Frequency range bounds must be nonnegative integers: {value}")
    if hi < lo:
        raise ValueError(f"Invalid range: {value}")
    return int(lo), int(hi)


def count_from_frequency(value, total, rng):
    if value is None:
        return 0
    if isinstance(value, str) and ("-" in value or ":" in value):
        return rng.randint(*parse_count_range(value))
    value = float(value)
    if value < 0:
        raise ValueError(f"Frequency cannot be negative: {value}")
    if value < 1:
        expected = total * value
        whole = int(expected)
        return whole + (1 if rng.random() < expected - whole else 0)
    return int(round(value))


def parse_frequency_range_arg(value):
    text = str(value).strip()
    if "-" not in text and ":" not in text:
        raise argparse.ArgumentTypeError("SNP frequency must be an integer count range, for example 0-1 or 10-50")
    try:
        parse_count_range(text)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return text


def parse_sv_frequency(value):
    """Normalize an SV scalar frequency or nonnegative integer count range."""
    if isinstance(value, bool):
        raise argparse.ArgumentTypeError("SV frequency must be a nonnegative scalar or integer count range")
    text = str(value).strip()
    try:
        scalar = float(text)
    except (TypeError, ValueError):
        try:
            lo, hi = parse_count_range(text)
        except (TypeError, ValueError) as exc:
            raise argparse.ArgumentTypeError(str(exc)) from exc
        return f"{lo}-{hi}"
    if not math.isfinite(scalar) or scalar < 0:
        raise argparse.ArgumentTypeError("SV frequency must be finite and nonnegative")
    return scalar


def parse_sv_spec(value):
    cross_count = None
    has_cross_count = False
    if isinstance(value, dict):
        subtype = SV_ALIASES.get(str(value["subtype"]).lower())
        size_range = value["size_range"]
        frequency = value["frequency"]
        has_cross_count = "cross_chromosome_count" in value
        cross_count = value.get("cross_chromosome_count")
    else:
        parts = str(value).replace(":", ",").split(",")
        if len(parts) not in (3, 4):
            raise argparse.ArgumentTypeError(
                "SV spec must be subtype,size_range,frequency[,cross_chromosome_count]"
            )
        subtype = SV_ALIASES.get(parts[0].strip().lower())
        size_range, frequency = parts[1].strip(), parts[2].strip()
        if len(parts) == 4:
            has_cross_count = True
            cross_count = parts[3].strip()
    if subtype is None:
        raise argparse.ArgumentTypeError("Unsupported SV subtype")
    lo, hi = parse_range(size_range, int)
    result = {
        "subtype": subtype,
        "size_range": f"{lo}-{hi}",
        "frequency": parse_sv_frequency(frequency),
    }
    if has_cross_count:
        if subtype != "translocation":
            raise argparse.ArgumentTypeError("cross_chromosome_count is only valid for translocation")
        try:
            count = float(cross_count)
            if isinstance(cross_count, bool) or not math.isfinite(count) or count < 0 or not count.is_integer():
                raise ValueError
        except (TypeError, ValueError, OverflowError) as exc:
            raise argparse.ArgumentTypeError("cross_chromosome_count must be a nonnegative integer") from exc
        result["cross_chromosome_count"] = int(count)
    return result


def normalize_ranges(values):
    return [parse_range(item, int) for item in (values or [])]


def validate_dna(label, seq, allowed=BASES):
    bad = sorted(set(seq.upper()) - set(allowed))
    if bad:
        raise ValueError(f"{label} contains non-{allowed} bases: {','.join(bad)}")



class SourceSequence(str):
    """A normal sequence string carrying a source/composition recipe during planning."""

    def __new__(cls, value, recipe):
        instance = super().__new__(cls, value)
        instance.recipe = recipe
        return instance

    def __reduce__(self):
        return type(self), (str(self), self.recipe)

    def __getitem__(self, item):
        value = super().__getitem__(item)
        if isinstance(item, slice):
            return SourceSequence(value, ["slice", self.recipe, item.start, item.stop, item.step])
        return value

    def __add__(self, other):
        if not isinstance(other, str):
            return NotImplemented
        return SourceSequence(str(self) + str(other), ["concat", self.recipe, sequence_recipe(other)])

    def __radd__(self, other):
        if not isinstance(other, str):
            return NotImplemented
        return SourceSequence(str(other) + str(self), ["concat", sequence_recipe(other), self.recipe])

    def upper(self):
        value = super().upper()
        return self if value == self else SourceSequence(value, ["upper", self.recipe])


def sequence_recipe(sequence):
    return sequence.recipe if isinstance(sequence, SourceSequence) else ["literal", str(sequence)]


def configure_sequence_sources(records, candidates, cfg, args, log):
    """Use symbolic source IDs so output tables do not depend on run-specific paths."""
    def describe(identifier, kind, file_key=None, index=None):
        key = "sequence_source." + identifier
        log.record(key + ".kind", kind)
        if file_key is not None:
            log.record(key + ".file", file_key)
        if index is not None:
            log.record(key + ".index", index)

    describe("reference", "fasta", "input.reference")
    for index, record in enumerate(records):
        record.seq = SourceSequence(record.seq, ["source", "reference", index])
    if candidates:
        describe("candidates", "fasta", "input.insertion_candidates")
        for index, sequence in enumerate(candidates):
            candidates[index] = SourceSequence(sequence, ["source", "candidates", index])
    for index, element in enumerate(cfg.get("is_elements") or [], 1):
        value = element.get("is_seq")
        if not value:
            continue
        identifier = f"is.{index}"
        element["_sequence_source"] = identifier
        try:
            is_file = Path(value).is_file()
        except OSError:
            is_file = False
        if is_file:
            describe(identifier, "fasta", f"input.is_element.{index}")
        elif args.config:
            describe(identifier, "config_literal", "input.config", index - 1)
        else:
            describe(identifier, "command_literal")


SEQUENCE_FIELDS = {"ref_sequence", "alt_sequence", "dest_alt_sequence", "apply_sequence", "inversion_prefix", "inversion_suffix"}


def write_novel_sequences(path, events, replay_events, removals, reference_records):
    """Store inserted fragments, including supplied IS/candidate sequences."""
    replay_by_id = {event["variant_id"]: (order, event) for order, event in enumerate(replay_events)}
    removal_by_id = {event["variant_id"]: event for event in removals}
    reference_indices = {record.name: index for index, record in enumerate(reference_records)}
    with open_sequence_text(path, "w") as handle:
        handle.write("variant_id\tnovel_sequences\tresume_data\n")
        for reported in sorted(events, key=lambda event: event["variant_id"]):
            novel, novel_ids = {}, {}

            def external_insertion(recipe):
                operation = recipe[0]
                if operation == "source":
                    return recipe[1] != "reference"
                if operation in ("slice", "rc", "upper"):
                    return external_insertion(recipe[1])
                if operation == "concat":
                    return any(external_insertion(part) for part in recipe[1:])
                return False

            def pack(recipe):
                operation = recipe[0]
                if operation == "literal":
                    sequence = recipe[1]
                    if not sequence:
                        return ["empty"]
                    if sequence not in novel_ids:
                        identifier = str(len(novel))
                        novel_ids[sequence] = identifier
                        novel[identifier] = sequence
                    return ["novel", novel_ids[sequence]]
                if operation == "source":
                    return recipe
                if operation in ("slice", "rc", "upper"):
                    return [operation, pack(recipe[1]), *recipe[2:]]
                if operation == "concat":
                    return [operation, *(pack(part) for part in recipe[1:])]
                raise ValueError(f"Unsupported sequence recipe: {operation}")

            def encode(event):
                # Recompose inversions on replay instead of saving expanded child sequences.
                event = planned_event(event)
                metadata = {}
                for key, value in event.items():
                    if key == "_before_inversion":
                        continue
                    if key in SEQUENCE_FIELDS:
                        if key == "ref_sequence" and value:
                            recipe = ["slice", ["source", "reference", reference_indices[event["chromosome"]]],
                                      event["pos"], event["pos"] + len(value), None]
                        else:
                            recipe = sequence_recipe(value)
                            # A supplied insertion is still new sequence in the
                            # generated genome. Store its actual bases, not
                            # just an external IS/candidate source reference.
                            if external_insertion(recipe):
                                if key == "alt_sequence" and event["variant_subtype"] == "flanking inversion":
                                    for fragment in (event.get("inversion_prefix", ""), event.get("inversion_suffix", "")):
                                        pack(["literal", str(fragment)])
                                else:
                                    pack(["literal", str(value)])
                        metadata[key] = pack(recipe)
                    else:
                        metadata[key] = value
                return metadata

            order, event = replay_by_id[reported["variant_id"]]
            owned = json.loads(reported.get("removed_reference_is") or "[]")
            payload = {
                "order": order, "event": encode(event),
                "removals": [encode(removal_by_id[item["variant_id"]]) for item in owned],
            }
            handle.write("\t".join((
                reported["variant_id"],
                json.dumps(novel, separators=(",", ":")),
                json.dumps(payload, separators=(",", ":")),
            )) + "\n")


def decode_sequence_recipe(recipe, novel, source_loader):
    operation = recipe[0]
    if operation == "empty":
        return ""
    if operation == "novel":
        return novel[recipe[1]]
    if operation == "source":
        if source_loader is None:
            raise ValueError("Source-derived replay requires a source loader")
        return source_loader(recipe[1], recipe[2])
    if operation == "slice":
        return decode_sequence_recipe(recipe[1], novel, source_loader)[slice(*recipe[2:])]
    if operation == "concat":
        return "".join(decode_sequence_recipe(part, novel, source_loader) for part in recipe[1:])
    if operation == "rc":
        return reverse_complement(decode_sequence_recipe(recipe[1], novel, source_loader))
    if operation == "upper":
        return decode_sequence_recipe(recipe[1], novel, source_loader).upper()
    raise ValueError(f"Unsupported sequence recipe: {operation}")


def reverse_complement(seq):
    result = seq.translate(COMPLEMENT)[::-1].upper()
    return SourceSequence(result, ["rc", seq.recipe]) if isinstance(seq, SourceSequence) else result


def random_dna(length, rng):
    return "".join(rng.choice(BASES) for _ in range(length))


def random_excerpt(seqs, length, rng):
    seqs = [seq for seq in seqs if seq]
    long = [seq for seq in seqs if len(seq) >= length]
    if long:
        seq = rng.choice(long)
        start = rng.randint(0, len(seq) - length)
        return seq[start:start + length]
    joined = "".join(seqs)
    return "".join(rng.choice(joined) for _ in range(length)) if joined else random_dna(length, rng)


def build_kmer_model(seq, k):
    seq = "".join(base for base in seq.upper() if base in BASES)
    model = {}
    for i in range(max(0, len(seq) - k)):
        key = seq[i:i + k]
        nxt = seq[i + k]
        model.setdefault(key, []).append(nxt)
    return model, seq


def generate_kmer(seq, length, k, rng):
    model, clean = build_kmer_model(seq, k)
    if not model:
        return random_dna(length, rng)
    key = rng.choice(list(model))
    out = list(key[:length])
    while len(out) < length:
        suffix = "".join(out[-k:])
        if suffix not in model:
            suffix = rng.choice(list(model))
        out.append(rng.choice(model[suffix]))
    return "".join(out[:length])


def generate_insertion(records, record, pos, length, cfg, candidate_seqs, rng):
    sv_cfg = cfg["sv"]
    if length < int(sv_cfg["long_insertion_min_size"]):
        return random_dna(length, rng)
    if candidate_seqs:
        return random_excerpt(candidate_seqs, length, rng)
    window = int(sv_cfg["long_insertion_local_window"])
    k = int(sv_cfg["long_insertion_k"])
    local = record.seq[max(0, pos - window):min(len(record.seq), pos + window)]
    # Match the model's filtering: a transition needs k + 1 usable bases.
    local = "".join(base for base in local.upper() if base in BASES)
    if len(local) > k:
        return generate_kmer(local, length, k, rng)
    return generate_kmer("".join(r.seq for r in records), length, k, rng)



IS_WEIGHT_EVENTS = {
    "simple insertion", "partial insertion", "nested", "consecutive", "nearby",
    "flanking deletion", "flanking insertion", "flanking inversion",
}


def normalize_chromosome_weights(value):
    """Validate the optional chromosome -> selector -> multiplier mapping."""
    if not isinstance(value, dict):
        raise ValueError("chromosome_weight: expected an object")
    result = {}
    for chromosome, selectors in value.items():
        if not isinstance(chromosome, str) or not chromosome.strip() or not isinstance(selectors, dict):
            raise ValueError("chromosome_weight: expected chromosome names mapped to objects")
        result[chromosome] = {}
        for selector, weight in selectors.items():
            if not isinstance(selector, str) or not selector.strip():
                raise ValueError("chromosome_weight: selectors must be nonempty strings")
            if isinstance(weight, bool) or not isinstance(weight, (int, float)):
                raise ValueError("chromosome_weight: weights must be finite nonnegative numbers")
            try:
                weight = float(weight)
            except (OverflowError, ValueError) as exc:
                raise ValueError("chromosome_weight: weights must be finite nonnegative numbers") from exc
            if not math.isfinite(weight) or weight < 0:
                raise ValueError("chromosome_weight: weights must be finite nonnegative numbers")
            result[chromosome][selector] = weight
    return result


def parse_chromosome_weight(value):
    """Parse a repeatable CLI argument containing comma-separated entries."""
    result = {}
    for entry in value.split(","):
        parts = [part.strip() for part in entry.split(":")]
        if len(parts) != 3 or not parts[0] or not parts[1]:
            raise argparse.ArgumentTypeError("Expected chromosome:selector:weight entries")
        chromosome, selector, weight = parts
        try:
            normalized = normalize_chromosome_weights({chromosome: {selector: float(weight)}})
        except (ValueError, OverflowError) as exc:
            raise argparse.ArgumentTypeError(str(exc)) from exc
        result.setdefault(chromosome, {}).update(normalized[chromosome])
    return result


def validate_chromosome_weight_targets(cfg, records):
    names = {record.name for record in records}
    selectors = {"SV", "IS", "SNP"} | set(SV_ALIASES.values()) | IS_WEIGHT_EVENTS
    selectors.update(
        element.get("raw_name", element.get("name"))
        for element in cfg.get("is_elements", []) if element.get("is_seq")
    )
    for chromosome, settings in cfg.get("chromosome_weight", {}).items():
        if chromosome not in names:
            raise ValueError(f"chromosome_weight: unknown chromosome {chromosome!r}")
        for selector in settings:
            if selector not in selectors:
                raise ValueError(f"chromosome_weight: unknown mode, event, or IS type {selector!r}")


def chromosome_weights(cfg, mode, event=None, is_type=None):
    """Specific IS type overrides event class, which overrides mode, per chromosome."""
    result = {}
    for chromosome, settings in cfg.get("chromosome_weight", {}).items():
        result[chromosome] = next(
            (settings[key] for key in (is_type, event, mode) if key in settings), 1.0
        )
    return result


def choose_record(records, rng, chromosome_weights=None):
    factors = [(chromosome_weights or {}).get(record.name, 1.0) for record in records]
    if any(factor != 1 for factor in factors):
        # Normalize before multiplying to avoid overflow from large valid weights.
        scale = max(factors, default=0)
        if scale == 0:
            return None
        lengths = [len(record.seq) for record in records]
        longest = max(lengths, default=0)
        if longest == 0:
            return None
        weights = [
            (factor / scale) * (length / longest)
            for factor, length in zip(factors, lengths)
        ]
        if not any(weights):
            return None
        return rng.choices(records, weights=weights, k=1)[0]
    # Preserve the existing random stream when all effective weights equal one.
    total = sum(len(r.seq) for r in records)
    pick, seen = rng.randrange(total), 0
    for record in records:
        seen += len(record.seq)
        if pick < seen:
            return record
    return records[-1]


def overlaps(intervals, chrom, start, end):
    return any(start < old_end and end > old_start for old_start, old_end in intervals.get(chrom, []))


class EventReservations(UserDict):
    """Track new event footprints separately from incoming hard exclusions."""

    def __init__(self, occupied):
        self.data = occupied
        self.spacing = {}


def placement_overlaps(intervals, chrom, start, end, cfg=None, target_offset=0):
    """Keep nearest occupied bases 1,000 bp apart, except at preferred BED placements."""
    distance = PLANNED_EVENT_MIN_DISTANCE - 1
    if cfg and cfg["preferred_bed_targeting"].get("preferred_is_bed"):
        pref = cfg["preferred_bed_targeting"]
        cache = cfg.setdefault("_spacing_bed_cache", {})
        bed = pref["preferred_is_bed"]
        if bed not in cache:
            cache[bed] = parse_bed(bed)
        if preferred_interval_hits(cache[bed], chrom, start + target_offset, int(pref["preferred_bed_window"])):
            distance = 0
    spaced = intervals.spacing if isinstance(intervals, EventReservations) else intervals
    return overlaps(intervals, chrom, start, end) or overlaps(spaced, chrom, start - distance, end + distance)


def reserve(intervals, chrom, start, end):
    intervals.setdefault(chrom, []).append((start, max(start + 1, end)))
    if isinstance(intervals, EventReservations):
        intervals.spacing.setdefault(chrom, []).append((start, max(start + 1, end)))


def reserve_event_neighborhoods(intervals, events, distance):
    """Reserve reference spans plus flanks around events for distance rules."""
    for event in events:
        start = event["pos"]
        end = start + max(1, len(event["ref_sequence"]))
        reserve(
            intervals,
            event["chromosome"],
            start - distance,
            end + distance,
        )


def intervals_on_reference(intervals, records):
    records_by_name = {record.name: record for record in records}
    valid = []
    for item in intervals:
        record = records_by_name.get(item["chromosome"])
        if record is None:
            continue
        if item["end"] > len(record.seq):
            raise ValueError(
                f"BED line {item['line']} interval "
                f"{item['chromosome']}:{item['start']}-{item['end']} exceeds "
                f"chromosome length {len(record.seq)}"
            )
        valid.append(item)
    return valid


def preferred_positions(records, cfg, rng, span=1, target_offset=0):
    """Return source ranges whose source + target_offset hits a BED window."""
    pref = cfg["preferred_bed_targeting"]
    bed = pref.get("preferred_is_bed")
    if not bed:
        return []
    intervals = intervals_on_reference(parse_bed(bed), records)
    by_name = {r.name: r for r in records}
    window = int(pref["preferred_bed_window"])
    state = cfg.setdefault("_preferred_targeting_state", {"covered": set()})
    state["intervals"] = intervals
    covered = state["covered"]
    required = math.ceil(
        len(intervals) * float(pref["preferred_bed_target_proportion"])
    )
    if len(covered) >= required:
        return []
    out = []
    for interval_idx, item in enumerate(intervals):
        record = by_name.get(item["chromosome"])
        if not record:
            continue
        start = max(0, item["start"] - window - target_offset)
        end = min(len(record.seq) - span, item["end"] + window - 1 - target_offset)
        if end >= start:
            out.append((record, start, end, interval_idx))
    rng.shuffle(out)
    out.sort(key=lambda candidate: candidate[3] in covered)
    return out


def event_attempts(count, events, failures, category):
    """Stop a category after three failures; appended event rows signal success."""
    for number in range(count):
        if failures.get(category, 0) >= MAX_CONSECUTIVE_EVENT_FAILURES:
            return
        previous_count = len(events)
        yield number
        failures[category] = 0 if len(events) > previous_count else failures.get(category, 0) + 1


def find_position(records, rng, used, span=1, cfg=None, validator=None, target_offset=0, breakpoints_only=False, chromosome_weights=None, spacing_used=None, attempt_budget=None):
    """Find an available position without reserving it or updating targeting state."""
    budget = [MAX_PLACEMENT_ATTEMPTS] if attempt_budget is None else attempt_budget
    def available(record, pos):
        spans = ((pos, pos + 1), (pos + span - 1, pos + span)) if breakpoints_only else ((pos, pos + span),)
        return (
            not any(overlaps(used, record.name, start, end) for start, end in spans)
            and not any(placement_overlaps(
                used if spacing_used is None else spacing_used,
                record.name, start, end, cfg, pos + target_offset - start,
            ) for start, end in spans)
            and (validator is None or validator(record, pos, span))
        )

    preferred = preferred_positions(records, cfg, rng, span, target_offset) if cfg and cfg["preferred_bed_targeting"].get("preferred_is_bed") else []
    for record, start, end, interval_idx in preferred:
        if budget[0] <= 0:
            break
        candidate_count = end - start + 1
        first_offset = rng.randrange(candidate_count)
        for step in range(min(candidate_count, budget[0])):
            budget[0] -= 1
            pos = start + (first_offset + step) % candidate_count
            if available(record, pos):
                return record, pos, interval_idx
    while budget[0] > 0:
        budget[0] -= 1
        record = choose_record(records, rng, chromosome_weights)
        if record is None:
            return None, None, None
        if len(record.seq) < span:
            continue
        pos = rng.randint(0, len(record.seq) - span)
        if available(record, pos):
            return record, pos, None
    return None, None, None


def preferred_interval_hits(intervals, chromosome, pos, window):
    """Credit each BED row whose expanded interval contains the event start."""
    return {
        idx for idx, item in enumerate(intervals)
        if item["chromosome"] == chromosome
        and item["start"] - window <= pos < item["end"] + window
    }


def credit_preferred_position(cfg, chromosome, pos, records=None):
    """Update coverage without changing reservations or consuming randomness."""
    if cfg is not None and cfg["preferred_bed_targeting"].get("preferred_is_bed"):
        state = cfg.setdefault("_preferred_targeting_state", {"covered": set()})
        # find_position initializes the reference-filtered interval list even
        # when the target has been reached or ordinary placement is used.
        intervals = state.get("intervals")
        if intervals is None:
            if records is None:
                raise ValueError("Preferred targeting must be initialized before committing a position")
            intervals = intervals_on_reference(
                parse_bed(cfg["preferred_bed_targeting"]["preferred_is_bed"]), records,
            )
            state["intervals"] = intervals
        state["covered"].update(preferred_interval_hits(
            intervals, chromosome, pos,
            int(cfg["preferred_bed_targeting"]["preferred_bed_window"]),
        ))


def commit_position(used, record, pos, span, cfg=None, interval_idx=None):
    reserve(used, record.name, pos, pos + span)
    credit_preferred_position(cfg, record.name, pos)


def pick_position(records, rng, used, span=1, cfg=None, validator=None, breakpoints_only=False, chromosome_weights=None, spacing_used=None):
    record, pos, interval_idx = find_position(
        records, rng, used, span, cfg, validator, breakpoints_only=breakpoints_only,
        chromosome_weights=chromosome_weights, spacing_used=spacing_used
    )
    if record is not None:
        if breakpoints_only:
            # Commit only this inversion's breakpoints. The interior is never
            # reserved, so reservations belonging to earlier or later events
            # require no removal or special restoration.
            reserve(used, record.name, pos, pos + 1)
            if span > 1:
                reserve(used, record.name, pos + span - 1, pos + span)
            credit_preferred_position(cfg, record.name, pos)
        else:
            commit_position(used, record, pos, span, cfg, interval_idx)
    return record, pos


def make_event(rep, idx, kind, subtype, record, pos, ref="", alt="", **extra):
    event = {
        "variant_id": f"var_{rep:03d}_{idx:06d}",
        "chromosome": record.name,
        "pos": pos,
        "event_type": kind,
        "variant_subtype": subtype,
        "strand": "",
        "size_bp": max(len(ref), len(alt)),
        "dest_chromosome": "",
        "dest_position_1based": "",
        "dest_start_1based": "",
        "source": "",
        "source_chromosome": "",
        "source_start_1based": "",
        "source_end_1based": "",
        "is_size_bp": "",
        "flanking_sv_size_bp": "",
        "ref_sequence": ref,
        "alt_sequence": alt,
        "dest_alt_sequence": "",
        "parent_variant_id": "",
        "parent_event_id": "",
        "placement_context": "",
        "inserted_sequence_offset": "",
        "apply": True,
    }
    event.update(extra)
    return event


def get_is_seq(element):
    value = element.get("is_seq")
    if not value:
        return None
    path = Path(value)
    try:
        is_file = path.is_file()
    except OSError:
        is_file = False
    seq = parse_fasta(path)[0].seq if is_file else "".join(str(value).split()).upper()
    validate_dna(element.get("name", "IS"), seq, element.get("_sequence_alphabet", BASES))
    if element.get("_sequence_source"):
        return SourceSequence(seq, ["source", element["_sequence_source"], 0])
    return seq


def normalize_is_elements(cfg):
    genome_name = str(cfg.get("name") or Path(cfg["input"]).stem).strip()
    if not genome_name:
        raise ValueError("Genome name cannot be empty")
    raw_names = set()
    output_names = set()
    for i, element in enumerate(cfg.get("is_elements") or [], 1):
        value = element.get("is_seq")
        if not value:
            continue
        path = Path(value)
        explicit_name = element.get("name")
        try:
            is_file = path.is_file()
        except OSError:
            is_file = False
        if cfg.get("_template_loaded"):
            if not isinstance(explicit_name, str) or not explicit_name.strip():
                raise ValueError(f"Config IS element {i} requires a nonempty name")
            raw_name = explicit_name.strip()
            element["_sequence_alphabet"] = BASES + "N-"
            if not is_file:
                text = "".join(str(value).split()).upper()
                if ":" in text:
                    raise ValueError(
                        f"Config IS element {i}: IS:sequence values are not allowed; "
                        "provide name separately and use a FASTA path or A/C/G/T/N/- sequence"
                    )
                validate_dna(raw_name, text, BASES + "N-")
                element["is_seq"] = text
            get_is_seq(element)
        elif is_file:
            raw_name = str(explicit_name or path.stem).strip()
        else:
            text = str(value).strip()
            if ":" in text:
                raw_name, seq = text.split(":", 1)
                raw_name = raw_name.strip()
                seq = "".join(seq.split()).upper()
                if not raw_name:
                    raise ValueError("Literal --is-seq must use IS_name:ATGC format")
                element["is_seq"] = seq
            elif explicit_name:
                raw_name = str(explicit_name).strip()
            else:
                raise ValueError(
                    "Literal --is-seq must use IS_name:ATGC format, for example --is-seq IS:ACGTACGT"
                )
        if not raw_name:
            raw_name = f"IS_{i}"
        if raw_name in raw_names:
            raise ValueError(
                f"Duplicate IS name for genome {genome_name}: {raw_name}"
            )
        output_name = f"{genome_name}_{raw_name}"
        if output_name in output_names:
            raise ValueError(f"Duplicate normalized IS name: {output_name}")
        raw_names.add(raw_name)
        output_names.add(output_name)
        element["raw_name"] = raw_name
        element["name"] = output_name


def oriented_is(seq, rng):
    strand = rng.choice(["+", "-"])
    return (seq if strand == "+" else reverse_complement(seq)), strand


def add_is_end_repeats(seq, rng):
    """Draw independent end lengths with weights 11, 10, ..., 1 for 0–10 bp."""
    left_len, right_len = rng.choices(range(11), weights=range(11, 0, -1), k=2)
    return seq[:left_len] + seq + (seq[-right_len:] if right_len else ""), left_len, right_len


def choose_other_is(elements, current, seqs, mixed_percentage, rng):
    if rng.random() < float(mixed_percentage) and len(elements) > 1:
        others = [e for e in elements if e is not current and seqs.get(e["name"])]
        if others:
            return rng.choice(others), True
    return current, False


def plan_reference_is(records, cfg, rng, rep, start_idx):
    ref_cfg = cfg["reference_is"]
    bed = ref_cfg.get("known_is_bed")
    if not bed:
        return [], [], start_idx
    intervals = parse_bed(bed)
    by_name = {r.name: r for r in records}
    candidates = intervals_on_reference(intervals, records)
    candidates_by_type = {}
    for item in candidates:
        candidates_by_type.setdefault(item["name"], []).append(item)
    removed = []
    removal_frequency = ref_cfg.get("reference_removal_frequency_range", 0.5)
    for type_candidates in candidates_by_type.values():
        removal_count = count_from_frequency(removal_frequency, len(type_candidates), rng)
        rng.shuffle(type_candidates)
        removed.extend(type_candidates[:removal_count])
    removed_ids = {(i["chromosome"], i["start"], i["end"], i["line"]) for i in removed}
    truth, removals, idx = [], [], start_idx
    removal_sites = {}
    removal_failures = {}
    for item in candidates:
        record = by_name[item["chromosome"]]
        if ((item["chromosome"], item["start"], item["end"], item["line"]) in removed_ids
                and removal_failures.get(item["name"], 0) < MAX_CONSECUTIVE_EVENT_FAILURES):
            left = rng.randint(0, min(10, item["start"]))
            right = rng.randint(0, min(10, len(record.seq) - item["end"]))
            start, end = item["start"] - left, item["end"] + right
            if placement_overlaps(removal_sites, record.name, start, end, cfg):
                removal_failures[item["name"]] = removal_failures.get(item["name"], 0) + 1
                truth.append(make_event(
                    rep, idx, "ReferenceIS", "reference existing", record, item["start"],
                    record.seq[item["start"]:item["end"]], record.seq[item["start"]:item["end"]],
                    source=item["name"], apply=False,
                ))
                idx += 1
                continue
            reserve(removal_sites, record.name, start, end)
            removal_failures[item["name"]] = 0
            removals.append(make_event(
                rep, idx, "ReferenceIS", "reference removed", record, start,
                record.seq[start:end], "",
                source=f"{item['name']};left_flank_removed={left};right_flank_removed={right}",
            ))
        else:
            truth.append(make_event(
                rep, idx, "ReferenceIS", "reference existing", record, item["start"],
                record.seq[item["start"]:item["end"]], record.seq[item["start"]:item["end"]],
                source=item["name"], apply=False,
            ))
        idx += 1
    return truth, removals, idx


def inversion_site_available(record, pos, span, events):
    """Reject any shared reference bases with an already planned inversion."""
    return not any(
        event["chromosome"] == record.name
        and event["variant_subtype"] in ("inversion", "flanking inversion")
        and pos < event["pos"] + len(event["ref_sequence"])
        and event["pos"] < pos + span
        for event in events
    )


def sv_site_compatible_with_is(subtype, record, pos, span, is_events):
    if subtype == "inversion" and not inversion_site_available(record, pos, span, is_events):
        return False
    candidate_end = pos + span
    for event in is_events:
        if event["chromosome"] != record.name:
            continue
        event_start = event["pos"]
        event_end = event_start + max(1, len(event["ref_sequence"]))
        # A right-side IS also occupies the boundary beyond the inverted span.
        occupied_end = event_end + int(
            event["variant_subtype"] == "flanking inversion"
            and bool(event.get("inversion_suffix"))
        )
        if candidate_end <= event_start or pos >= occupied_end:
            continue
        if (
            event["variant_subtype"] == "flanking inversion"
            and event_start < pos
            and candidate_end < event_end
        ):
            # Its interior is available, but both inversion breakpoints remain
            # occupied. Keep checking for other IS components inside the span.
            continue
        if subtype == "inversion":
            if pos < event_start and occupied_end < candidate_end:
                continue
            return False
        return False
    return True


def plan_sv(records, cfg, rng, rep, start_idx, used, candidate_seqs, enabled, is_events=()):
    if not enabled:
        return [], start_idx
    if not isinstance(used, EventReservations):
        used = EventReservations(used)
    # IS inversions expose their interior; only their breakpoints and IS
    # insertion boundary constrain later SV placements.
    for event in is_events:
        if event["variant_subtype"] == "flanking inversion":
            start = event["pos"]
            end = start + len(event["ref_sequence"])
            reserve(used, event["chromosome"], start, start + 1)
            reserve(used, event["chromosome"], end - 1, end)
            if event.get("inversion_suffix"):
                reserve(used, event["chromosome"], end, end + 1)
        else:
            reserve_event_neighborhoods(used.spacing, [event], 0)
    events, idx = [], start_idx
    total = sum(len(r.seq) for r in records)
    failures = {}
    for raw_spec in cfg["sv"]["specs"]:
        spec = parse_sv_spec(raw_spec)
        weights = chromosome_weights(cfg, "SV", spec["subtype"])
        lo, hi = parse_range(spec["size_range"], int)
        event_count = count_from_frequency(spec["frequency"], total, rng)
        cross_count = spec.get("cross_chromosome_count")
        cross_target = min(cross_count or 0, event_count) if len({r.name for r in records}) > 1 else 0
        for event_number in event_attempts(event_count, events, failures, spec["subtype"]):
            size = rng.randint(lo, hi)
            span = 1 if spec["subtype"] == "insertion" else size
            validator = lambda record, pos, span, subtype=spec["subtype"]: (
                sv_site_compatible_with_is(subtype, record, pos, span, is_events)
                and (subtype != "inversion" or inversion_site_available(record, pos, span, events))
            )
            if spec["subtype"] == "translocation":
                # Legacy three-field specs leave the chromosome relationship unrestricted.
                cross_chromosome = None if cross_count is None else event_number < cross_target
                attempt_budget = [MAX_PLACEMENT_ATTEMPTS]
                for _attempt in range(MAX_PLACEMENT_ATTEMPTS):
                    if attempt_budget[0] <= 0:
                        break
                    record, pos, source_interval_idx = find_position(
                        records, rng, used, span, cfg, validator=validator, chromosome_weights=weights,
                        attempt_budget=attempt_budget,
                    )
                    if record is None:
                        break
                    tentative_used = {
                        chrom: intervals.copy()
                        for chrom, intervals in used.items()
                    }
                    tentative_used = EventReservations(tentative_used)
                    tentative_used.spacing = {chrom: spans.copy() for chrom, spans in used.spacing.items()}
                    reserve(
                        tentative_used, record.name, pos, pos + span
                    )
                    source_chromosome = record.name

                    def destination_compatible(candidate, candidate_pos, candidate_span):
                        if cross_chromosome is not None:
                            if (candidate.name != source_chromosome) != cross_chromosome:
                                return False
                        return sv_site_compatible_with_is(
                            "translocation", candidate, candidate_pos, candidate_span, is_events
                        )

                    # Keep the full reference list so preferred BED coverage IDs stay stable.
                    destination_weights = weights if cross_chromosome is None else {
                        candidate.name: (
                            weights.get(candidate.name, 1.0)
                            if (candidate.name != source_chromosome) == cross_chromosome else 0.0
                        )
                        for candidate in records
                    }
                    dest_record, dest_pos, dest_interval_idx = find_position(
                        records,
                        rng,
                        tentative_used,
                        1,
                        cfg,
                        validator=destination_compatible,
                        chromosome_weights=destination_weights,
                        attempt_budget=attempt_budget,
                    )
                    if dest_record is None:
                        continue
                    commit_position(
                        used, record, pos, span, cfg, source_interval_idx
                    )
                    commit_position(
                        used,
                        dest_record,
                        dest_pos,
                        1,
                        cfg,
                        dest_interval_idx,
                    )
                    ref = record.seq[pos:pos + size]
                    parent_event_id = f"translocation_{rep:03d}_{idx:06d}"
                    deletion = make_event(
                        rep,
                        idx,
                        "SV",
                        "translocation_del",
                        record,
                        pos,
                        ref,
                        "",
                        dest_chromosome=dest_record.name,
                        dest_position_1based=dest_pos + 1,
                        dest_start_1based=dest_pos + 1,
                        source_chromosome=record.name,
                        source_start_1based=pos + 1,
                        source_end_1based=pos + size,
                        parent_event_id=parent_event_id,
                        placement_context="translocation_source",
                    )
                    insertion = make_event(
                        rep,
                        idx + 1,
                        "SV",
                        "translocation_ins",
                        dest_record,
                        dest_pos,
                        "",
                        ref,
                        source_chromosome=record.name,
                        source_start_1based=pos + 1,
                        source_end_1based=pos + size,
                        parent_variant_id=deletion["variant_id"],
                        parent_event_id=parent_event_id,
                        placement_context="translocation_destination",
                    )
                    events.extend([deletion, insertion])
                    idx += 2
                    break
                continue
            if spec["subtype"] == "duplication":
                def duplication_site_available(record, pos, span):
                    destination = pos + span
                    return (
                        validator(record, pos, span)
                        and not placement_overlaps(used, record.name, destination, destination + 1, cfg)
                        and sv_site_compatible_with_is(
                            "duplication", record, destination, 1, is_events
                        )
                    )

                # Check both sites before reserving either. Preferred targeting
                # and coverage use the destination at source + size.
                record, pos, _ = find_position(
                    records, rng, used, span, cfg, validator=duplication_site_available,
                    target_offset=size, chromosome_weights=weights,
                )
                if record is not None:
                    reserve(used, record.name, pos, pos + span)
                    reserve(used, record.name, pos + span, pos + span + 1)
            elif spec["subtype"] == "inversion":
                record, pos = pick_position(
                    records, rng, used, span, cfg, validator=validator,
                    breakpoints_only=True, chromosome_weights=weights,
                )
            else:
                record, pos = pick_position(
                    records, rng, used, span, cfg, validator=validator,
                    chromosome_weights=weights,
                )
            if record is None:
                continue
            ref = "" if spec["subtype"] == "insertion" else record.seq[pos:pos + size]
            alt, extra = "", {}
            if spec["subtype"] == "insertion":
                alt = generate_insertion(records, record, pos, size, cfg, candidate_seqs, rng)
            elif spec["subtype"] == "deletion":
                alt = ""
            elif spec["subtype"] == "inversion":
                alt = reverse_complement(ref)
            elif spec["subtype"] == "duplication":
                source_pos = pos
                pos += size
                credit_preferred_position(cfg, record.name, pos, records)
                alt = ref
                ref = ""
                extra = {
                    "source_chromosome": record.name,
                    "source_start_1based": source_pos + 1,
                    "source_end_1based": source_pos + size,
                    "placement_context": "duplication_destination",
                }
            events.append(make_event(rep, idx, "SV", spec["subtype"], record, pos, ref, alt, **extra))
            idx += 1
    return events, idx


def plan_is(records, cfg, rng, rep, start_idx, used, candidate_seqs, enabled, reference_events=()):
    if not enabled:
        return [], start_idx
    if not isinstance(used, EventReservations):
        used = EventReservations(used)
    # Retained reference IS copies are hard exclusions, including in preferred
    # BED windows. Their own adjacent events may exempt only the parent span.
    reserve_event_neighborhoods(used.data, reference_events, 0)
    elements = cfg.get("is_elements") or []
    seqs = {e.get("name", f"IS_{i+1}"): get_is_seq(e) for i, e in enumerate(elements)}
    events, idx = [], start_idx
    adjacent_reference_ids = set()
    failures = {}
    total = sum(len(r.seq) for r in records)
    for e in elements:
        name = e.get("name", "IS")
        raw_name = e.get("raw_name", name)
        seq = seqs.get(name)
        if not seq:
            continue
        for _ in event_attempts(count_from_frequency(e.get("simple_insertion_frequency_range"), total, rng),
                                events, failures, (name, "simple insertion")):
            record, pos = pick_position(records, rng, used, 1, cfg,
                chromosome_weights=chromosome_weights(cfg, "IS", "simple insertion", raw_name))
            if record is None:
                continue
            oriented, strand = oriented_is(seq, rng)
            alt, left, right = add_is_end_repeats(oriented, rng)
            events.append(make_event(rep, idx, "IS", "simple insertion", record, pos, "", alt, strand=strand, source=f"{name};left_repeat={left};right_repeat={right}"))
            idx += 1
        prop = parse_range(e.get("partial_insertion_proportion_range"), float)
        for _ in event_attempts(count_from_frequency(e.get("partial_insertion_frequency_range"), total, rng),
                                events, failures, (name, "partial insertion")):
            record, pos = pick_position(records, rng, used, 1, cfg,
                chromosome_weights=chromosome_weights(cfg, "IS", "partial insertion", raw_name))
            if record is None or not prop:
                continue
            length = max(1, min(len(seq), round(len(seq) * rng.uniform(*prop))))
            start = rng.randint(0, len(seq) - length)
            oriented, strand = oriented_is(seq[start:start + length], rng)
            alt, left, right = add_is_end_repeats(oriented, rng)
            events.append(make_event(rep, idx, "IS", "partial insertion", record, pos, "", alt, strand=strand, source=f"{name}:{start + 1}-{start + length};left_repeat={left};right_repeat={right}"))
            idx += 1
        adjacent_percentage = float(e.get("reference_adjacent_insertion_percentage") or 0)
        if not 0 <= adjacent_percentage <= 1:
            raise ValueError("reference_adjacent_insertion_percentage must be between 0 and 1")
        matching_references = [
            event for event in reference_events
            if event.get("source") == raw_name
        ]
        adjacent_expected = adjacent_percentage * len(matching_references)
        adjacent_count = int(adjacent_expected)
        if rng.random() < adjacent_expected - adjacent_count:
            adjacent_count += 1
        by_record = {record.name: record for record in records}
        for _ in event_attempts(adjacent_count, events, failures, (name, "reference adjacent")):
            available = []
            for candidate in matching_references:
                if candidate["variant_id"] in adjacent_reference_ids:
                    continue
                candidate_record = by_record[candidate["chromosome"]]
                adjacent_used = EventReservations({
                    chrom: spans.copy() for chrom, spans in used.items()
                })
                adjacent_used.spacing = used.spacing
                parent_span = (candidate["pos"], candidate["pos"] + len(candidate["ref_sequence"]))
                adjacent_used[candidate_record.name].remove(parent_span)
                for side, candidate_pos in (
                    ("left", candidate["pos"]),
                    ("right", candidate["pos"] + len(candidate["ref_sequence"])),
                ):
                    if not placement_overlaps(adjacent_used, candidate_record.name, candidate_pos, candidate_pos + 1, cfg):
                        available.append((candidate, candidate_record, side, candidate_pos))
            if not available:
                continue
            reference, record, side, pos = rng.choice(available)
            inserted_element, mixed = choose_other_is(
                elements,
                e,
                seqs,
                e.get("reference_adjacent_mixed_percentage", 0),
                rng,
            )
            inserted_name = inserted_element.get("name", name)
            inserted_seq = seqs.get(inserted_name) or seq
            adjacent_reference_ids.add(reference["variant_id"])
            reference["variant_subtype"] = (
                "reference_existing_adjacent_mix" if mixed
                else "reference_existing_adjacent_same"
            )
            reserve(used, record.name, pos, pos + 1)
            credit_preferred_position(cfg, record.name, pos, records)
            oriented, strand = oriented_is(inserted_seq, rng)
            alt, left, right = add_is_end_repeats(oriented, rng)
            subtype = "reference_adjacent_mix" if mixed else "reference_adjacent_same"
            events.append(make_event(
                rep, idx, "IS", subtype, record, pos, "", alt,
                strand=strand,
                source=f"{inserted_name};reference_is={reference['source']};left_repeat={left};right_repeat={right}",
                parent_variant_id=reference["variant_id"],
                placement_context=f"is_reference_adjacent_{side}",
            ))
            idx += 1
        paired = [
            ("nested", "nested_insertion_frequency_range", "nested_mixed_percentage", None),
            ("consecutive", "consecutive_insertion_frequency_range", "consecutive_mixed_percentage", "0-0"),
            ("nearby", "nearby_insertion_frequency_range", "nearby_mixed_percentage", e.get("nearby_insertion_gap_range", "1-100")),
        ]
        for label, freq_key, mix_key, gap_range in paired:
            for _ in event_attempts(count_from_frequency(e.get(freq_key), total, rng),
                                    events, failures, (name, label)):
                # Reject impossible nesting before reserving a genome site.
                if label == "nested" and len(seq) < 60:
                    continue
                second, mixed = choose_other_is(elements, e, seqs, e.get(mix_key, 0), rng)
                seq2 = seqs.get(second.get("name")) or seq
                gap = None
                if label == "nearby":
                    lo, hi = parse_range(gap_range, int)
                    gap = rng.randint(lo, hi)
                record, pos = pick_position(records, rng, used, (gap + 1) if gap is not None else 1, cfg,
                    chromosome_weights=chromosome_weights(cfg, "IS", label, raw_name))
                if record is None:
                    continue
                o1, strand1 = oriented_is(seq, rng)
                o2, strand2 = oriented_is(seq2, rng)
                subtype = f"{label} {'mixed' if mixed else 'same'}"
                parent_event_id = f"is_pair_{rep:03d}_{idx:06d}"
                if label == "nested":
                    offset = rng.randint(30, len(o1) - 30)
                    outer = make_event(
                        rep, idx, "IS", f"{subtype} outer", record, pos, "", o1,
                        strand=strand1, source=name, placement_context="is_nested_outer",
                        inserted_sequence_offset=offset, parent_event_id=parent_event_id,
                        apply_sequence=o1[:offset] + o2 + o1[offset:],
                    )
                    inner = make_event(
                        rep, idx + 1, "IS", f"{subtype} inner", record, pos, "", o2,
                        strand=strand2, source=second.get("name"), parent_variant_id=outer["variant_id"],
                        parent_event_id=parent_event_id,
                        placement_context="is_nested_inner",
                        inserted_sequence_offset=offset, apply=False,
                    )
                    events.extend([outer, inner])
                    idx += 2
                    continue
                if label == "consecutive":
                    first = make_event(
                        rep, idx, "IS", f"{subtype} first", record, pos, "", o1,
                        strand=strand1, source=name, parent_event_id=parent_event_id,
                        placement_context="is_consecutive_first",
                    )
                    second_event = make_event(
                        rep, idx + 1, "IS", f"{subtype} second", record, pos, "", o2,
                        strand=strand2, source=second.get("name"), parent_variant_id=first["variant_id"],
                        parent_event_id=parent_event_id,
                        placement_context="is_consecutive_second",
                    )
                    events.extend([first, second_event])
                    idx += 2
                    continue
                first = make_event(
                    rep, idx, "IS", f"{subtype} first", record, pos, "", o1,
                    strand=strand1, source=f"{name};reference_gap={gap}",
                    parent_event_id=parent_event_id, placement_context="is_nearby_first",
                )
                second_event = make_event(
                    rep, idx + 1, "IS", f"{subtype} second", record, pos + gap, "", o2,
                    strand=strand2, source=f"{second.get('name')};reference_gap={gap}",
                    parent_variant_id=first["variant_id"], parent_event_id=parent_event_id,
                    placement_context="is_nearby_second",
                )
                credit_preferred_position(cfg, record.name, second_event["pos"], records)
                events.extend([first, second_event])
                idx += 2
        companions = [
            ("flanking deletion", e.get("flanking_deletion_size_ranges"), e.get("flanking_deletion_frequency_ranges")),
            ("flanking insertion", e.get("flanking_insertion_size_ranges"), e.get("flanking_insertion_frequency_ranges")),
            ("flanking inversion", e.get("flanking_inversion_size_ranges"), e.get("flanking_inversion_frequency_ranges")),
        ]
        for subtype, size_ranges, freq_ranges in companions:
            for size_range, freq_range in zip(normalize_ranges(size_ranges), freq_ranges or []):
                for _ in event_attempts(count_from_frequency(freq_range, total, rng),
                                        events, failures, (name, subtype)):
                    size = rng.randint(*size_range)
                    span = 1 if subtype == "flanking insertion" else size
                    side = rng.choice(["left", "right"]) if subtype == "flanking inversion" else None
                    validator = None
                    if subtype == "flanking inversion":
                        validator = lambda record, pos, span: (
                            inversion_site_available(record, pos, span, events)
                            and (side != "right" or not placement_overlaps(
                                used, record.name, pos + span, pos + span + 1, cfg
                            ))
                        )
                    record, pos = pick_position(
                        records, rng, used, span, cfg,
                        validator=validator,
                        breakpoints_only=subtype == "flanking inversion",
                        chromosome_weights=chromosome_weights(cfg, "IS", subtype, raw_name),
                    )
                    if record is None:
                        continue
                    oriented, strand = oriented_is(seq, rng)
                    if subtype == "flanking insertion":
                        extra = generate_insertion(records, record, pos, size, cfg, candidate_seqs, rng)
                        ref, alt = "", oriented + extra
                    else:
                        ref = record.seq[pos:pos + size]
                        alt = oriented if subtype == "flanking deletion" else oriented + reverse_complement(ref)
                    event_extra = {}
                    if subtype == "flanking inversion":
                        if side == "right":
                            reserve(used, record.name, pos + size, pos + size + 1)
                        prefix = oriented if side == "left" else ""
                        suffix = oriented if side == "right" else ""
                        alt = prefix + reverse_complement(ref) + suffix
                        event_extra = {
                            "inversion_prefix": prefix,
                            "inversion_suffix": suffix,
                            "placement_context": f"is_flanking_inversion_{side}",
                        }
                    events.append(make_event(
                        rep, idx, "IS", subtype, record, pos, ref, alt,
                        strand=strand, source=f"{name};companion_size={size}",
                        size_bp=len(oriented) + size,
                        is_size_bp=len(oriented),
                        flanking_sv_size_bp=size,
                        **event_extra,
                    ))
                    idx += 1
    return events, idx


def generated_event_span(event, sequence):
    """Return the reported half-open span within an event's applied sequence."""
    if event["event_type"] != "IS":
        return 0, len(sequence)
    if event.get("inversion_suffix"):
        return len(sequence) - len(event["inversion_suffix"]), len(sequence)
    if event.get("inversion_prefix"):
        return 0, len(event["inversion_prefix"])
    fields = dict(part.split("=", 1) for part in event.get("source", "").split(";") if "=" in part)
    if "left_repeat" in fields:
        left, right = int(fields["left_repeat"]), int(fields.get("right_repeat", 0))
        length = len(event["alt_sequence"])
        copy_length = next(n for n in range(length + 1)
                           if n + min(n, left) + min(n, right) == length)
        start = min(left, copy_length)
        return start, start + copy_length
    if event["variant_subtype"] == "flanking insertion":
        return 0, int(event.get("is_size_bp") or len(sequence) - int(event.get("flanking_sv_size_bp") or 0))
    return 0, len(sequence)


def generated_component_spans(event, sequence, events):
    """Include composed descendants so their coordinates survive reversals."""
    spans = dict(event.get("_generated_spans", {}))
    spans[event["variant_id"]] = generated_event_span(event, sequence)
    for inner in events:
        if (inner.get("placement_context") == "is_nested_inner"
                and inner.get("parent_variant_id") == event["variant_id"]):
            start = int(inner["inserted_sequence_offset"])
            spans[inner["variant_id"]] = (start, start + len(inner["alt_sequence"]))
    return spans


def integrate_inversion_children(is_events, snp_events, sv_events, reference_events=()):
    """Compose events contained in SV or IS flanking inversions, inner first."""
    for event in is_events + snp_events + sv_events:
        # Keep the plan available so application can recompose against a
        # reference shortened by reference-IS removals.
        if "_before_inversion" not in event:
            event["_before_inversion"] = event.copy()
    inversions = [
        event for event in sv_events + is_events
        if event["variant_subtype"] in ("inversion", "flanking inversion")
    ]
    inversion_children = {}

    def original_span(event):
        return event.get(
            "_original_span",
            (event["pos"], event["pos"] + len(event["ref_sequence"])),
        )

    def containing_inversion(chromosome, start, end, excluded=None):
        candidates = [
            parent for parent in inversions
            if parent is not excluded
            and parent["chromosome"] == chromosome
            and original_span(parent)[0] < start
            and end < original_span(parent)[1]
        ]
        return min(candidates, key=lambda parent: original_span(parent)[1] - original_span(parent)[0], default=None)

    def annotate_inversion(event, parent, position, prefix="inversion"):
        left, right = original_span(parent)
        event[f"{prefix}_variant_id"] = parent["variant_id"]
        event[f"{prefix}_chromosome"] = parent["chromosome"]
        event[f"{prefix}_start_1based"] = left + 1
        event[f"{prefix}_end_1based"] = right
        event[f"{prefix}_offset"] = position - left

    for event in is_events + snp_events + sv_events:
        start = event["pos"]
        original_start, original_end = original_span(event)
        inversion_parent = containing_inversion(
            event["chromosome"], original_start, original_end, excluded=event
        )
        if inversion_parent is not None:
            annotate_inversion(event, inversion_parent, original_start)
            # Non-applied components (e.g. nested inner insertions) still
            # belong to the inversion, but are already composed by their pair.
            if not event.get("apply", True):
                continue
            inversion_children.setdefault(inversion_parent["variant_id"], []).append(
                (event, "source")
            )
            if not event.get("parent_variant_id"):
                event["parent_variant_id"] = inversion_parent["variant_id"]
            if not event.get("placement_context"):
                event["placement_context"] = "inversion_interior"
            if event.get("inserted_sequence_offset") in (None, ""):
                event["inserted_sequence_offset"] = start - inversion_parent["pos"]
            event["apply"] = False

    for event in sv_events:
        if not event.get("dest_alt_sequence"):
            continue
        dest_pos = event.get("_original_dest_pos", int(event["dest_position_1based"]) - 1)
        inversion_parent = containing_inversion(
            event["dest_chromosome"], dest_pos, dest_pos
        )
        if inversion_parent is not None:
            annotate_inversion(event, inversion_parent, dest_pos, "destination_inversion")
            inversion_children.setdefault(inversion_parent["variant_id"], []).append(
                (event, "destination")
            )
            event["apply_destination"] = False

    # Retained reference copies are annotations, not sequence replacements.
    # Track their boundaries through the innermost containing inversion.
    reference_children = {}
    for event in reference_events:
        if event.get("apply", True) or not event["ref_sequence"]:
            continue
        start, end = original_span(event)
        candidates = [parent for parent in inversions
                      if parent["chromosome"] == event["chromosome"]
                      and original_span(parent)[0] <= start
                      and end <= original_span(parent)[1]]
        parent = min(candidates, key=lambda item: original_span(item)[1] - original_span(item)[0], default=None)
        if parent is not None:
            reference_children.setdefault(parent["variant_id"], []).append(event)

    inversion_by_id = {event["variant_id"]: event for event in inversions}
    parent_ids = sorted(
        inversion_children.keys() | reference_children.keys(),
        key=lambda parent_id: original_span(inversion_by_id[parent_id])[1] - original_span(inversion_by_id[parent_id])[0],
    )
    for parent_id in parent_ids:
        children = inversion_children.get(parent_id, [])
        parent = inversion_by_id[parent_id]
        parent_start = parent["pos"]
        replacements = {}
        insertions_by_pos = {}
        for child, role in children:
            if role == "destination":
                offset = int(child["dest_position_1based"]) - 1 - parent_start
                insertions_by_pos.setdefault(offset, []).append(
                    (child["dest_alt_sequence"], None)
                )
                continue
            offset = child["pos"] - parent_start
            sequence = child.get("apply_sequence", child["alt_sequence"])
            if child["ref_sequence"]:
                replacements[offset] = (
                    offset + len(child["ref_sequence"]), sequence, child
                )
            else:
                insertions_by_pos.setdefault(offset, []).append((sequence, child))
        source = parent["ref_sequence"]
        pieces, offset = [], 0
        spans, length = {}, 0
        reference_starts, reference_ends = {}, {}
        for retained in reference_children.get(parent_id, []):
            start = retained["pos"] - parent_start
            reference_starts.setdefault(start, []).append(retained["variant_id"])
            reference_ends.setdefault(start + len(retained["ref_sequence"]), []).append(retained["variant_id"])
        retained_starts, retained_ends = {}, {}

        def emit(sequence, child=None):
            nonlocal length
            if child is not None:
                for variant_id, (start, end) in generated_component_spans(child, sequence, is_events).items():
                    spans[variant_id] = (length + start, length + end)
            pieces.append(sequence)
            length += len(sequence)

        while offset < len(source):
            for variant_id in reference_ends.get(offset, []):
                retained_ends[variant_id] = length
            if offset in insertions_by_pos:
                for sequence, child in insertions_by_pos[offset]:
                    emit(sequence, child)
            for variant_id in reference_starts.get(offset, []):
                retained_starts[variant_id] = length
            if offset in replacements:
                end, sequence, child = replacements[offset]
                emit(sequence, child)
                offset = end
            else:
                emit(source[offset])
                offset += 1
        for variant_id in reference_ends.get(len(source), []):
            retained_ends[variant_id] = length
        for variant_id in retained_starts.keys() & retained_ends.keys():
            spans[variant_id] = (retained_starts[variant_id], retained_ends[variant_id])
        if len(source) in insertions_by_pos:
            for sequence, child in insertions_by_pos[len(source)]:
                emit(sequence, child)
        prefix_length = len(parent.get("inversion_prefix", ""))
        parent["_generated_spans"] = {
            variant_id: (prefix_length + length - end, prefix_length + length - start)
            for variant_id, (start, end) in spans.items()
        }
        parent["apply_sequence"] = (
            parent.get("inversion_prefix", "")
            + reverse_complement("".join(pieces))
            + parent.get("inversion_suffix", "")
        )


def plan_snps(records, cfg, rng, rep, start_idx, used, enabled, spacing_used=None):
    if not enabled:
        return [], start_idx
    if not isinstance(used, EventReservations):
        used = EventReservations(used)
    events, idx = [], start_idx
    total = sum(len(r.seq) for r in records)
    snp_cfg = cfg["snp"]
    target = count_from_frequency(snp_cfg["frequency"], total, rng)
    weights = float(snp_cfg["transition"]) + float(snp_cfg["transversion"])
    for _ in event_attempts(target, events, {}, "SNP"):
        record, pos = pick_position(
            records,
            rng,
            used,
            1,
            validator=lambda record, pos, span: record.seq[pos].upper() in BASES,
            chromosome_weights=chromosome_weights(cfg, "SNP"),
            spacing_used=spacing_used,
        )
        if record is None:
            continue
        if spacing_used is not None:
            reserve(spacing_used, record.name, pos, pos + 1)
        ref = record.seq[pos].upper()
        if rng.random() < float(snp_cfg["transition"]) / weights:
            alt, subtype = TRANSITIONS[ref], "transition"
        else:
            alt, subtype = rng.choice([b for b in BASES if b not in (ref, TRANSITIONS[ref])]), "transversion"
        events.append(make_event(rep, idx, "SNP", subtype, record, pos, ref, alt))
        idx += 1
    return events, idx


def variant_boundaries(events):
    boundaries = []
    for event in events:
        if event["event_type"] == "SNP":
            continue
        pos = event["pos"]
        ref_len = len(event["ref_sequence"])
        boundaries.append((event["chromosome"], pos))
        if ref_len > 1:
            boundaries.append((event["chromosome"], pos + ref_len - 1))
        if event.get("dest_chromosome") and event.get("dest_position_1based"):
            boundaries.append((event["dest_chromosome"], int(event["dest_position_1based"]) - 1))
    return boundaries


def plan_boundary_snps(records, cfg, rng, rep, start_idx, used, source_events):
    boundaries = variant_boundaries(source_events)
    if not boundaries:
        return [], start_idx
    record_by_name = {record.name: record for record in records}
    snp_cfg = cfg["snp"]
    target = count_from_frequency(snp_cfg.get("boundary_snp_frequency"), sum(len(r.seq) for r in records), rng)
    window = int(snp_cfg["boundary_snp_window"])
    events, idx = [], start_idx
    weights = float(snp_cfg["transition"]) + float(snp_cfg["transversion"])
    candidates = []
    candidate_sites = set()
    boundary_sites = set(boundaries)
    for chrom, boundary in boundaries:
        record = record_by_name.get(chrom)
        if not record:
            continue
        for pos in range(max(0, boundary - window), min(len(record.seq) - 1, boundary + window) + 1):
            site = (chrom, pos)
            if (
                site not in boundary_sites
                and site not in candidate_sites
                and record.seq[pos].upper() in BASES
                and not overlaps(used, chrom, pos, pos + 1)
            ):
                candidates.append((record, pos))
                candidate_sites.add(site)
    rng.shuffle(candidates)
    for record, pos in candidates[:target]:
        reserve(used, record.name, pos, pos + 1)
        ref = record.seq[pos].upper()
        if rng.random() < float(snp_cfg["transition"]) / weights:
            alt, subtype = TRANSITIONS[ref], "boundary transition"
        else:
            alt, subtype = rng.choice([b for b in BASES if b not in (ref, TRANSITIONS[ref])]), "boundary transversion"
        events.append(make_event(rep, idx, "SNP", subtype, record, pos, ref, alt, source=f"boundary_window={window}"))
        idx += 1
    return events, idx


class RemovalCoordinates:
    """Map original reference boundaries through the union of deletions."""

    def __init__(self, intervals):
        self.intervals = []
        for start, end in sorted(intervals):
            if self.intervals and start <= self.intervals[-1][1]:
                old_start, old_end = self.intervals[-1]
                self.intervals[-1] = (old_start, max(old_end, end))
            else:
                self.intervals.append((start, end))
        self.starts = [start for start, _ in self.intervals]
        self.removed_before = []
        total = 0
        for start, end in self.intervals:
            self.removed_before.append(total)
            total += end - start

    def map(self, pos):
        index = bisect_right(self.starts, pos) - 1
        if index < 0:
            return pos
        start, end = self.intervals[index]
        return pos - self.removed_before[index] - (min(pos, end) - start)

    def remove(self, seq):
        pieces, pos = [], 0
        for start, end in self.intervals:
            pieces.append(seq[pos:start])
            pos = end
        pieces.append(seq[pos:])
        return "".join(pieces)


def remap_after_reference_removals(records, events, removals):
    """Create application copies; keep the original-coordinate truth plan intact."""
    intervals = {}
    for event in removals:
        if event.get("apply", True) and event["ref_sequence"]:
            start = event["pos"]
            intervals.setdefault(event["chromosome"], []).append(
                (start, start + len(event["ref_sequence"]))
            )
    coordinates = {
        record.name: RemovalCoordinates(intervals.get(record.name, []))
        for record in records
    }
    reduced = [
        Record(record.name, record.desc, coordinates[record.name].remove(record.seq))
        for record in records
    ]
    sequences = {record.name: record.seq for record in reduced}

    def surviving_sequence(chrom, start, end):
        mapping = coordinates[chrom]
        return sequences[chrom][mapping.map(start):mapping.map(end)]

    adjusted = []
    for original in events:
        event = original.get("_before_inversion", original).copy()
        chrom, start = event["chromosome"], event["pos"]
        end = start + len(event["ref_sequence"])
        event["_original_span"] = (start, end)
        event["pos"] = coordinates[chrom].map(start)
        if event["ref_sequence"]:
            event["ref_sequence"] = surviving_sequence(chrom, start, end)
            if event["event_type"] == "SNP" and not event["ref_sequence"]:
                event["apply"] = False
        subtype = event["variant_subtype"]
        if subtype in ("inversion", "flanking inversion"):
            event["alt_sequence"] = (
                event.get("inversion_prefix", "")
                + reverse_complement(event["ref_sequence"])
                + event.get("inversion_suffix", "")
            )
            event.pop("apply_sequence", None)
        elif subtype in ("duplication", "translocation_ins"):
            event["alt_sequence"] = surviving_sequence(
                event["source_chromosome"],
                int(event["source_start_1based"]) - 1,
                int(event["source_end_1based"]),
            )
        if event.get("dest_chromosome") and event.get("dest_position_1based"):
            dest_pos = int(event["dest_position_1based"]) - 1
            event["_original_dest_pos"] = dest_pos
            event["dest_position_1based"] = coordinates[event["dest_chromosome"]].map(dest_pos) + 1
            event["dest_start_1based"] = event["dest_position_1based"]
            # Support the older single-record translocation representation too.
            if event.get("dest_alt_sequence"):
                event["dest_alt_sequence"] = surviving_sequence(chrom, start, end)
        adjusted.append(event)

    integrate_inversion_children(
        [event for event in adjusted if event["event_type"] == "IS"],
        [event for event in adjusted if event["event_type"] == "SNP"],
        [event for event in adjusted if event["event_type"] == "SV"],
        [event for event in adjusted if event["event_type"] == "ReferenceIS"],
    )
    return reduced, adjusted



def planned_event(event):
    """Recover the uncomposed component for application and complex-SV reporting."""
    return event.get("_before_inversion", event)


def sv_reference_sites(event):
    """Reference spans or insertion sites used by a single SV component."""
    start, end = event.get(
        "_original_span", (event["pos"], event["pos"] + len(event["ref_sequence"]))
    )
    sites = [(event["chromosome"], start, end, "event")]
    if event["variant_subtype"] == "flanking inversion" and event.get("inversion_suffix"):
        sites.append((event["chromosome"], end, end + 1, "is_insertion"))
    if event["variant_subtype"] == "duplication" and event.get("source_chromosome"):
        sites.append((
            event["source_chromosome"], int(event["source_start_1based"]) - 1,
            int(event["source_end_1based"]), "duplication_source",
        ))
    # Modern translocations have independent source and destination rows.
    if event.get("dest_alt_sequence"):
        sites.append((
            event["dest_chromosome"],
            event.get("_original_dest_pos", int(event["dest_position_1based"]) - 1),
            event.get("_original_dest_pos", int(event["dest_position_1based"]) - 1),
            "translocation_destination",
        ))
    return sites


def inversion_breakpoint_overlaps(events):
    """Describe SV interactions touching an SV inversion's first or last base."""
    components = [planned_event(event) for event in events if event["event_type"] == "SV"]
    hits = []
    for inversion in components:
        if inversion["variant_subtype"] != "inversion":
            continue
        chrom, left, right, _ = sv_reference_sites(inversion)[0]
        if right <= left:
            continue
        for other in components:
            if other["variant_id"] == inversion["variant_id"]:
                continue
            for other_chrom, start, end, role in sv_reference_sites(other):
                if other_chrom != chrom:
                    continue
                if other["variant_subtype"] == "inversion" and (
                    left < start and end < right or start < left and right < end
                ):
                    # Strict nesting already has an established composition path.
                    continue
                touched = [
                    side for side, boundary in (("left", left), ("right", right - 1))
                    if (start <= boundary < end if end > start else boundary == start)
                ]
                if touched:
                    hits.append({
                        "inversion_variant_id": inversion["variant_id"],
                        "overlapping_variant_id": other["variant_id"],
                        "overlapping_subtype": other["variant_subtype"],
                        "site_role": role,
                        "breakpoints": touched,
                    })
    return hits


def annotate_complex_svs(events):
    """Annotate reporting copies without changing component application plans."""
    hits = inversion_breakpoint_overlaps(events)
    if not hits:
        return
    by_id = {event["variant_id"]: event for event in events if event["event_type"] == "SV"}
    roots = {identifier: identifier for identifier in by_id}

    def root(identifier):
        while roots[identifier] != identifier:
            roots[identifier] = roots[roots[identifier]]
            identifier = roots[identifier]
        return identifier

    def join(first, second):
        roots[root(second)] = root(first)

    for hit in hits:
        join(hit["inversion_variant_id"], hit["overlapping_variant_id"])
    # Keep both components of a translocation in the same reporting group.
    pairs = {}
    for identifier, event in by_id.items():
        pair = event.get("parent_event_id")
        if pair:
            if pair in pairs:
                join(pairs[pair], identifier)
            else:
                pairs[pair] = identifier
    active = {root(hit["inversion_variant_id"]) for hit in hits}
    groups = {}
    for identifier in by_id:
        if root(identifier) in active:
            groups.setdefault(root(identifier), []).append(identifier)
    for identifiers in groups.values():
        identifiers.sort()
        components = [planned_event(by_id[identifier]) for identifier in identifiers]
        chromosome_spans = {}
        for component in components:
            for chrom, start, end, _ in sv_reference_sites(component):
                left, right = chromosome_spans.get(chrom, (start, max(start + 1, end)))
                chromosome_spans[chrom] = (min(left, start), max(right, start + 1, end))
        subtypes = sorted({component["variant_subtype"] for component in components})
        other_types = sorted({
            "translocation" if subtype.startswith("translocation") else subtype
            for subtype in subtypes if subtype != "inversion"
        })
        label = "inversion + " + (" + ".join(other_types) if other_types else "inversion")
        details = []
        for component in components:
            chrom, start, end, _ = sv_reference_sites(component)[0]
            details.append({
                "variant_id": component["variant_id"],
                "variant_subtype": component["variant_subtype"],
                "chromosome": chrom,
                "start_1based": start + 1,
                "end_1based": max(start + 1, end),
                "size_bp": component["size_bp"],
                "source_chromosome": component.get("source_chromosome", ""),
                "source_start_1based": component.get("source_start_1based", ""),
                "source_end_1based": component.get("source_end_1based", ""),
                "dest_chromosome": component.get("dest_chromosome", ""),
                "dest_position_1based": component.get("dest_position_1based", ""),
                "parent_event_id": component.get("parent_event_id", ""),
            })
        metadata = json.dumps({
            "components": details,
            "breakpoint_overlaps": [
                hit for hit in hits if hit["inversion_variant_id"] in identifiers
            ],
            "application_order": "non-inversions first; inversions shortest original span first, then variant ID",
        }, separators=(",", ":"))
        for identifier in identifiers:
            event = by_id[identifier]
            event["_complex_truth_span"] = chromosome_spans[event["chromosome"]]
            event.setdefault("original_variant_subtype", planned_event(event)["variant_subtype"])
            event["variant_subtype"] = label + (
                " with reference IS removal" if event.get("reference_removal_count") else ""
            )
            event["complex_event_id"] = "complex_" + identifiers[0]
            event["complex_sv_subtypes"] = json.dumps(subtypes, separators=(",", ":"))
            event["complex_sv_details"] = metadata


def apply_breakpoint_complex_events(records, events):
    """Apply non-inversion edits, then invert between mapped reference boundaries.

    Only boundary coordinates are tracked, avoiding per-base coordinate objects.
    Inversions are applied shortest original span first, breaking ties by variant ID.
    Later boundaries strictly inside an earlier inversion move with that inversion.
    """
    components = [planned_event(event) for event in events]
    inversions = [
        event for event in components
        if event["variant_subtype"] in ("inversion", "flanking inversion")
        and event.get("apply", True)
    ]
    inversion_ids = {event["variant_id"] for event in inversions}
    replacements, insertions = {}, {}
    for event in components:
        if event["variant_id"] in inversion_ids:
            continue
        chrom, pos = event["chromosome"], event["pos"]
        alt = event.get("apply_sequence", event["alt_sequence"])
        if event.get("apply", True):
            if event["ref_sequence"]:
                replacements.setdefault(chrom, {})[pos] = (
                    pos + len(event["ref_sequence"]), alt,
                )
            else:
                insertions.setdefault(chrom, {}).setdefault(pos, []).append(alt)
        if event.get("dest_alt_sequence") and event.get("apply_destination", True):
            insertions.setdefault(event["dest_chromosome"], {}).setdefault(
                int(event["dest_position_1based"]) - 1, []
            ).append(event["dest_alt_sequence"])
    out = []
    for record in records:
        local = [event for event in inversions if event["chromosome"] == record.name]
        boundaries = sorted({
            boundary for event in local
            for boundary in (event["pos"], event["pos"] + len(event["ref_sequence"]))
        })
        mapped = {}
        pieces, pos, length, boundary_index = [], 0, 0, 0
        repl, ins = replacements.get(record.name, {}), insertions.get(record.name, {})

        def emit(sequence):
            nonlocal length
            pieces.append(sequence)
            length += len(sequence)

        while pos < len(record.seq):
            while boundary_index < len(boundaries) and boundaries[boundary_index] == pos:
                mapped[boundaries[boundary_index]] = length
                boundary_index += 1
            for sequence in ins.get(pos, []):
                emit(sequence)
            if pos in repl:
                end, alt = repl[pos]
                for inner in sorted(point for point in ins if pos < point < end):
                    for sequence in ins[inner]:
                        emit(sequence)
                # Boundaries deleted by a replacement collapse to its surviving start.
                while boundary_index < len(boundaries) and boundaries[boundary_index] < end:
                    mapped[boundaries[boundary_index]] = length
                    boundary_index += 1
                emit(alt)
                pos = end
            else:
                emit(record.seq[pos])
                pos += 1
        while boundary_index < len(boundaries):
            mapped[boundaries[boundary_index]] = length
            boundary_index += 1
        for sequence in ins.get(len(record.seq), []):
            emit(sequence)
        sequence = "".join(pieces)

        def inversion_order(event):
            start, end = event.get(
                "_original_span", (event["pos"], event["pos"] + len(event["ref_sequence"]))
            )
            return end - start, event["variant_id"]

        for event in sorted(local, key=inversion_order):
            left, right = sorted((
                mapped[event["pos"]],
                mapped[event["pos"] + len(event["ref_sequence"])],
            ))
            prefix = event.get("inversion_prefix", "")
            suffix = event.get("inversion_suffix", "")
            sequence = sequence[:left] + prefix + reverse_complement(sequence[left:right]) + suffix + sequence[right:]
            for boundary, offset in mapped.items():
                if left < offset < right:
                    offset = left + right - offset
                suffix_shift = len(suffix) if offset >= right else 0
                if offset > left:
                    offset += len(prefix)
                offset += suffix_shift
                mapped[boundary] = offset
        out.append(Record(record.name, record.desc, sequence))
    return out


def apply_events(records, events, removals):
    if removals:
        records, events = remap_after_reference_removals(records, events, removals)
    if inversion_breakpoint_overlaps(events):
        return apply_breakpoint_complex_events(records, events)
    replacements, insertions = {}, {}
    for event in events:
        chrom, pos = event["chromosome"], event["pos"]
        ref = event["ref_sequence"]
        alt = event.get("apply_sequence", event["alt_sequence"])
        if event.get("apply", True):
            if event["event_type"] == "SNP":
                replacements.setdefault(chrom, {})[pos] = (pos + 1, alt)
            elif ref:
                replacements.setdefault(chrom, {})[pos] = (pos + len(ref), alt)
            else:
                insertions.setdefault(chrom, {}).setdefault(pos, []).append(alt)
        if event.get("dest_alt_sequence") and event.get("apply_destination", True):
            dest_pos0 = int(event["dest_position_1based"]) - 1
            insertions.setdefault(event["dest_chromosome"], {}).setdefault(dest_pos0, []).append(event["dest_alt_sequence"])
    out = []
    for record in records:
        pieces, pos = [], 0
        repl = replacements.get(record.name, {})
        ins = insertions.get(record.name, {})
        while pos < len(record.seq):
            if pos in ins:
                pieces.extend(ins[pos])
            if pos in repl:
                end, alt = repl[pos]
                for inner in sorted(p for p in ins if pos < p < end):
                    pieces.extend(ins[inner])
                pieces.append(alt)
                pos = end
            else:
                pieces.append(record.seq[pos])
                pos += 1
        if len(record.seq) in ins:
            pieces.extend(ins[len(record.seq)])
        out.append(Record(record.name, record.desc, "".join(pieces)))
    return out


def annotate_generated_starts(records, events, removals):
    """Add generated starts, locating the IS copy itself for IS events."""
    application_records, application_events = records, events
    if removals:
        application_records, application_events = remap_after_reference_removals(
            records, events, removals
        )
    else:
        # Recompose reporting copies with retained reference annotations too;
        # these are absent from the sequence-editing IS event list.
        application_events = [planned_event(event).copy() for event in events]
        integrate_inversion_children(
            [event for event in application_events if event["event_type"] == "IS"],
            [event for event in application_events if event["event_type"] == "SNP"],
            [event for event in application_events if event["event_type"] == "SV"],
            [event for event in application_events if event["event_type"] == "ReferenceIS"],
        )

    originals = {event["variant_id"]: event for event in events}
    for event in events:
        event.pop("generated_start_1based", None)

    def annotate_component(event, sequence, output_pos):
        for variant_id, (start, end) in generated_component_spans(
            event, sequence, application_events
        ).items():
            originals[variant_id]["generated_start_1based"] = output_pos + start + 1

    destinations = {}
    for event in application_events:
        if event.get("dest_alt_sequence") and event.get("apply_destination", True):
            destinations.setdefault(event["dest_chromosome"], {}).setdefault(
                int(event["dest_position_1based"]) - 1, []
            ).append(event["dest_alt_sequence"])

    for record in application_records:
        local = [event for event in application_events if event["chromosome"] == record.name]
        replacements, insertions, boundaries = {}, {}, {}
        for event in local:
            boundaries.setdefault(event["pos"], []).append(event)
            if not event.get("apply", True):
                continue
            alt = event.get("apply_sequence", event["alt_sequence"])
            if event["ref_sequence"]:
                replacements[event["pos"]] = (
                    event["pos"] + len(event["ref_sequence"]), alt, event
                )
            else:
                insertions.setdefault(event["pos"], []).append((alt, event))

        output_pos = pos = 0
        while pos <= len(record.seq):
            for sequence, event in insertions.get(pos, []):
                annotate_component(event, sequence, output_pos)
                output_pos += len(sequence)
            for sequence in destinations.get(record.name, {}).get(pos, []):
                output_pos += len(sequence)
            replacement = replacements.get(pos)
            if replacement:
                end, alt, event = replacement
                annotate_component(event, alt, output_pos)
                # Starts swallowed by a deletion/replacement map to its output
                # junction, matching the way application collapses those
                # reference boundaries.
                for boundary in (point for point in boundaries if pos <= point < end):
                    for planned in boundaries[boundary]:
                        originals[planned["variant_id"]].setdefault(
                            "generated_start_1based", output_pos + 1
                        )
                output_pos += len(alt)
                pos = end
                continue
            for planned in boundaries.get(pos, []):
                originals[planned["variant_id"]].setdefault(
                    "generated_start_1based", output_pos + 1
                )
            if pos == len(record.seq):
                break
            output_pos += 1
            pos += 1



def combine_reference_removal_records(events, removals):
    """Build reporting copies, assigning each overlapping removal to one SV.

    Overlap uses original-reference half-open spans (or insertion points),
    including duplication/translocation source spans. Application plans remain
    untouched. The smallest matching span owns a removal in nested events.
    """
    reported = [event.copy() for event in events]
    standalone = []
    for event in reported:
        event["sv_event_count"] = int(event["event_type"] == "SV")
        event["reference_removal_count"] = 0
    for removal in removals:
        start = removal["pos"]
        end = start + len(removal["ref_sequence"])
        matches = []
        for event in reported:
            if event["event_type"] != "SV":
                continue
            spans = [(event["chromosome"], event["pos"],
                      event["pos"] + len(event["ref_sequence"]))]
            if event.get("source_chromosome") and event.get("source_end_1based"):
                spans.append((event["source_chromosome"],
                              int(event["source_start_1based"]) - 1,
                              int(event["source_end_1based"])))
            for chrom, left, right in spans:
                if chrom == removal["chromosome"] and (
                    (left < end and start < right) if left != right
                    else start <= left < end
                ):
                    matches.append((right - left, event["variant_id"], event))
        if not matches:
            standalone.append({**removal, "sv_event_count": 0,
                               "reference_removal_count": 1})
            continue
        owner = min(matches, key=lambda match: match[:2])[2]
        if not owner["reference_removal_count"]:
            owner["original_variant_subtype"] = owner["variant_subtype"]
            owner["variant_subtype"] += " with reference IS removal"
            owner["removed_reference_is"] = []
        owner["reference_removal_count"] += 1
        owner["removed_reference_is"].append({
            "variant_id": removal["variant_id"],
            "chromosome": removal["chromosome"],
            "start_1based": start + 1,
            "end_1based": end,
            "source": removal.get("source", ""),
            "ref_sequence": removal["ref_sequence"],
        })
    for event in reported:
        if "removed_reference_is" in event:
            event["removed_reference_is"] = json.dumps(event["removed_reference_is"], separators=(",", ":"))
    return reported, standalone


TRUTH_FIELDS = [
    "variant_id", "chromosome", "position_1based", "start_1based", "end_1based",
    "event_type", "variant_subtype", "strand", "size_bp", "dest_chromosome",
    "dest_position_1based", "dest_start_1based", "source", "parent_variant_id", "parent_event_id",
    "placement_context", "inserted_sequence_offset", "source_chromosome",
    "source_start_1based", "source_end_1based", "is_size_bp",
    "flanking_sv_size_bp",
    "original_variant_subtype", "sv_event_count", "reference_removal_count",
    "removed_reference_is",
    "inversion_variant_id", "inversion_chromosome", "inversion_start_1based",
    "inversion_end_1based", "inversion_offset",
    "destination_inversion_variant_id", "destination_inversion_chromosome",
    "destination_inversion_start_1based", "destination_inversion_end_1based",
    "destination_inversion_offset",
    "complex_event_id", "complex_sv_subtypes", "complex_sv_details",
    "generated_start_1based",
]


def truth_row(event):
    """Serialize one event with the same reporting coordinates as truth TSV."""
    # Complex groups report their pooled span without changing application coordinates.
    start0, end = event.get("_complex_truth_span", (
        event["pos"], event["pos"] + max(1, len(event["ref_sequence"]))
    ))
    if event["variant_subtype"] == "flanking inversion":
        # Report the IS insertion boundary on the original reference.
        start0 = event["pos"]
        if event.get("inversion_suffix"):
            start0 += len(event["ref_sequence"])
        end = start0 + 1
    row = {
        **event,
        "position_1based": start0 + 1,
        "start_1based": start0 + 1,
        "end_1based": end,
    }
    if (event.get("placement_context") in ("is_nested_outer", "is_nested_inner")
            and event.get("inversion_variant_id")
            and event.get("inserted_sequence_offset") not in (None, "")):
        # Sign denotes containment in an inversion; keep application/replay
        # offsets positive so they remain valid sequence-slicing boundaries.
        row["inserted_sequence_offset"] = f'-{abs(int(event["inserted_sequence_offset"]))}'
    return "\t".join(str(row.get(field, "")) for field in TRUTH_FIELDS) + "\n"


def write_truth(path, events, preserve_order=False):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(f"{index}.{field}" for index, field in enumerate(TRUTH_FIELDS, 1)) + "\n")
        ordered = events if preserve_order else sorted(
            events, key=lambda e: (e["chromosome"], e["pos"], e["variant_id"])
        )
        for event in ordered:
            handle.write(truth_row(event))
            if preserve_order:
                # Debug rows become visible one by one in injection order.
                handle.flush()


def write_sequences(path, events, replay_events=None, removals=(), reference_records=None):
    if reference_records is None:
        raise ValueError("v3 sequence output requires reference_records")
    return write_novel_sequences(
        path, events, events if replay_events is None else replay_events, removals, reference_records,
    )


def target_log(outdir, prefix, sample, cfg, events, records):
    bed = cfg["preferred_bed_targeting"].get("preferred_is_bed")
    if not bed:
        return
    intervals = intervals_on_reference(parse_bed(bed), records)
    window = int(cfg["preferred_bed_targeting"]["preferred_bed_window"])
    targeted = set()
    for event in events:
        if event["event_type"] in ("SV", "IS"):
            targeted.update(preferred_interval_hits(
                intervals, event["chromosome"], event["pos"], window,
            ))
    required = math.ceil(len(intervals) * float(cfg["preferred_bed_targeting"]["preferred_bed_target_proportion"]))
    path = outdir / f"{prefix}.preferred_bed_targeting.tsv"
    new = not path.exists()
    with open(path, "a", encoding="utf-8") as handle:
        if new:
            handle.write("sample\trequired_intervals\tobserved_intervals\ttotal_intervals\tstatus\n")
        handle.write(f"{sample}\t{required}\t{len(targeted)}\t{len(intervals)}\t{'ok' if len(targeted) >= required else 'shortfall'}\n")
    if len(targeted) < required:
        print(f"WARNING: preferred BED targeting shortfall for {sample}: {len(targeted)}/{required}", file=sys.stderr)


def resolve_config(args):
    cfg = copy.deepcopy(DEFAULTS)
    cfg["dryrun"] = False
    cfg["batch_start"] = DEFAULT_BATCH_START
    cfg["fasta_compression_threshold"] = DEFAULT_FASTA_COMPRESSION_THRESHOLD
    cfg["chromosome_weight"] = {}
    if args.config:
        with open(args.config, encoding="utf-8") as handle:
            config = json.load(handle)
        validate_config_completeness(config)
        deep_update(cfg, config)
        # A dry run plans but never injects events, so its debug event log is
        # intentionally disabled. The CLI flag may enable config-based runs.
        cfg["debug"] = bool(cfg.get("debug") or args.debug) and not cfg["dryrun"]
        cfg["batch_start"] = parse_batch_start(cfg["batch_start"])
        cfg["chromosome_weight"] = normalize_chromosome_weights(cfg["chromosome_weight"])
        cfg["_template_loaded"] = True
        if cfg["modes"].get("complex_mode"):
            cfg["modes"]["snp_enabled"] = True
            cfg["modes"]["sv_enabled"] = True
        return cfg
    cfg["dryrun"] = args.dryrun
    cfg["debug"] = bool(args.debug) and not cfg["dryrun"]
    if args.fasta_compression_threshold is not None:
        cfg["fasta_compression_threshold"] = args.fasta_compression_threshold
    for settings in args.chromosome_weight or []:
        deep_update(cfg["chromosome_weight"], settings)
    if args.input:
        cfg["input"] = args.input
    if args.out:
        cfg["outdir"] = args.out
    if args.prefix:
        cfg["prefix"] = args.prefix
    if args.batch is not None:
        cfg["batch"] = args.batch
    if args.batch_start is not None:
        cfg["batch_start"] = parse_batch_start(args.batch_start)
    if args.seed is not None:
        cfg["seed"] = args.seed
    if args.complex_mode:
        cfg["modes"]["complex_mode"] = True
    if args.snp is not None:
        cfg["modes"]["snp_enabled"] = True
        cfg["snp"]["frequency"] = args.snp
    if args.transition is not None:
        cfg["snp"]["transition"] = args.transition
    if args.transversion is not None:
        cfg["snp"]["transversion"] = args.transversion
    if args.boundary_snp_frequency is not None:
        cfg["snp"]["boundary_snp_frequency"] = args.boundary_snp_frequency
    if args.boundary_snp_window is not None:
        cfg["snp"]["boundary_snp_window"] = args.boundary_snp_window
    if args.sv is not None:
        cfg["modes"]["sv_enabled"] = True
        explicit_specs = [value for value in args.sv if value is not None]
        if explicit_specs:
            cfg["sv"]["specs"] = [parse_sv_spec(value) for value in explicit_specs]
    for key in ("long_insertion_min_size", "long_insertion_candidates_fasta", "long_insertion_local_window", "long_insertion_k"):
        val = getattr(args, key)
        if val is not None:
            cfg["sv"][key] = val
    if args.known_is_bed:
        cfg["reference_is"]["known_is_bed"] = args.known_is_bed
    if args.reference_removal_frequency_range is not None:
        cfg["reference_is"]["reference_removal_frequency_range"] = args.reference_removal_frequency_range
    if args.preferred_is_bed:
        cfg["preferred_bed_targeting"]["preferred_is_bed"] = args.preferred_is_bed
    if args.preferred_bed_target_proportion is not None:
        cfg["preferred_bed_targeting"]["preferred_bed_target_proportion"] = args.preferred_bed_target_proportion
    if args.preferred_bed_window is not None:
        cfg["preferred_bed_targeting"]["preferred_bed_window"] = args.preferred_bed_window
    if args.is_seq:
        if not cfg.get("is_elements"):
            cfg["is_elements"] = copy.deepcopy(DEFAULTS["is_elements"])
        cfg["is_elements"][0]["is_seq"] = args.is_seq
        cfg["_is_cli_enabled"] = True
    is_cli_fields = [
        "simple_insertion_frequency_range",
        "partial_insertion_frequency_range",
        "partial_insertion_proportion_range",
        "nested_insertion_frequency_range",
        "nested_mixed_percentage",
        "reference_adjacent_insertion_percentage",
        "reference_adjacent_mixed_percentage",
        "consecutive_insertion_frequency_range",
        "consecutive_mixed_percentage",
        "nearby_insertion_frequency_range",
        "nearby_insertion_gap_range",
        "nearby_mixed_percentage",
        "flanking_deletion_size_ranges",
        "flanking_deletion_frequency_ranges",
        "flanking_insertion_size_ranges",
        "flanking_insertion_frequency_ranges",
        "flanking_inversion_size_ranges",
        "flanking_inversion_frequency_ranges",
    ]
    if any(getattr(args, field) is not None for field in is_cli_fields):
        if not cfg.get("is_elements"):
            cfg["is_elements"] = copy.deepcopy(DEFAULTS["is_elements"])
        elements = cfg["is_elements"]
        if len(elements) != 1:
            raise ValueError(
                "IS-specific CLI overrides are ambiguous with multiple IS elements; "
                "edit the corresponding is_elements entries in the config file"
            )
        target = elements[0]
        for field in is_cli_fields:
            value = getattr(args, field)
            if value is not None:
                target[field] = value
        cfg["_is_cli_enabled"] = True
    if cfg["modes"].get("complex_mode"):
        cfg["modes"]["snp_enabled"] = True
        cfg["modes"]["sv_enabled"] = True
    return cfg


def has_is_event_settings(element):
    keys = [
        "simple_insertion_frequency_range",
        "partial_insertion_frequency_range",
        "nested_insertion_frequency_range",
        "reference_adjacent_insertion_percentage",
        "consecutive_insertion_frequency_range",
        "nearby_insertion_frequency_range",
        "flanking_deletion_frequency_ranges",
        "flanking_insertion_frequency_ranges",
        "flanking_inversion_frequency_ranges",
    ]
    return any(element.get(key) for key in keys)


def is_mode_enabled(cfg):
    """Return whether usable IS events are enabled for this configuration."""
    requested = (
        cfg["modes"].get("complex_mode")
        or cfg.get("_is_cli_enabled")
        or cfg["modes"].get("is_enabled")
    )
    return bool(requested and any(
        element.get("is_seq") and has_is_event_settings(element)
        for element in cfg.get("is_elements", [])
    ))


class ConfigArgumentParser(argparse.ArgumentParser):
    def parse_args(self, args=None, namespace=None):
        arguments = list(sys.argv[1:] if args is None else args)
        probe = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        probe.add_argument("--config")
        probe.add_argument("--simulate")
        selected, _ = probe.parse_known_args(arguments)
        if selected.simulate is not None:
            if selected.config is not None:
                self.error("--simulate cannot be combined with --config")
            parsed = super().parse_args(arguments, namespace)
            if arguments not in (["--simulate", selected.simulate], ["--simulate=" + selected.simulate]):
                self.error("--simulate must be used alone; inputs and outputs come from the log")
            return parsed
        if selected.config is not None:
            option_arguments = arguments[:arguments.index("--")] if "--" in arguments else arguments
            if any(arg == "--is-seq" or arg.startswith("--is-seq=") for arg in option_arguments):
                print(
                    "WARNING: --is-seq is ignored when --config is provided; "
                    "using the config's is_elements sequences.",
                    file=sys.stderr,
                )
            # Ignore other options before their type/arity validation, too.
            arguments = ["--config", selected.config]
            if "--debug" in option_arguments:
                arguments.append("--debug")
        return super().parse_args(arguments, namespace)


def build_parser():
    parser = ConfigArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              Single-IS mode:
                python3 random_genome_variator.py -i reference.fasta --is-seq IS.fasta

              Mixed IS and SV mode with default values:
                python3 random_genome_variator.py -i reference.fasta --is-seq IS.fasta --sv

              Complex mode:
                python3 random_genome_variator.py -i reference.fasta --is-seq IS.fasta --known-is-bed known_is.bed --complex-mode

              Multiple-IS config:
                python3 random_genome_variator.py --config multi_is_config.json
            """),
    )
    is_defaults = DEFAULTS["is_elements"][0]
    sv_defaults = DEFAULTS["sv"]
    snp_defaults = DEFAULTS["snp"]
    ref_defaults = DEFAULTS["reference_is"]
    preferred_defaults = DEFAULTS["preferred_bed_targeting"]
    general = parser.add_argument_group("General input/output")
    general.add_argument("-i", "--input", help="MANDATORY unless provided in --config: input reference FASTA. Default: not set.")
    general.add_argument("--config", help="JSON template/config file; all other command-line parameters are ignored. Default: not set.")
    general.add_argument("--simulate", metavar="LOG", help="Resume a completed dry run from its log, truth TSVs, and sequence tables. Use alone; writes the recorded FASTA outputs.")
    general.add_argument("--dryrun", action="store_true", help="Write truth TSVs, gzip sequence tables (*.sequences.txt.gz), and logs; skip genome construction and FASTA files. Ignored with --config; use config dryrun instead. Default: false.")
    general.add_argument("--debug", action="store_true", help="After successful genome construction, write each injected event sequentially to SAMPLE.debug.tsv using the truth TSV structure. No debug event file is produced by a dry run.")
    general.add_argument("--fasta-compression-threshold", type=int, metavar="BASES",
                         help=f"Gzip output FASTA when total reference length is at least BASES (0: always). Default: {DEFAULT_FASTA_COMPRESSION_THRESHOLD}. With --config, use fasta_compression_threshold instead.")
    general.add_argument("--chromosome_weight", "--chromosome-weight", action="append",
                         type=parse_chromosome_weight, metavar="CHROM:SELECTOR:WEIGHT",
                         help="Whole-genome chromosome length multiplier, default 1. Repeat or comma-separate entries; selectors: SV, IS, SNP, SV subtype, IS event class, or raw IS name. Specific settings override mode settings. BED and reference-adjacent placement take priority.")
    general.add_argument("--out", help=f"Output directory. Default: {DEFAULTS['outdir']}.")
    general.add_argument("--prefix", default=None, help=f"Output prefix. Default: {DEFAULTS['prefix']}.")
    general.add_argument("--batch", type=int, default=None, help=f"Number of replicates. Default: {DEFAULTS['batch']}.")
    general.add_argument("--batch-start", type=parse_batch_start, metavar="INDEX", default=None,
                         help=f"First replicate index, padded to at least three digits. Default: {DEFAULT_BATCH_START:03d}. Seeds follow the index so split batches match one larger run. With --config, use batch_start.")
    general.add_argument("--seed", type=int, default=None, help=f"Random seed used to generate simulation seeds for batch replicates or a single simulation. Default: {DEFAULTS['seed']}.")

    complex_group = parser.add_argument_group("Complex mode")
    complex_group.add_argument("--complex-mode", action="store_true", help=f"Activate complex mode. Default: {DEFAULTS['modes']['complex_mode']}.")

    snp_group = parser.add_argument_group("SNP mode")
    snp_group.add_argument("--snp", nargs="?", const=snp_defaults["frequency"], type=parse_frequency_range_arg, metavar="RANGE", help=f"Activate SNP mode; optionally provide a nonnegative integer count range such as 0-1 or 10-50. Default rate when the flag has no range: {snp_defaults['frequency']}.")
    snp_group.add_argument("--transition", type=float, help=f"Relative transition SNP weight. Default: {snp_defaults['transition']}.")
    snp_group.add_argument("--transversion", type=float, help=f"Relative transversion SNP weight. Default: {snp_defaults['transversion']}.")
    snp_group.add_argument("--boundary-snp-frequency", help=f"SV/IS-boundary SNP count/rate. This process runs whenever SV or IS mode is active, independently of genome-wide SNP mode. Default: {snp_defaults['boundary_snp_frequency']}.")
    snp_group.add_argument("--boundary-snp-window", type=int, help=f"Window around SV/IS event boundaries for boundary SNPs. Default: {snp_defaults['boundary_snp_window']}.")

    sv_group = parser.add_argument_group("SV mode")
    sv_group.add_argument("--sv", action="append", nargs="?", const=None, metavar="SPEC", help="Activate SV mode; optionally provide a repeatable subtype,size_range,frequency spec, with an optional fourth cross_chromosome_count for translocation. With no spec, use defaults: deletion,100-500,5; insertion,100-500,5; inversion,100-500,5; duplication,100-500,5; translocation,100-500,5,2.")
    sv_group.add_argument("--long-insertion-min-size", type=int, help=f"Minimum size using long-insertion sequence sourcing. Default: {sv_defaults['long_insertion_min_size']}.")
    sv_group.add_argument("--long-insertion-candidates-fasta", help="FASTA source for candidate insertion sequences. Default: not set.")
    sv_group.add_argument("--long-insertion-local-window", type=int, help=f"Local reference window for k-mer insertion generation. Default: {sv_defaults['long_insertion_local_window']}.")
    sv_group.add_argument("--long-insertion-k", type=int, help=f"K-mer size for generated insertion sequence. Default: {sv_defaults['long_insertion_k']}.")

    is_group = parser.add_argument_group("IS mode")
    is_group.add_argument("--is-seq", help="Activate single-IS insertion mode using an IS FASTA file or a literal sequence as IS_name:ATGC. FASTA input uses the filename prefix as the IS name. Default: not set.")
    is_group.add_argument("--simple-insertion-frequency-range", help=f"Full-length IS insertion count/rate. Default: {is_defaults['simple_insertion_frequency_range']}.")
    is_group.add_argument("--partial-insertion-frequency-range", help=f"Partial IS insertion count/rate. Default: {is_defaults['partial_insertion_frequency_range']}.")
    is_group.add_argument("--partial-insertion-proportion-range", help=f"Partial IS length proportion range. Default: {is_defaults['partial_insertion_proportion_range']}.")
    is_group.add_argument("--nested-insertion-frequency-range", help=f"Nested IS insertion count/rate. Default: {is_defaults['nested_insertion_frequency_range']}.")
    is_group.add_argument("--nested-mixed-percentage", type=float, help=f"Fraction of nested events using different IS types. Default: {is_defaults['nested_mixed_percentage']}.")
    is_group.add_argument("--reference-adjacent-insertion-percentage", type=float, help=f"Fraction of preserved known IS records receiving an adjacent IS insertion. Default: {is_defaults['reference_adjacent_insertion_percentage']}.")
    is_group.add_argument("--reference-adjacent-mixed-percentage", type=float, help=f"Fraction of reference-adjacent events using a different IS type. Default: {is_defaults['reference_adjacent_mixed_percentage']}.")
    is_group.add_argument("--consecutive-insertion-frequency-range", help=f"Consecutive IS insertion count/rate. Default: {is_defaults['consecutive_insertion_frequency_range']}.")
    is_group.add_argument("--consecutive-mixed-percentage", type=float, help=f"Fraction of consecutive events using different IS types. Default: {is_defaults['consecutive_mixed_percentage']}.")
    is_group.add_argument("--nearby-insertion-frequency-range", help=f"Nearby IS insertion count/rate. Default: {is_defaults['nearby_insertion_frequency_range']}.")
    is_group.add_argument("--nearby-insertion-gap-range", help=f"Gap range for nearby IS insertions. Default: {is_defaults['nearby_insertion_gap_range']}.")
    is_group.add_argument("--nearby-mixed-percentage", type=float, help=f"Fraction of nearby events using different IS types. Default: {is_defaults['nearby_mixed_percentage']}.")
    is_group.add_argument("--flanking-deletion-size-ranges", nargs="+", help=f"IS-flanking deletion size ranges. Default: {','.join(is_defaults['flanking_deletion_size_ranges'])}.")
    is_group.add_argument("--flanking-deletion-frequency-ranges", nargs="+", help=f"IS-flanking deletion frequency ranges. Default: {','.join(is_defaults['flanking_deletion_frequency_ranges'])}.")
    is_group.add_argument("--flanking-insertion-size-ranges", nargs="+", help=f"IS-flanking insertion size ranges. Default: {','.join(is_defaults['flanking_insertion_size_ranges'])}.")
    is_group.add_argument("--flanking-insertion-frequency-ranges", nargs="+", help=f"IS-flanking insertion frequency ranges. Default: {','.join(is_defaults['flanking_insertion_frequency_ranges'])}.")
    is_group.add_argument("--flanking-inversion-size-ranges", nargs="+", help=f"IS-flanking inversion size ranges. Default: {','.join(is_defaults['flanking_inversion_size_ranges'])}.")
    is_group.add_argument("--flanking-inversion-frequency-ranges", nargs="+", help=f"IS-flanking inversion frequency ranges. Default: {','.join(is_defaults['flanking_inversion_frequency_ranges'])}.")

    reference_group = parser.add_argument_group("Reference IS mode")
    reference_group.add_argument("--known-is-bed", help="Required for reference IS mode: 0-based half-open reference IS BED. Default: not set.")
    reference_group.add_argument("--reference-removal-frequency-range", help=f"Reference IS removal count/rate, applied independently to each IS type in BED column 4. Default: {ref_defaults['reference_removal_frequency_range']}.")

    preferred_group = parser.add_argument_group("Preferred BED targeting")
    preferred_group.add_argument("--preferred-is-bed", help="Required to activate preferred BED targeting: 0-based half-open preferred placement BED. Default: not set.")
    preferred_group.add_argument("--preferred-bed-target-proportion", type=float, help=f"Minimum preferred interval target proportion per replicate. Default: {preferred_defaults['preferred_bed_target_proportion']}.")
    preferred_group.add_argument("--preferred-bed-window", type=int, help=f"Window around preferred intervals for placement. Default: {preferred_defaults['preferred_bed_window']}.")
    return parser


def run_simulation(log):
    parser = build_parser()
    args = parser.parse_args()
    if args.simulate:
        return resume_simulation(args.simulate, log)
    cfg = resolve_config(args)
    if not cfg["input"]:
        raise SystemExit("Input FASTA is required via --input or config input")
    if not cfg["outdir"]:
        raise SystemExit("Output directory is required via --out or config outdir")
    outdir = Path(cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    log.open(outdir / f"{cfg['prefix']}.log")
    log.record("mode", "dryrun" if cfg["dryrun"] else "simulation")
    threshold = cfg["fasta_compression_threshold"]
    if isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 0:
        raise ValueError("fasta_compression_threshold must be a nonnegative integer")
    log.record("compression.fasta_threshold_bases", threshold)
    log.record("compression.sequences", "gzip")
    log.record("seed.master", int(cfg["seed"]))
    log.record("batch", int(cfg["batch"]))
    batch_start = cfg["batch_start"]
    log.record("batch_start", f"{batch_start:03d}")
    log.file("input.reference", cfg["input"])
    if args.config:
        log.file("input.config", args.config)
    for key, value in (
        ("input.known_is_bed", cfg["reference_is"].get("known_is_bed")),
        ("input.preferred_is_bed", cfg["preferred_bed_targeting"].get("preferred_is_bed")),
        ("input.insertion_candidates", cfg["sv"].get("long_insertion_candidates_fasta")),
    ):
        if value:
            log.file(key, value)
    for index, element in enumerate(cfg.get("is_elements") or [], 1):
        value = element.get("is_seq")
        if value:
            try:
                is_file = Path(value).is_file()
            except OSError:
                is_file = False
            if is_file:
                log.file(f"input.is_element.{index}", value)
            else:
                log.record(f"input.is_element.{index}.kind", "literal sequence; replay reads the config file or logged command")
    if cfg["preferred_bed_targeting"].get("preferred_is_bed"):
        log.file("output.preferred_targeting", outdir / f"{cfg['prefix']}.preferred_bed_targeting.tsv", checksum=False)
    records = parse_fasta(cfg["input"])
    reference_bases = sum(len(record.seq) for record in records)
    fasta_extension = "fasta.gz" if reference_bases >= threshold else "fasta"
    log.record("reference.total_bases", reference_bases)
    log.record("compression.fasta", "gzip" if fasta_extension.endswith(".gz") else "none")
    candidate_seqs = []
    cand_fa = cfg["sv"].get("long_insertion_candidates_fasta")
    if cand_fa:
        candidate_seqs = [r.seq for r in parse_fasta(cand_fa)]
        for i, seq in enumerate(candidate_seqs, 1):
            validate_dna(f"long insertion candidate {i}", seq)
    normalize_is_elements(cfg)
    configure_sequence_sources(records, candidate_seqs, cfg, args, log)
    validate_chromosome_weight_targets(cfg, records)

    complex_mode = bool(cfg["modes"].get("complex_mode"))
    known_is_bed = cfg["reference_is"].get("known_is_bed")
    ref_enabled = bool(known_is_bed)
    warnings = []
    if complex_mode and not known_is_bed:
        warnings.append(
            "Complex mode: reference_is.known_is_bed is missing; "
            "reference IS simulation will be skipped."
        )
    is_enabled = is_mode_enabled(cfg)
    if complex_mode and not is_enabled:
        warnings.append(
            "Complex mode: no usable is_seq is configured; "
            "the IS component will be skipped."
        )
    reference_adjacent_enabled = any(
        e.get("reference_adjacent_insertion_percentage") not in (None, 0, 0.0, "0")
        for e in cfg.get("is_elements", [])
    )
    if is_enabled and reference_adjacent_enabled and not known_is_bed:
        warnings.append(
            "IS mode: reference_is.known_is_bed is missing; no "
            "reference_adjacent_same or reference_adjacent_mix events will be generated."
        )
    sv_enabled = complex_mode or bool(cfg["modes"].get("sv_enabled"))
    snp_enabled = complex_mode or bool(cfg["modes"].get("snp_enabled"))
    boundary_snp_enabled = sv_enabled or is_enabled
    seed_rng = random.Random(int(cfg["seed"]))
    # Keep each index's seed identical whether produced in one run or a later chunk.
    for _ in range(batch_start - 1):
        seed_rng.randrange(2**32)

    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)

    for rep in range(batch_start, batch_start + int(cfg["batch"])):
        replicate_seed = seed_rng.randrange(2**32)
        rng = random.Random(replicate_seed)
        sample = f"{cfg['prefix']}_{rep:03d}" if int(cfg["batch"]) > 1 or batch_start != DEFAULT_BATCH_START else cfg["prefix"]
        key = f"replicate.{rep}"
        log.record(f"{key}.sample", sample)
        log.record(f"{key}.seed", replicate_seed)
        log.record(f"{key}.started", dt.datetime.now().astimezone().isoformat())
        fasta_path = outdir / f"{sample}.{fasta_extension}"
        sequences_path = outdir / f"{sample}.sequences.txt.gz"
        for suffix, extension in (("truth", "truth.tsv"), ("fasta", fasta_extension), ("sequences", "sequences.txt.gz")):
            log.file(f"{key}.{suffix}", outdir / f"{sample}.{extension}", checksum=False)
        cfg["_preferred_targeting_state"] = {"covered": set()}
        idx = 1
        ref_truth, removals, idx = plan_reference_is(records, cfg, rng, rep, idx) if ref_enabled else ([], [], idx)
        is_used, sv_used, snp_used = EventReservations({}), EventReservations({}), {}
        reserve_event_neighborhoods(is_used, removals, 0)
        is_events, idx = plan_is(records, cfg, rng, rep, idx, is_used, candidate_seqs, is_enabled, ref_truth)
        reserve_event_neighborhoods(sv_used, removals, 0)
        reserve_event_neighborhoods(sv_used.data, ref_truth, 0)
        reserve_event_neighborhoods(snp_used, ref_truth, 0)
        sv_events, idx = plan_sv(
            records, cfg, rng, rep, idx, sv_used, candidate_seqs,
            sv_enabled, is_events,
        )
        for ev in removals + is_events + sv_events:
            if ev["ref_sequence"] and ev["event_type"] in ("ReferenceIS", "SV", "IS"):
                start = ev["pos"]
                end = start + len(ev["ref_sequence"])
                if ev["variant_subtype"] in ("inversion", "flanking inversion"):
                    # Boundary SNPs may occur within an inversion, but not
                    # directly on either inversion breakpoint.
                    reserve(snp_used, ev["chromosome"], start, start + 1)
                    if end - start > 1:
                        reserve(snp_used, ev["chromosome"], end - 1, end)
                    if ev.get("inversion_suffix"):
                        reserve(snp_used, ev["chromosome"], end, end + 1)
                else:
                    reserve(snp_used, ev["chromosome"], start, end)
        # Boundary SNPs use only physical occupancy. Ordinary SNP spacing includes
        # all planned reference footprints, including insertion/duplication sites.
        snp_spacing = {}
        for event in removals + is_events + sv_events:
            for chrom, start, end, _ in sv_reference_sites(event):
                reserve(snp_spacing, chrom, start, end)
        boundary_events, idx = plan_boundary_snps(
            records, cfg, rng, rep, idx, snp_used, is_events + sv_events
        ) if boundary_snp_enabled else ([], idx)
        snp_events, idx = plan_snps(records, cfg, rng, rep, idx, snp_used, snp_enabled, snp_spacing)
        integrate_inversion_children(
            is_events, boundary_events + snp_events, sv_events
        )
        truth = ref_truth + is_events + sv_events + boundary_events + snp_events
        replay_events = truth
        if not cfg["dryrun"]:
            mutated = apply_events(records, truth, removals)
        annotate_generated_starts(records, truth, removals)
        truth, removal_records = combine_reference_removal_records(truth, removals)
        annotate_complex_svs(truth)
        if not cfg["dryrun"]:
            write_fasta(mutated, fasta_path)
            if cfg["debug"]:
                # Do not create the injection log until construction and FASTA
                # output have both succeeded. Thus every logged row is backed by
                # the completed genome, rather than merely by a planned event.
                debug_path = outdir / f"{sample}.debug.tsv"
                write_truth(debug_path, truth, preserve_order=True)
                log.file(f"{key}.debug", debug_path)
        write_truth(outdir / f"{sample}.truth.tsv", truth)
        log.file(f"{key}.truth", outdir / f"{sample}.truth.tsv")
        write_sequences(
            sequences_path, truth + removal_records,
            replay_events + removal_records, removals, reference_records=records,
        )
        log.file(f"{key}.sequences", sequences_path)
        if removals:
            removal_path = outdir / f"{sample}.reference_is_removal.tsv"
            log.file(f"{key}.removals", removal_path, checksum=False)
            write_truth(removal_path, removal_records)
            log.file(f"{key}.removals", removal_path)
        target_log(outdir, cfg["prefix"], sample, cfg, truth, records)
        if not cfg["dryrun"]:
            log.file(f"{key}.fasta", fasta_path)
        log.record(f"{key}.event_count", len(truth))
        log.record(f"{key}.removal_count", len(removals))
        log.record(f"{key}.status", "complete")
        if cfg["dryrun"]:
            print(f"Dry run: wrote {sample}.truth.tsv, {sequences_path.name}, and logs (genome construction skipped)")
        else:
            debug_message = f", {sample}.debug.tsv" if cfg["debug"] else ""
            print(f"Wrote {fasta_path.name}, {sample}.truth.tsv, and {sequences_path.name}{debug_message}")


def file_checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RunLog:
    """Timestamped TSV diagnostics and a versioned resume manifest, without config dumps."""

    def __init__(self):
        self.handle = None
        self.stderr = sys.stderr
        self.pending = []
        self.started = time.monotonic()

    def record(self, key, value, level="INFO"):
        row = (dt.datetime.now().astimezone().isoformat(timespec="milliseconds"), level, key, str(value))
        if self.handle is None:
            self.pending.append(row)
        else:
            csv.writer(self.handle, delimiter="\t", lineterminator="\n").writerow(row)
            self.handle.flush()

    def open(self, path):
        self.handle = open(path, "w", encoding="utf-8", newline="")
        self.record("format", "random_genome_variator_resume_v3")
        self.record("output.log", Path(path).resolve())
        self.record("status", "running")
        self.record("working_directory", Path.cwd())
        self.record("command", shlex.join(sys.argv))
        self.record("python", platform.python_version())
        self.record("platform", platform.platform())
        self.file("input.script", __file__)
        for row in self.pending:
            csv.writer(self.handle, delimiter="\t", lineterminator="\n").writerow(row)
        self.pending.clear()
        self.handle.flush()

    def file(self, key, path, checksum=True):
        path = Path(path).resolve()
        self.record(key + ".path", path)
        if checksum:
            self.record(key + ".sha256", file_checksum(path))

    def write(self, text):
        self.stderr.write(text)
        if text.strip():
            level = "WARNING" if "WARNING" in text else "ERROR" if "ERROR" in text else "INFO"
            self.record("diagnostic", text.rstrip(), level)
        return len(text)

    def flush(self):
        self.stderr.flush()
        if self.handle:
            self.handle.flush()

    def close(self):
        if self.handle:
            self.record("elapsed_seconds", f"{time.monotonic() - self.started:.6f}")
            self.handle.close()
            self.handle = None


def read_run_manifest(path):
    manifest = {}
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 4:
                manifest[row[2]] = row[3]
    if manifest.get("format") != "random_genome_variator_resume_v3":
        raise ValueError("--simulate requires a v3 replay manifest; generate a new dry run with this version")
    if manifest.get("mode") != "dryrun" or manifest.get("status") != "complete":
        raise ValueError("--simulate requires a completed dry-run log")
    return manifest


def read_replay_sequences(path, source_loader=None):
    """Load sequence columns and application metadata keyed by truth variant ID."""
    replay = {}
    with open_sequence_text(path) as handle:
        fields = handle.readline().rstrip("\n").split("\t")
        required = {"variant_id", "novel_sequences", "resume_data"}
        if not required.issubset(fields):
            raise ValueError(f"{path}: not a v3 sequence table; generate a new dry run")
        for line_number, line in enumerate(handle, 2):
            values = line.rstrip("\n").split("\t")
            if len(values) != len(fields):
                raise ValueError(f"{path}:{line_number}: malformed sequence row")
            row = dict(zip(fields, values))
            try:
                payload = json.loads(row["resume_data"])
                event = payload["event"]
                identifier = row["variant_id"]
                if event["variant_id"] != identifier or identifier in replay:
                    raise ValueError("inconsistent or duplicate replay data")
                novel = json.loads(row["novel_sequences"])
                for component in [event] + payload["removals"]:
                    for key in SEQUENCE_FIELDS & component.keys():
                        component[key] = decode_sequence_recipe(component[key], novel, source_loader)
                replay[identifier] = payload
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number}: invalid replay data: {exc}") from exc
    return replay


def read_replay_truth(path, replay):
    """Restore exact application plans; displayed pooled coordinates are reporting-only."""
    events, removals, rows = [], [], []
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        fields = [field.split(".", 1)[-1] for field in header]
        for line_number, line in enumerate(handle, 2):
            values = line.rstrip("\n").split("\t")
            if len(values) != len(fields):
                raise ValueError(f"{path}:{line_number}: malformed truth row")
            row = dict(zip(fields, values))
            try:
                payload = replay.pop(row["variant_id"])
                event = payload["event"]
                if event["variant_id"] != row["variant_id"]:
                    raise ValueError("inconsistent replay data")
                events.append((int(payload["order"]), event))
                removals.extend(payload["removals"])
                rows.append({**event, "variant_id": row["variant_id"]})
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number}: invalid replay data: {exc}") from exc
    events.sort(key=lambda item: item[0])
    return [event for _, event in events], removals, rows


def resume_simulation(log_path, log):
    log_path = Path(log_path).resolve()
    log.open(log_path.with_suffix(".simulate.log"))
    log.record("mode", "resume")
    log.file("input.dryrun_log", log_path)
    manifest = read_run_manifest(log_path)

    def verified_file(key):
        try:
            path = Path(manifest[key + ".path"])
            expected = manifest[key + ".sha256"]
        except KeyError as exc:
            raise ValueError(f"Dry-run log is missing {key} path or checksum") from exc
        log.file("input." + key, path)
        if file_checksum(path) != expected:
            raise ValueError(f"Checksum mismatch for {key}: {path}; restore the dry-run input/output")
        return path

    reference = verified_file("input.reference")
    records = parse_fasta(reference)
    script_hash = manifest.get("input.script.sha256")
    if script_hash and script_hash != file_checksum(__file__):
        print("WARNING: simulator script differs from the dry-run version", file=sys.stderr)
    batch = int(manifest.get("batch", "0"))
    if batch < 1:
        raise ValueError("Dry-run log contains no replicates")
    batch_start = parse_batch_start(manifest.get("batch_start", DEFAULT_BATCH_START))
    log.record("batch", batch)
    log.record("batch_start", f"{batch_start:03d}")
    log.record("seed.master", manifest.get("seed.master", "unknown"))
    # Validate all replicates before creating any genome outputs.
    plans = []
    inputs = {reference, log_path}
    source_cache = {"reference": [record.seq for record in records]}

    def load_source(identifier, index):
        if identifier not in source_cache:
            key = "sequence_source." + identifier
            kind = manifest.get(key + ".kind")
            if kind in ("fasta", "config_literal"):
                source_path = verified_file(manifest[key + ".file"])
                inputs.add(source_path)
                if kind == "fasta":
                    sequences = [record.seq for record in parse_fasta(source_path)]
                else:
                    with open(source_path, encoding="utf-8") as handle:
                        config = json.load(handle)
                    value = config["is_elements"][int(manifest[key + ".index"])]["is_seq"]
                    sequences = ["".join(value.split()).upper()]
            elif kind == "command_literal":
                arguments = build_parser().parse_args(shlex.split(manifest["command"])[1:])
                value = arguments.is_seq.split(":", 1)[1]
                sequences = ["".join(value.split()).upper()]
            else:
                raise ValueError(f"Missing or unsupported replay source: {identifier}")
            source_cache[identifier] = sequences
        if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(source_cache[identifier]):
            raise ValueError(f"Invalid record index for replay source {identifier}: {index}")
        return source_cache[identifier][index]

    outputs = set()
    for rep in range(batch_start, batch_start + batch):
        key = f"replicate.{rep}"
        if manifest.get(key + ".status") != "complete":
            raise ValueError(f"Incomplete dry-run replicate: {rep}")
        truth_path = verified_file(key + ".truth")
        inputs.add(truth_path)
        sequences = Path(manifest[key + ".sequences.path"]).resolve()
        sequences = verified_file(key + ".sequences")
        inputs.add(sequences)
        replay = read_replay_sequences(sequences, load_source)
        events, removals, _ = read_replay_truth(truth_path, replay)
        if key + ".removals.path" in manifest:
            removal_path = verified_file(key + ".removals")
            inputs.add(removal_path)
            standalone, extra, _ = read_replay_truth(removal_path, replay)
            removals.extend(standalone)
            removals.extend(extra)
        if replay:
            raise ValueError(f"{sequences}: sequence rows have no matching truth variant")
        integrate_inversion_children(
            [event for event in events if event["event_type"] == "IS"],
            [event for event in events if event["event_type"] == "SNP"],
            [event for event in events if event["event_type"] == "SV"],
        )
        fasta = Path(manifest[key + ".fasta.path"]).resolve()
        new_outputs = {fasta}
        if new_outputs & outputs or fasta == sequences:
            raise ValueError("Dry-run log contains conflicting output paths")
        outputs.update(new_outputs)
        plans.append((key, events, removals, fasta))
    if outputs & inputs:
        raise ValueError("Resume output paths would overwrite dry-run inputs")
    for key, events, removals, fasta in plans:
        log.record(key + ".seed", manifest.get(key + ".seed", "unknown"))
        log.record(key + ".random_draws", 0)
        log.file(key + ".fasta", fasta, checksum=False)
        mutated = apply_events(records, events, removals)
        fasta.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(mutated, fasta)
        log.file(key + ".fasta", fasta)
        log.record(key + ".status", "complete")
        print(f"Resumed {key}: wrote {fasta}")


def main():
    log = RunLog()
    with contextlib.redirect_stderr(log):
        try:
            run_simulation(log)
        except BaseException as exc:
            if not isinstance(exc, SystemExit) or exc.code not in (None, 0):
                log.record("exception", f"{type(exc).__name__}: {exc}", "ERROR")
                log.record("status", "failed")
            raise
        else:
            log.record("status", "complete")
        finally:
            log.close()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
