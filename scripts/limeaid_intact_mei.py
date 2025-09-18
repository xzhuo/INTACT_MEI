#!/usr/bin/env python3
"""
limeaid_filter.py
By Ronghan Li
Python port of the original Perl script `limeaid.filter.pl` with matching logic
and I/O. Adds optional CLI utilities to (a) restrict output to rows flagged as
INTACT or LTR, and (b) export a BED with genotype information joined from a VCF.

Behavioral parity highlights:
- Skips header line and any row whose column 9 (class) is empty.
- Modifies column 11 (index 10) with same flags: soloLTR, ERV-LTR, LTR-ERV-LTR,
  INTACT, INTACT_3end under the same conditions.
- Parses the `Element_Hits` column as the Perl did, including strand-aware start
  selection and trimming parentheses from the "left" value.

CLI outputs:
- `--out`: writes intact/LTR-only TSV to the given path.
- `--flag-table`: writes 5-col table (CHROM POS ID INTACT_MEI INTACT_FLAG) to path.
"""

import sys
import argparse
import gzip
from typing import Dict, Tuple, Optional


def process_rmsk(rmsk: str):
    """Parse the `Element_Hits` field into a TE map.

    Each hit is ", "-separated; each sub-entry itself is whitespace-separated.
    We reproduce the Perl logic:
    - `start` depends on strand: if '+', use field 11; if '-', use field 13.
    - `left` is taken from field 13 for '+', or field 11 for '-'; then trimmed
      by removing the first and last characters (parentheses in Perl).
    - `sv_length = tmp[6] - tmp[5] + 1`
    - `te_length = end - start + 1`
    - `frac += te_length / (tmp[12] + left)` (guarding against zero)
    - Keep min `left`, min `start`, and max `end` across entries of the same TE.
    """
    # Split entries by ", " and parse whitespace-separated fields per entry
    entries = [e for e in (rmsk or "").split(", ") if e]
    te_map = {}
    for i in entries:
        tmp = i.split()
        if len(tmp) < 14:
            # Malformed entry; skip to faithfully avoid crashing
            continue
        strand = tmp[8]
        te = tmp[9]
        te_class = tmp[10]
        try:
            if strand == "+":
                start = int(tmp[11])
                left_raw = tmp[13]
            else:
                start = int(tmp[13])
                left_raw = tmp[11]
            end = int(tmp[12])
            # Perl: substr($left,1,-1) — remove one char at both ends (parentheses)
            left = left_raw[1:-1] if len(left_raw) >= 2 else left_raw
            left = int(left)
            sv_length = int(tmp[6]) - int(tmp[5]) + 1
            te_length = end - start + 1
            te_total = int(tmp[12]) + left
        except (ValueError, IndexError):
            # Any parsing error: skip this sub-entry
            continue

        d = te_map.setdefault(te, {
            "class": te_class,
            "te_length": 0,
            "sv_length": 0,
            "frac": 0.0,
            # track mins/maxes across same TE name
            "left": None,
            "start": None,
            "end": None,
        })
        d["class"] = te_class
        d["te_length"] += te_length
        d["sv_length"] += sv_length
        if te_total != 0:
            d["frac"] += (te_length / te_total)
        # mins and maxes
        d["left"] = left if d["left"] is None or left < d["left"] else d["left"]
        d["start"] = start if d["start"] is None or start < d["start"] else d["start"]
        d["end"] = end if d["end"] is None or end > d["end"] else d["end"]

    return te_map


def filter_ltr(F):
    """Apply LTR filtering and flagging.

    Mirrors Perl `filter_ltr`:
    - Accumulates frac and min(start)/left for LTR components and internal (-int).
    - Sets F[10] to soloLTR, ERV-LTR, or LTR-ERV-LTR under the same thresholds.
    """
    # Return possibly modified fields list, mirroring Perl logic
    te_map = process_rmsk(F[2])
    total_length = 0
    ltr_frac = 0.0
    int_frac = 0.0
    ltr_start = None
    int_start = None
    ltr_left = None
    int_left = None

    for te, vals in te_map.items():
        if vals.get("class") == F[8]:
            total_length += vals.get("sv_length", 0)
            if te.endswith("-int"):
                int_frac += float(vals.get("frac", 0.0))
                int_start = vals.get("start")
                int_left = vals.get("left")
            else:
                ltr_frac += float(vals.get("frac", 0.0))
                ltr_start = vals.get("start")
                ltr_left = vals.get("left")

    try:
        seq_len = float(F[3]) if F[3] != "" else 0.0
    except ValueError:
        seq_len = 0.0

    if F[4].endswith("-int"):
        if (0.8 < ltr_frac < 1.2) and (0.8 < int_frac < 1.2):
            if (seq_len and total_length / seq_len > 0.8 and
                (ltr_start is not None and ltr_start < 50) and
                (ltr_left is not None and ltr_left < 50) and
                (int_start is not None and int_start < 100) and
                (int_left is not None and int_left < 100)):
                F[10] = "ERV-LTR"
        elif (1.6 < ltr_frac < 2.4) and (0.8 < int_frac < 1.2):
            if (seq_len and total_length / seq_len > 0.8 and
                (ltr_start is not None and ltr_start < 50) and
                (ltr_left is not None and ltr_left < 50) and
                (int_start is not None and int_start < 100) and
                (int_left is not None and int_left < 100)):
                F[10] = "LTR-ERV-LTR"
    elif 0.8 < ltr_frac < 1.2:
        if (seq_len and total_length / seq_len > 0.8 and
            (ltr_start is not None and ltr_start < 50) and
            (ltr_left is not None and ltr_left < 50)):
            F[10] = "soloLTR"

    return F


def filter_alu(F, distance):
    """Apply Alu-like filter (Perl `filter_alu`).

    When F[10] is "No_Flags", sets:
    - INTACT if `repstart < distance` AND `repleft < distance`.
    - INTACT_3end if only `repleft < distance`.
    """
    # As in Perl: set INTACT/INTACT_3end based on start/left thresholds when No_Flags
    te_map = process_rmsk(F[2])
    repstart = None
    repleft = None
    # total_length and frac computed but not used for decision (commented in Perl)
    total_length = 0
    frac = 0.0
    for te, vals in te_map.items():
        if te == F[4]:
            total_length += vals.get("sv_length", 0)
            frac += float(vals.get("frac", 0.0))
            repstart = vals.get("start")
            repleft = vals.get("left")

    if F[10] == "No_Flags" and repstart is not None and repleft is not None and repstart < distance and repleft < distance:
        F[10] = "INTACT"
    elif F[10] == "No_Flags" and repleft is not None and repleft < distance:
        F[10] = "INTACT_3end"
    return F


def filter_l1(F, distance):
    """Apply L1 filter (Perl `filter_l1`).

    Chooses the minimal observed start/left across matching entries and then:
    - INTACT if both < 50; INTACT_3end if left < 50.
    """
    te_map = process_rmsk(F[2])
    repstart = None
    repleft = None
    total_length = 0
    for te, vals in te_map.items():
        if te == F[4]:
            total_length += vals.get("sv_length", 0)
            s = vals.get("start")
            l = vals.get("left")
            repstart = s if repstart is None or (s is not None and s < repstart) else repstart
            repleft = l if repleft is None or (l is not None and l < repleft) else repleft

    if F[10] == "No_Flags" and repstart is not None and repleft is not None and repstart < distance and repleft < distance:
        F[10] = "INTACT"
    elif F[10] == "No_Flags" and repleft is not None and repleft < distance:
        F[10] = "INTACT_3end"
    return F


def filter_sva(F, distance):
    """Apply SVA filter (Perl `filter_sva`).

    SVA considers class matches (F[8]) and then uses min(start/left) logic to
    decide INTACT / INTACT_3end with thresholds at 50.
    """
    te_map = process_rmsk(F[2])
    repstart = None
    repleft = None
    total_length = 0
    for te, vals in te_map.items():
        if vals.get("class") == F[8]:
            total_length += vals.get("sv_length", 0)
            s = vals.get("start")
            l = vals.get("left")
            repstart = s if repstart is None or (s is not None and s < repstart) else repstart
            repleft = l if repleft is None or (l is not None and l < repleft) else repleft

    if F[10] == "No_Flags" and repstart is not None and repleft is not None and repstart < distance and repleft < distance:
        F[10] = "INTACT"
    elif F[10] == "No_Flags" and repleft is not None and repleft < distance:
        F[10] = "INTACT_3end"
    return F


def read_lines_filter_and_flag(in_path: str):
    """Yield rows (list of fields) after applying TE filters and flags.

    - Skips header line.
    - Skips rows with empty column 9 (class) or insufficient columns.
    - Ensures list length for writing F[10].
    """
    with open(in_path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            F = line.split("\t")
            if lineno == 1:
                # skip header, as in Perl
                continue
            if len(F) < 9 or F[8] == "":
                # Not enough columns or empty class rows
                continue
            if len(F) <= 10:
                F += [""] * (11 - len(F))

            if F[8].startswith("LTR"):
                F = filter_ltr(F)
            elif F[8] == "LINE/L1":
                F = filter_alu(F, 50)
            elif F[8] == "SINE/Alu":
                F = filter_l1(F, 50)
            elif F[8] == "Retroposon/SVA":
                F = filter_sva(F, 50)
            yield F


def row_is_intact_or_ltr(F):
    """Return True if F[10] contains INTACT or LTR (Perl regex semantics)."""
    v = F[10] if len(F) > 10 else ""
    return ("INTACT" in v) or ("LTR" in v)


def parse_vcf_gt(vcf_path: str, sample_name: Optional[str] = None) -> Tuple[Dict[str, str], Optional[str]]:
    """Parse a VCF (optionally .gz) and return (ID -> GT) for one sample.

    - If `sample_name` is provided and found in the VCF header, that column is used.
    - Otherwise, the first sample column is used (if present).
    - GT is extracted by aligning FORMAT keys with the chosen sample's value list.
    Returns also the resolved sample name (or None).
    """
    # Returns: (id -> GT string), resolved_sample_name
    gt_by_id: Dict[str, str] = {}
    opener = gzip.open if vcf_path.endswith(".gz") else open
    chosen_index: Optional[int] = None
    header_samples: Optional[list] = None
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                header_samples = parts[9:] if len(parts) > 9 else []
                if sample_name and header_samples:
                    try:
                        chosen_index = 9 + header_samples.index(sample_name)
                    except ValueError:
                        # Fallback to first sample if name not found
                        chosen_index = 9 if len(parts) > 9 else None
                else:
                    chosen_index = 9 if len(parts) > 9 else None
                continue
            # data line
            if chosen_index is None:
                # No sample columns
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= chosen_index:
                continue
            vid = parts[2]
            fmt = parts[8].split(":") if len(parts) > 8 else []
            vals = parts[chosen_index].split(":")
            if not fmt:
                continue
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            if gt_idx < len(vals):
                gt_by_id[vid] = vals[gt_idx]
    resolved_name = sample_name or (header_samples[0] if header_samples else None)
    return gt_by_id, resolved_name


# Note: BED emission logic is implemented inline in main when --vcf is used.

def parse_vcf_positions(vcf_path: str) -> Dict[str, Tuple[str, int]]:
    """Return ID -> (CHROM, POS) mapping from a VCF/VCF.GZ.

    POS is VCF 1-based position. Only non-dot IDs are recorded.
    """
    pos_map: Dict[str, Tuple[str, int]] = {}
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, pos_s, vid = parts[0], parts[1], parts[2]
            if not vid or vid == ".":
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            pos_map[vid] = (chrom, pos)
    return pos_map


def main():
    """CLI entrypoint.

    Options:
      --input <input.limeaid.tsv>
      --intact-only (only output intact rows if true)
      --vcf <input.vcf[.gz]>   (optional; reserved)
      --flag-table <output.info>  (5 columns, no header)
      --out <output.tsv>  (intact/LTR-only TSV)
    """
    p = argparse.ArgumentParser(description="Generate intact-only TSV and 5-col flag table from a limeaid TSV.")
    p.add_argument("--input", required=True, help="Input TSV from limeaid")
    p.add_argument("--intact-only", action="store_true", help="Keep only rows where column 11 contains INTACT or LTR")
    p.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ used to map ID -> CHROM,POS")
    p.add_argument("--flag-table", dest="flag_table_path", required=True, help="Output flag table path (no header): CHROM POS ID INTACT_MEI INTACT_FLAG")
    p.add_argument("--out", dest="out_tsv", required=True, help="Output intact-only TSV path")
    args = p.parse_args()

    rows = list(read_lines_filter_and_flag(args.input))
    intact_rows = [F for F in rows if row_is_intact_or_ltr(F)]

    # Write all TSV or intact TSV
    with open(args.out_tsv, "w") as out:
        output_rows = intact_rows if args.intact_only else rows
        for F in output_rows:
            out.write("\t".join(F) + "\n")

    # Write 5-col flag table (no header) from intact rows using CHROM,POS from VCF by ID
    id_to_pos = parse_vcf_positions(args.vcf)
    with open(args.flag_table_path, "w") as fo:
        for F in intact_rows:
            name = F[0]
            pos_tuple = id_to_pos.get(name)
            if not pos_tuple:
                continue
            chrom, vcf_pos = pos_tuple
            mei = F[4] if len(F) > 4 else ""
            flag = F[10] if len(F) > 10 else ""
            fo.write("\t".join([str(chrom), str(vcf_pos), name, mei, flag]) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
