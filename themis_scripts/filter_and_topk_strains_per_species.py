#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PanTax-DBG helper:
Filter strain_group_abundance table and keep top-K strains per species.

Equivalent to:
  awk '{if($3>min_reads)}' |
  grep -v "0.00000000" |
  grep -v "," |
  tail -n +2 |
  filter_topk_strains_per_species.py

But:
- header-aware
- column-name aware
- safer and reusable
"""

import argparse
import csv
import gzip
from collections import defaultdict
from typing import Dict, List, Tuple


# -------------------------
# IO helpers
# -------------------------
def smart_open(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def to_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("-inf")

def run(
    *,
    input: str,
    out: str,
    min_reads: float = 5.0,
    topk: int = 5,
    keep_ties: bool = False,
    delim: str = "\t",
) -> str:
    if topk <= 0:
        raise SystemExit("[Error] topk must be >= 1")

    groups: Dict[str, List[Tuple[float, dict]]] = defaultdict(list)

    with smart_open(input, "rt") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise RuntimeError("[Error] Input file has no header")

        required = {
            "species",
            "assigned_reads",
            "rel_abundance_within_species",
            "strains",
        }
        missing = required.difference(reader.fieldnames)
        if missing:
            raise RuntimeError(f"[Error] Missing required columns: {sorted(missing)}")

        for row in reader:
            if to_float(row["assigned_reads"]) <= float(min_reads):
                continue
            if to_float(row["rel_abundance_within_species"]) <= 0:
                continue

            strains = (row["strains"] or "").strip()
            if not strains or "," in strains:
                continue

            species = row["species"]
            score = to_float(row["rel_abundance_within_species"])
            groups[species].append((score, row))

    # if not groups:
    #     raise SystemExit("[Error] No rows left after filtering")
    #20260518add
    if not groups:
        with open(out, "wt", newline="") as fo:
            fo.write("species\tgroup_id\tassigned_reads\trel_abundance_within_species\tstrains\n")
        print("[Warning] No rows left after filtering; wrote empty strain table.")
        return out    

    with open(out, "wt", newline="") as fo:
        writer = None

        for sp in sorted(groups.keys()):
            rows = groups[sp]
            rows.sort(key=lambda x: x[0], reverse=True)

            if len(rows) <= topk:
                selected = rows
            else:
                kth_score = rows[topk - 1][0]
                if keep_ties:
                    selected = [(s, r) for s, r in rows if s >= kth_score]
                else:
                    selected = rows[:topk]

            for _, row in selected:
                if writer is None:
                    writer = csv.DictWriter(
                        fo,
                        fieldnames=row.keys(),
                        delimiter=delim,
                    )
                    writer.writeheader()
                writer.writerow(row)

    return out

# -------------------------
# Main logic
# -------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Filter strain groups and keep top-K strains per species.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True,
                    help="strain_group_abundance TSV (or .gz)")
    ap.add_argument("-o", "--out", required=True,
                    help="Output TSV")
    ap.add_argument("--min-reads", type=float, default=5,
                    help="Minimum assigned_reads to keep a row")
    ap.add_argument("--topk", type=int, required=True,
                    help="Keep top-K strains per species")
    ap.add_argument("--keep-ties", action="store_true",
                    help="Keep all rows tied with the K-th value")
    ap.add_argument("--delim", default="\t",
                    help="Input delimiter")
    args = ap.parse_args()

    run(
        input=args.input,
        out=args.out,
        min_reads=args.min_reads,
        topk=args.topk,
        keep_ties=args.keep_ties,
        delim=args.delim,
    )
   

if __name__ == "__main__":
    main()
