#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PanTax-DBG helper (taxid-keyed):
Select top-K strains per species based on ganon predictions,
with singleton abundance filtering.

Inputs:
- ganon species abundance file: first column is species_taxid
- ganon strain abundance file: first column is strain_taxid (NCBI taxid)
- ref_info file: must include columns:
    strain_name, strain_taxid, species_taxid, genome_path

Output:
TSV with columns:
species_taxid, strain_taxid, strain_name, new_strain_taxid, genome_path
"""

import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict


def read_species_set(ganon_species_file: Path):
    """Read species_taxid set from ganon species file (use column 1)."""
    species_set = set()
    with open(ganon_species_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # If header is not a header (data line), treat it as first row
        if header and header[0].strip():
            h0 = header[0].strip().lower()
            if not ("species" in h0 or "taxid" in h0):
                species_set.add(header[0].strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()
            if sid:
                species_set.add(sid)

    return species_set


def read_ref_mapping_taxid(ref_info_file: Path):
    """
    Build mapping keyed by strain_taxid:
      strain_taxid -> (species_taxid, genome_path, strain_name)
    """
    with open(ref_info_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit("[Error] Ref information file is empty or missing a header.")

        cols = {name: i for i, name in enumerate(header)}
        for need in ("strain_name", "strain_taxid", "species_taxid", "genome_path"):
            if need not in cols:
                sys.exit(f"[Error] Ref information file is missing a necessary column: {need}")

        i_name = cols["strain_name"]
        i_taxid = cols["strain_taxid"]
        i_species = cols["species_taxid"]
        i_path = cols["genome_path"]

        mapping = {}
        dup_taxid = 0
        for row in reader:
            if not row or len(row) <= max(i_name, i_taxid, i_species, i_path):
                continue

            strain_name = str(row[i_name]).strip()
            strain_taxid = str(row[i_taxid]).strip()
            species_taxid = str(row[i_species]).strip()
            genome_path = str(row[i_path]).strip()

            if not strain_taxid or not species_taxid or not genome_path:
                continue

            # If duplicates exist, keep the first and count duplicates (report later)
            if strain_taxid in mapping:
                dup_taxid += 1
                continue

            mapping[strain_taxid] = (species_taxid, genome_path, strain_name)

    return mapping, dup_taxid


def read_ganon_strains_with_abund(ganon_strain_file: Path):
    """
    Read (strain_taxid, abundance) pairs from ganon strain file.
    First column is interpreted as strain_taxid (string), second as float abundance.
    """
    items = []
    with open(ganon_strain_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # If header looks like data, accept it
        if header and len(header) >= 2:
            h0 = header[0].strip().lower()
            h1 = header[1].strip().lower()
            header_is_header = ("strain" in h0) or ("taxid" in h0) or ("abund" in h1)
            if not header_is_header:
                st = str(header[0]).strip()
                try:
                    abund = float(header[1])
                except Exception:
                    abund = None
                if st and abund is not None:
                    items.append((st, abund))

        for row in reader:
            if not row or len(row) < 2:
                continue
            st = str(row[0]).strip()
            try:
                abund = float(row[1])
            except Exception:
                continue
            if st:
                items.append((st, abund))

    return items


def run(ref_info,
        ganon_species,
        ganon_strain,
        top_k,
        singleton_min_abund,
        out_tsv):
    if top_k <= 0:
        raise SystemExit("[Error] top_k must be a positive integer.")

    ref_info_file = Path(ref_info)
    ganon_species_file = Path(ganon_species)
    ganon_strain_file = Path(ganon_strain)
    out_tsv_path = Path(out_tsv)

    # 1) species set
    species_set = read_species_set(ganon_species_file)
    if not species_set:
        raise SystemExit("[Error] No species_taxid read from ganon species file.")

    # 2) ref mapping keyed by strain_taxid
    ref_map, dup_taxid = read_ref_mapping_taxid(ref_info_file)
    if not ref_map:
        raise SystemExit("[Error] Ref mapping is empty. Please check the ref_info file.")

    # 3) predicted strain items (strain_taxid, abund)
    predicted_items = read_ganon_strains_with_abund(ganon_strain_file)
    if not predicted_items:
        raise SystemExit("[Error] No (strain_taxid, abundance) was read from the ganon strain file.")

    # 4) group strains per species using ref_map
    per_species_strains = defaultdict(list)  # species_taxid -> list[(strain_taxid, abund, genome_path, strain_name)]
    not_in_ref = 0
    not_in_species = 0

    for strain_taxid, abund in predicted_items:
        # Here strain_taxid is expected to be NCBI taxid (string)
        if strain_taxid not in ref_map:
            not_in_ref += 1
            continue
        species_taxid, genome_path, strain_name = ref_map[strain_taxid]
        if species_taxid not in species_set:
            not_in_species += 1
            continue
        per_species_strains[species_taxid].append((strain_taxid, abund, genome_path, strain_name))

    # 5) singleton filtering
    filtered_species = {}
    thr = float(singleton_min_abund)
    dropped_singletons = 0

    for sp, arr in per_species_strains.items():
        if len(arr) == 1:
            (_, abund, _, _) = arr[0]
            if not (abund > thr):
                dropped_singletons += 1
                continue
        filtered_species[sp] = arr

    if not filtered_species:
        raise SystemExit("[Error] No species were found after filtering. Check threshold or ID consistency (taxid↔ref_info).")

    # 6) select top-K per species
    selected_records = []  # (species_taxid, strain_taxid, strain_name, new_strain_taxid, genome_path)
    for sp, arr in filtered_species.items():
        # sort by abundance desc, then taxid for stability
        arr_sorted = sorted(arr, key=lambda x: (-x[1], x[0]))
        top_arr = arr_sorted[:top_k]
        for idx, (strain_taxid, abund, genome_path, strain_name) in enumerate(top_arr, start=1):
            new_name = f"{sp}_{idx}"
            selected_records.append((sp, strain_taxid, strain_name, new_name, genome_path))

    if not selected_records:
        raise SystemExit("[Error] No records after Top-K selection. Please check top_k and inputs.")

    # 7) write output
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species_taxid", "strain_taxid", "strain_name", "new_strain_taxid", "genome_path"])
        w.writerows(selected_records)

    # summary
    print(f"[Completed] Top-K mapping table written: {out_tsv_path}")
    print(f"[Statistics] species passed: {len(filtered_species)}; final strains: {len(selected_records)}")
    if dropped_singletons:
        print(f"[Hint] Dropped singleton-species with abundance ≤ threshold: {dropped_singletons}", file=sys.stderr)
    if not_in_ref:
        print(f"[Note] {not_in_ref} predicted strain_taxid not found in ref_info; skipped.", file=sys.stderr)
    if not_in_species:
        print(f"[Note] {not_in_species} strains mapped to species not present in ganon species set; skipped.", file=sys.stderr)
    if dup_taxid:
        print(f"[Warning] {dup_taxid} duplicated strain_taxid rows in ref_info were ignored (kept first).", file=sys.stderr)

    return out_tsv_path


def parse_args():
    p = argparse.ArgumentParser(
        description="Select top-K strains per species (taxid-keyed) with singleton abundance filtering."
    )
    p.add_argument("--ref_info", required=True,
                   help="Ref info TSV (must include strain_name, strain_taxid, species_taxid, genome_path).")
    p.add_argument("--ganon_species", required=True,
                   help="ganon_species_abundance.txt (species_taxid in column 1).")
    p.add_argument("--ganon_strain", required=True,
                   help="ganon_strain_abundance.txt (strain_taxid in column 1, abundance in column 2).")
    p.add_argument("--top_k", type=int, required=True,
                   help="Number of Top-K strains per species.")
    p.add_argument("--singleton_min_abund", type=float,
                   default=1e-7,
                   help="Min abundance threshold for singleton-species (strictly greater than).")
    p.add_argument("--out_tsv", default="ganon_species_strain_topk.tsv",
                   help="Output TSV filename.")
    return p.parse_args()


def main():
    args = parse_args()
    run(
        ref_info=args.ref_info,
        ganon_species=args.ganon_species,
        ganon_strain=args.ganon_strain,
        top_k=args.top_k,
        singleton_min_abund=args.singleton_min_abund,
        out_tsv=args.out_tsv,
    )


if __name__ == "__main__":
    main()
