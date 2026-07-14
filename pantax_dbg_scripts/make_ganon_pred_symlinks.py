#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PanTax-DBG helper (taxid-keyed):
Filter predicted strains by species and Ref database,
and output a species–strain mapping table (no symlinks).
"""

import csv
import sys
import argparse
from pathlib import Path
from collections import defaultdict


def resolve_genome_path(path_text: str, ref_info_file: Path) -> str:
    """Return an absolute genome path for downstream ccDBG construction.

    Relative paths in a ref-info table are interpreted from the current
    PanTax-DBG working directory, with a fallback to the ref-info directory.
    This prevents later graph-building steps, which run inside output
    subdirectories, from resolving example genome paths incorrectly.
    """
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return str(path)

    cwd_path = (Path.cwd() / path).resolve()
    if cwd_path.exists():
        return str(cwd_path)

    return str((ref_info_file.parent / path).resolve())


# ----------------------
# Read species set
# ----------------------
def read_species_set(ganon_species_file: Path):
    species_set = set()
    with open(ganon_species_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        if header:
            h0 = str(header[0]).strip().lower()
        else:
            h0 = ""

        # header is real header → skip
        if not (h0 and ("species_taxid" in h0 or "species" in h0)):
            if header and str(header[0]).strip():
                species_set.add(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()
            if sid:
                species_set.add(sid)

    return species_set


# ----------------------
# Read Ref mapping (taxid-keyed)
# ----------------------
def read_ref_mapping_taxid(ref_info_file: Path):
    """
    Build mapping:
      strain_taxid -> (species_taxid, genome_path, strain_name)
    """
    mapping = {}
    dup_taxid = 0

    with open(ref_info_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit("[Error] Ref information file is empty or missing a header.")

        cols = {str(name).strip(): i for i, name in enumerate(header)}
        for need in ("strain_name", "strain_taxid", "species_taxid", "genome_path"):
            if need not in cols:
                sys.exit(f"[Error] Ref information file is missing a required column: {need}")

        i_name = cols["strain_name"]
        i_taxid = cols["strain_taxid"]
        i_species = cols["species_taxid"]
        i_path = cols["genome_path"]

        for row in reader:
            if not row or len(row) <= max(i_name, i_taxid, i_species, i_path):
                continue

            strain_name = str(row[i_name]).strip()
            strain_taxid = str(row[i_taxid]).strip()
            species_taxid = str(row[i_species]).strip()
            genome_path = str(row[i_path]).strip()

            if not strain_taxid or not species_taxid or not genome_path:
                continue

            genome_path = resolve_genome_path(genome_path, ref_info_file)

            if strain_taxid in mapping:
                dup_taxid += 1
                continue

            mapping[strain_taxid] = (species_taxid, genome_path, strain_name)

    return mapping, dup_taxid


# ----------------------
# Read ganon strain list (taxid)
# ----------------------
def read_ganon_strains_taxid(ganon_strain_file: Path):
    strains = []
    with open(ganon_strain_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        if header:
            h0 = str(header[0]).strip().lower()
        else:
            h0 = ""

        # header is data
        if not h0 or ("strain" not in h0 and "taxid" not in h0):
            if header and str(header[0]).strip():
                strains.append(str(header[0]).strip())

        for row in reader:
            if not row:
                continue
            sid = str(row[0]).strip()
            if sid:
                strains.append(sid)

    return strains


# ----------------------
# Main logic
# ----------------------
def run(ref_info,
        ganon_species,
        ganon_strain,
        out_tsv,
        tmp_dir=None):

    ref_info_file = Path(ref_info)
    ganon_species_file = Path(ganon_species)
    ganon_strain_file = Path(ganon_strain)
    out_tsv_path = Path(out_tsv)

    # 1) species set
    species_set = read_species_set(ganon_species_file)
    if not species_set:
        raise SystemExit("[Error] No species_taxid read from ganon species file.")

    # 2) ref mapping (taxid-keyed)
    ref_map, dup_taxid = read_ref_mapping_taxid(ref_info_file)
    if not ref_map:
        raise SystemExit("[Error] Ref mapping is empty.")

    # 3) predicted strain taxids
    predicted_strains = read_ganon_strains_taxid(ganon_strain_file)
    if not predicted_strains:
        raise SystemExit("[Error] No strain_taxid read from ganon strain file.")

    # 4) filter + assign new ids
    per_species_counter = defaultdict(int)
    records = []  # (species_taxid, strain_taxid, strain_name, new_strain_taxid, genome_path)
    not_in_ref = 0

    for strain_taxid in predicted_strains:
        if strain_taxid not in ref_map:
            not_in_ref += 1
            continue

        species_taxid, genome_path, strain_name = ref_map[strain_taxid]
        if species_taxid not in species_set:
            continue

        per_species_counter[species_taxid] += 1
        new_strain_taxid = f"{species_taxid}_{per_species_counter[species_taxid]}"
        records.append((species_taxid, strain_taxid, strain_name, new_strain_taxid, genome_path))

    if not records:
        raise SystemExit("[Error] No records were found after filtering.")

    # 5) write output
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "species_taxid",
            "strain_taxid",
            "strain_name",
            "new_strain_taxid",
            "genome_path",
        ])
        writer.writerows(records)

    if tmp_dir:
        print(
            f"[Note] PanTax-DBG: tmp_dir={tmp_dir} has been ignored; symbolic links are no longer created.",
            file=sys.stderr
        )

    print(f"[Completed] Mapping table written: {out_tsv_path}")
    print(f"[Statistics] species: {len(per_species_counter)}; strains: {len(records)}")
    if not_in_ref:
        print(f"[Note] {not_in_ref} predicted strain_taxid not found in Ref and skipped.", file=sys.stderr)
    if dup_taxid:
        print(f"[Warning] {dup_taxid} duplicated strain_taxid rows in ref_info were ignored.", file=sys.stderr)

    return out_tsv_path


# ----------------------
# CLI
# ----------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Filter predicted strains using Ref database (taxid-keyed, no symlinks)."
    )
    p.add_argument("--ref_info", required=True,
                   help="Ref info TSV (strain_name, strain_taxid, species_taxid, genome_path).")
    p.add_argument("--ganon_species", required=True,
                   help="Species abundance file.")
    p.add_argument("--ganon_strain", required=True,
                   help="Strain abundance file (first column = strain_taxid).")
    p.add_argument("--out_tsv", default="ganon_species_strain_selected.taxid.tsv",
                   help="Output TSV filename.")
    p.add_argument("--tmp_dir", default=None,
                   help="Deprecated; ignored.")
    return p.parse_args()


def main():
    args = parse_args()
    run(
        ref_info=args.ref_info,
        ganon_species=args.ganon_species,
        ganon_strain=args.ganon_strain,
        out_tsv=args.out_tsv,
        tmp_dir=args.tmp_dir,
    )


if __name__ == "__main__":
    main()
