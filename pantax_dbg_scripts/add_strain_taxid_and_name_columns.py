#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PanTax-DBG helper:
Append strain_taxid and strain_name columns to strain_group_abundance.tsv
by mapping "new_strain_taxid" (e.g., 122_1) from ganon_species_strain_topk.taxid.tsv.

Mapping file example (TSV):
species_taxid  strain_taxid  strain_name  new_strain_taxid  genome_path
303            4111914       GCF_...       303_1             /path/...

Input abundance example (TSV):
species  group_id  assigned_reads  rel_abundance_within_species  strains
122      G1        ...            ...                           122_1
139      G3        ...            ...                           139_1,139_3

Output: append two columns:
strain_taxid   strain_name
(If multiple strains in a row, join with comma in the same order.)
"""

import argparse
import csv
import gzip
from typing import Dict, List, Tuple


def smart_open(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def split_ids(s: str) -> List[str]:
    """Split strains field by comma. If no comma, returns single element."""
    s = (s or "").strip()
    if not s:
        return []
    parts = [x.strip() for x in s.split(",")]
    return [p for p in parts if p]


def read_new2info(map_path: str, delim: str = "\t") -> Dict[str, Tuple[str, str]]:
    """
    Read mapping file and build:
      new_strain_taxid -> (strain_taxid, strain_name)
    """
    with smart_open(map_path, "rt") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise RuntimeError(f"Empty header in map file: {map_path}")

        fn = [x.strip() for x in reader.fieldnames]
        need = {"new_strain_taxid", "strain_taxid", "strain_name"}
        missing = need.difference(set(fn))
        if missing:
            raise RuntimeError(
                f"Map file must contain columns: {sorted(need)}. Missing: {sorted(missing)}. Got: {fn}"
            )

        m: Dict[str, Tuple[str, str]] = {}
        dup = 0
        for row in reader:
            new_id = (row.get("new_strain_taxid") or "").strip()
            taxid = (row.get("strain_taxid") or "").strip()
            name = (row.get("strain_name") or "").strip()
            if not new_id:
                continue
            if new_id in m:
                dup += 1
                continue
            m[new_id] = (taxid, name)

        # 不强制报错，但给使用者一个提示会更安全
        if dup:
            # 不用 print 到 stdout，以免影响上游捕获
            import sys
            print(f"[Warning] {dup} duplicated new_strain_taxid in map were ignored (kept first).", file=sys.stderr)

        return m

def run(
    *,
    input: str,
    map: str,
    out: str,
    new_col: str = "strains",
    add_taxid_col: str = "strain_taxid",
    add_name_col: str = "strain_name",
    missing: str = "NA",
):
    new2info = read_new2info(map, "\t")

    with smart_open(input, "rt") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError(f"Empty header in input file: {input}")

        if new_col not in reader.fieldnames:
            raise RuntimeError(
                f"Input file missing column '{new_col}'. Header: {reader.fieldnames}"
            )

        out_fields = list(reader.fieldnames)
        if add_taxid_col not in out_fields:
            out_fields.append(add_taxid_col)
        if add_name_col not in out_fields:
            out_fields.append(add_name_col)

        with open(out, "wt", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t")
            writer.writeheader()

            for row in reader:
                strains_field = row.get(new_col, "")
                new_ids = split_ids(strains_field)

                if not new_ids:
                    row[add_taxid_col] = missing
                    row[add_name_col] = missing
                    writer.writerow(row)
                    continue

                taxids: List[str] = []
                names: List[str] = []
                for nid in new_ids:
                    taxid, name = new2info.get(nid, (missing, missing))
                    taxids.append(taxid)
                    names.append(name)

                row[add_taxid_col] = ",".join(taxids)
                row[add_name_col] = ",".join(names)
                writer.writerow(row)

    return out

def main():
    ap = argparse.ArgumentParser(
        description="Append strain_taxid and strain_name columns to strain_group_abundance.tsv "
                    "using ganon_species_strain_topk.taxid.tsv mapping.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-i", "--input", required=True,
                    help="query_*.strain_group_abundance.tsv (or .gz)")
    ap.add_argument("-m", "--map", required=True,
                    help="ganon_species_strain_topk.taxid.tsv (or .gz)")
    ap.add_argument("-o", "--out", required=True,
                    help="Output TSV path")
    ap.add_argument("--new-col", default="strains",
                    help="Column in input containing new_strain_taxid list (comma-separated if multiple).")
    ap.add_argument("--add-taxid-col", default="strain_taxid",
                    help="New column name to append for strain taxids.")
    ap.add_argument("--add-name-col", default="strain_name",
                    help="New column name to append for strain names (e.g., GCF...).")
    ap.add_argument("--missing", default="NA",
                    help="Value to use when mapping is missing.")
    args = ap.parse_args()

    run(
        input=args.input,
        map=args.map,
        out=args.out,
        new_col=args.new_col,
        add_taxid_col=args.add_taxid_col,
        add_name_col=args.add_name_col,
        missing=args.missing,
    )



if __name__ == "__main__":
    main()
