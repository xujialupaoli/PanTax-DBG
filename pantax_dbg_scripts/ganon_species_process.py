#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse


def run(report_file, output="species_abundance.txt"):
    """
    Extract species-level relative abundance from a ganon .tre file.

    Expected TRE columns:
      rank, taxid, lineage, name, matches, unique, lca, reads, pct

    Output:
      speciesID <tab> abundance
    """
    report_file = Path(report_file)
    output = Path(output)

    tax_profile_dict = {}
    with report_file.open("r", encoding="utf-8", errors="replace") as f_in:
        for line in f_in:
            if not line.strip():
                continue
            if not line.startswith("species\t"):
                continue
            tokens = line.rstrip("\n").split("\t")
            if len(tokens) < 9:
                continue
            species_taxid = tokens[1].strip()
            if not species_taxid:
                continue
            try:
                abundance = float(tokens[8]) / 100.0
            except ValueError:
                continue
            tax_profile_dict[species_taxid] = abundance

    sorted_items = sorted(tax_profile_dict.items(), key=lambda item: item[1], reverse=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as f_out:
        f_out.write("speciesID\tabundance\n")
        for k, v in sorted_items:
            f_out.write(f"{k}\t{v}\n")

    return str(output)


def main():
    ap = argparse.ArgumentParser(description="Extract species abundance from a ganon TRE file.")
    ap.add_argument("report_file")
    ap.add_argument("-o", "--output", default="species_abundance.txt")
    args = ap.parse_args()
    run(args.report_file, args.output)


if __name__ == "__main__":
    main()
