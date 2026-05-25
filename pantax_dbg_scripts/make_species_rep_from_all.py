#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict


def load_tax_name_rank(tax_path):
    rank = {}
    name = {}
    with open(tax_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                node = fields[0]
                rank[node] = fields[2]
                name[node] = fields[3]
    return rank, name


def run(species_all, tax, out_rep):
    rank_map, name_map = load_tax_name_rank(tax)

    matches = defaultdict(int)
    unique_reads = defaultdict(int)

    cur_read = None
    cur_species = []

    def flush():
        nonlocal cur_read, cur_species
        if cur_read is None or not cur_species:
            return
        for sp in cur_species:
            matches[sp] += 1
        if len(cur_species) == 1:
            unique_reads[cur_species[0]] += 1

    with open(species_all, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            rid, sp, kc = line.rstrip("\n").split("\t")[:3]

            if cur_read is None:
                cur_read = rid

            if rid != cur_read:
                flush()
                cur_read = rid
                cur_species = []

            cur_species.append(sp)

        flush()

    with open(out_rep, "w", encoding="utf-8") as out:
        for sp in matches.keys():
            r = rank_map.get(sp, "no rank")
            n = name_map.get(sp, sp)
            out.write(f"H1\t{sp}\t{matches[sp]}\t{unique_reads.get(sp,0)}\t0\t{r}\t{n}\n")

    print(f"[PanTax-DBG] wrote {out_rep} (targets={len(matches)})")
    return out_rep


def main():
    ap = argparse.ArgumentParser(description="Create species-level ganon rep from results.species.all.")
    ap.add_argument("species_all")
    ap.add_argument("tax")
    ap.add_argument("out_rep")
    args = ap.parse_args()
    run(args.species_all, args.tax, args.out_rep)


if __name__ == "__main__":
    main()
