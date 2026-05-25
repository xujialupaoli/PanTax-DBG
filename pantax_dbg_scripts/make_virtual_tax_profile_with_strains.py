#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict


def load_tax(tax_path):
    tax = {}
    with open(tax_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            tid, parent, rank, name = parts[0], parts[1], parts[2], parts[3]
            tax[tid] = (parent, rank, name)
    return tax


def load_topk(tsv_path, topk=10):
    d = defaultdict(list)
    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            sp, st, u, m = parts[:4]
            try:
                d[sp].append((st, int(u), int(m)))
            except ValueError:
                continue
    for sp in list(d.keys()):
        d[sp] = d[sp][:topk]
    return d


def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def allocate_reads(total_reads, weights):
    K = len(weights)
    if K == 0:
        return []
    if total_reads <= 0:
        return [0] * K

    s = sum(weights)
    if s <= 0:
        base = total_reads // K
        rem = total_reads - base * K
        out = [base] * K
        for i in range(rem):
            out[i] += 1
        return out

    raw = [total_reads * (w / s) for w in weights]
    floor = [int(x) for x in raw]
    rem = total_reads - sum(floor)
    frac = sorted([(raw[i] - floor[i], i) for i in range(K)], reverse=True)
    out = floor[:]
    for j in range(rem):
        out[frac[j][1]] += 1
    return out


def run(species_tre, species_topk_tsv, tax, out="tax_profile.tre", topk=10, weight="sum"):
    tax_map = load_tax(tax)
    top = load_topk(species_topk_tsv, topk=topk)

    with open(species_tre, "r", encoding="utf-8", errors="replace") as f:
        lines = [line.rstrip("\n") for line in f]

    root_total_reads = None
    for ln in lines:
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < 9:
            continue
        if parts[0] == "root" and parts[1] == "1":
            root_total_reads = safe_int(parts[7], None)
            break
    if root_total_reads is None or root_total_reads <= 0:
        root_total_reads = 0

    out_lines = []
    inserted = 0

    for ln in lines:
        out_lines.append(ln)
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < 9:
            continue

        rk, tid = parts[0], parts[1]
        lineage = parts[2]
        reads = safe_int(parts[7], 0)

        if rk != "species" or tid not in top:
            continue

        strains = top[tid]
        weights = []
        for st, u, m in strains:
            if weight == "unique":
                w = u
            elif weight == "multi":
                w = m
            elif weight == "unique_plus_0.1multi":
                w = u + 0.1 * m
            else:
                w = u + m
            weights.append(float(w))

        alloc = allocate_reads(reads, weights)

        for (st, u, m), st_reads in zip(strains, alloc):
            st_name = tax_map.get(st, ("", "", st))[2] if st in tax_map else st
            st_lineage = f"{lineage}|{st}"
            pct = (st_reads / root_total_reads) * 100.0 if root_total_reads > 0 else 0.0
            out_lines.append("\t".join([
                "strain", st, st_lineage, st_name, "0", "0", "0", str(st_reads), f"{pct:.5f}",
            ]))
            inserted += 1

    with open(out, "w", encoding="utf-8") as f:
        for ln in out_lines:
            f.write(ln + "\n")

    print(f"[PanTax-DBG] wrote virtual strain-augmented TRE: {out} (inserted_strains={inserted})")
    return str(out)


def main():
    ap = argparse.ArgumentParser(description="Insert top strains into species_profile.tre.")
    ap.add_argument("species_tre")
    ap.add_argument("species_top10_tsv")
    ap.add_argument("tax")
    ap.add_argument("-o", "--out", default="tax_profile.tre")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--weight", choices=["unique", "multi", "sum", "unique_plus_0.1multi"], default="sum")
    args = ap.parse_args()
    run(args.species_tre, args.species_top10_tsv, args.tax, args.out, args.topk, args.weight)


if __name__ == "__main__":
    main()
