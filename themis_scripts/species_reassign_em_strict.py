#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict
from pathlib import Path
import sys


def load_tax_name_rank(db_tax):
    name = {}
    rank = {}
    if db_tax is None:
        return name, rank
    with open(db_tax, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                node, parent, r, n = fields[0], fields[1], fields[2], fields[3]
                name[node] = n
                rank[node] = r
    return name, rank


def get_top_match(matches, prob):
    t0, k0 = matches[0]
    best_t, best_k = t0, k0
    best_p = prob.get(t0, 0.0)
    for t, k in matches[1:]:
        p = prob.get(t, 0.0)
        if p > best_p:
            best_p = p
            best_t = t
            best_k = k
    return best_t, best_k


def run(all_path, out_prefix, threshold=1e-7, max_iter=0, tax_path=None):
    all_path = Path(all_path)
    out_prefix = str(out_prefix)

    print("[PanTax-DBG] Reassigning reads with species-level hard-EM", file=sys.stderr)

    targets = defaultdict(lambda: len(targets))
    read_matches = {}
    init_weight = defaultdict(int)

    with all_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rid, sp, kc = parts[0], parts[1], parts[2]
            tid = targets[sp]
            read_matches.setdefault(rid, []).append((tid, int(kc)))

    targets_rev = {v: k for k, v in targets.items()}
    total_reads = len(read_matches)
    if total_reads == 0:
        raise SystemExit("[PanTax-DBG][error] no reads in species.all")

    total_unique = 0
    for matches in read_matches.values():
        if len(matches) == 1:
            total_unique += 1
            init_weight[matches[0][0]] += 1

    if total_unique == 0:
        prob = {tid: 1.0 / len(targets) for tid in targets_rev.keys()}
    else:
        prob = {tid: init_weight.get(tid, 0) / total_unique for tid in targets_rev.keys()}

    reassigned = defaultdict(int)
    em_it = 0
    while True:
        reassigned = defaultdict(int)
        for tid, c in init_weight.items():
            reassigned[tid] += c

        for matches in read_matches.values():
            if len(matches) > 1:
                tid, _ = get_top_match(matches, prob)
                reassigned[tid] += 1

        diff = 0.0
        for tid in targets_rev.keys():
            newp = reassigned.get(tid, 0) / total_reads
            diff += abs(prob.get(tid, 0.0) - newp)
            prob[tid] = newp

        em_it += 1
        print(f"[PanTax-DBG] hard-EM iteration {em_it}: diff={diff:.6g}", file=sys.stderr)

        if diff <= float(threshold):
            break
        if int(max_iter) > 0 and em_it >= int(max_iter):
            break

    direct = defaultdict(int)
    for matches in read_matches.values():
        for tid, _ in matches:
            direct[tid] += 1

    tax_name, tax_rank = load_tax_name_rank(tax_path)
    out_rep = out_prefix + ".rep"
    with open(out_rep, "w", encoding="utf-8") as out:
        written = 0
        for tid, sp in targets_rev.items():
            assigned = reassigned.get(tid, 0)
            if assigned <= 0:
                continue
            r = tax_rank.get(sp, "species")
            n = tax_name.get(sp, sp)
            out.write(f"H1\t{sp}\t{direct.get(tid,0)}\t{assigned}\t0\t{r}\t{n}\n")
            written += 1
        out.write(f"#total_classified\t{total_reads}\n")
        out.write("#total_unclassified\t0\n")

    print(f"[PanTax-DBG] wrote {out_rep} (targets={written})", file=sys.stderr)
    return out_rep


def main():
    ap = argparse.ArgumentParser(description="Species-level hard-EM reassignment for PanTax-DBG.")
    ap.add_argument("species_all")
    ap.add_argument("out_prefix")
    ap.add_argument("threshold", type=float)
    ap.add_argument("max_iter", type=int, nargs="?", default=0)
    ap.add_argument("tax", nargs="?", default=None)
    args = ap.parse_args()
    run(args.species_all, args.out_prefix, args.threshold, args.max_iter, args.tax)


if __name__ == "__main__":
    main()
