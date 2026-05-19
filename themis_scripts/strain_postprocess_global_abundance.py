#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PanTax-DBG strain postprocess:
- Filter ggcat strain-group table by assigned_reads cutoff (threshold/div)
- Keep only candidate species from mix_abundance_prediction.tsv
- Renormalize within-species rel_abundance_within_species
- Compute global strain abundance = species_abundance * rel_renorm
- Write:
  1) <out_prefix>.with_global.tsv  (original columns + 2 new cols)
  2) <out_prefix>.abundance.tsv    (strain_taxid \\t abundance)  [optionally with strain_name]
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import gzip
from collections import defaultdict
from typing import Dict, List, Tuple, Optional


def smart_open(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def is_number(x: str) -> bool:
    try:
        float(x)
        return True
    except Exception:
        return False


def read_mix_species_abundance(mix_path: str) -> Dict[str, float]:
    """
    mix_abundance_prediction.tsv can be:
      (a) with header: speciesID\\tabundance
      (b) no header:  species\\tabundance
    Return dict: species -> abundance
    """
    sp2ab: Dict[str, float] = {}
    with smart_open(mix_path, "rt") as f:
        first = f.readline()
        if not first:
            return sp2ab
        first = first.rstrip("\n")
        p = first.split("\t")

        # header detect: 2nd col not numeric
        if len(p) >= 2 and (not is_number(p[1])):
            pass
        else:
            if len(p) >= 2 and is_number(p[1]):
                sp2ab[p[0].strip()] = float(p[1])

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 2:
                continue
            sp = p[0].strip()
            if not sp:
                continue
            try:
                ab = float(p[1])
            except Exception:
                continue
            sp2ab[sp] = ab
    return sp2ab


def run(
    *,
    mix: str,
    input: str,
    threshold: float,
    out_prefix: str,
    div: float = 5.0,
    sum0: str = "uniform",
    precision: int = 10,
    emit_name: bool = False,
    delim: str = "\t",
) -> Tuple[str, str]:
    """
    Module entry for themis/profile.py

    Writes:
      <out_prefix>.with_global.tsv
      <out_prefix>.abundance.tsv

    Returns:
      (out_with_global, out_abund)
    """
    sp2ab = read_mix_species_abundance(mix)
    if not sp2ab:
        raise SystemExit("[Error] mix file produced no species abundances.")

    cutoff = float(threshold) / float(div)

    with smart_open(input, "rt") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise SystemExit("[Error] input file has empty header.")

        need = {"species", "assigned_reads", "rel_abundance_within_species", "strains", "strain_taxid", "strain_name"}
        missing = need.difference(set(reader.fieldnames))
        if missing:
            raise SystemExit(f"[Error] input missing columns: {sorted(missing)}")

        kept: List[dict] = []
        by_sp: Dict[str, List[dict]] = defaultdict(list)

        for row in reader:
            sp = (row.get("species") or "").strip()
            if not sp or sp not in sp2ab:
                continue

            # assigned_reads > cutoff
            try:
                ar = float(row["assigned_reads"])
            except Exception:
                continue
            if not (ar > cutoff):
                continue

            # rel_abundance_within_species > 0
            try:
                rel = float(row["rel_abundance_within_species"])
            except Exception:
                continue
            if not (rel > 0):
                continue

            # single strain only (no comma)
            strains = (row.get("strains") or "").strip()
            if (not strains) or ("," in strains):
                continue

            kept.append(row)
            by_sp[sp].append(row)

    if not kept:
        raise SystemExit("[Error] No rows left after filtering. Try lowering cutoff or check inputs.")

    # per-species sum(rel)
    sp_sum: Dict[str, float] = {}
    for sp, lst in by_sp.items():
        s = 0.0
        for r in lst:
            try:
                s += float(r["rel_abundance_within_species"])
            except Exception:
                pass
        sp_sum[sp] = s

    out_with_global = f"{out_prefix}.with_global.tsv"
    out_abund = f"{out_prefix}.abundance.tsv"

    # keep stable columns
    fieldnames = list(kept[0].keys())
    if "rel_within_species_renorm" not in fieldnames:
        fieldnames.append("rel_within_species_renorm")
    if "global_abundance" not in fieldnames:
        fieldnames.append("global_abundance")

    with open(out_with_global, "wt", encoding="utf-8", newline="") as fo, \
         open(out_abund, "wt", encoding="utf-8", newline="") as fa:

        w = csv.DictWriter(fo, fieldnames=fieldnames, delimiter=delim)
        w.writeheader()

        if emit_name:
            fa.write("strain_taxid\tstrain_name\tabundance\n")
        else:
            fa.write("strain_taxid\tabundance\n")

        # species sorted; within species sort by global abundance desc
        for sp in sorted(by_sp.keys()):
            lst = by_sp[sp]
            denom = sp_sum.get(sp, 0.0)
            n = len(lst) if lst else 1

            tmp: List[Tuple[float, dict]] = []
            for r in lst:
                rel = float(r["rel_abundance_within_species"])
                if denom > 0:
                    rel_new = rel / denom
                else:
                    if sum0 == "uniform":
                        rel_new = 1.0 / n
                    elif sum0 == "zero":
                        rel_new = 0.0
                    else:
                        rel_new = rel

                global_ab = sp2ab[sp] * rel_new
                r2 = dict(r)
                r2["rel_within_species_renorm"] = f"{rel_new:.{precision}f}"
                r2["global_abundance"] = f"{global_ab:.{precision}f}"
                tmp.append((global_ab, r2))

            tmp.sort(key=lambda x: x[0], reverse=True)

            for global_ab, r2 in tmp:
                w.writerow(r2)
                taxid = (r2.get("strain_taxid") or "").strip()
                name = (r2.get("strain_name") or "").strip()
                if not taxid:
                    continue
                if emit_name:
                    fa.write(f"{taxid}\t{name}\t{global_ab:.{precision}f}\n")
                else:
                    fa.write(f"{taxid}\t{global_ab:.{precision}f}\n")

    # optional log
    print(f"[Completed] wrote: {out_with_global}")
    print(f"[Completed] wrote: {out_abund}")
    print(f"[Info] cutoff = threshold/div = {threshold}/{div} = {cutoff}")

    return out_with_global, out_abund


def main():
    ap = argparse.ArgumentParser(
        description="PanTax-DBG: strain filter + global abundance postprocess.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-m", "--mix", required=True,
                    help="mix_abundance_prediction.tsv (species, abundance)")
    ap.add_argument("-i", "--input", required=True,
                    help="Input strain-group table (TOPK result, with strain_taxid & strain_name)")
    ap.add_argument("--threshold", type=float, required=True,
                    help="threshold computed in themis/profile.py")
    ap.add_argument("--div", type=float, default=5.0,
                    help="cutoff = threshold/div (default=5.0)")
    ap.add_argument("-o", "--out_prefix", required=True,
                    help="Output prefix (writes <prefix>.with_global.tsv and <prefix>.abundance.tsv)")
    ap.add_argument("--sum0", choices=["uniform", "zero", "keep"], default="uniform",
                    help="If a species sum(rel)=0: uniform|zero|keep.")
    ap.add_argument("--precision", type=int, default=10,
                    help="Decimal places for numeric outputs.")
    ap.add_argument("--emit-name", action="store_true",
                    help="If set, output abundance TSV: strain_taxid\\tstrain_name\\tabundance.")
    ap.add_argument("--delim", default="\t", help="TSV delimiter")
    args = ap.parse_args()

    run(
        mix=args.mix,
        input=args.input,
        threshold=args.threshold,
        out_prefix=args.out_prefix,
        div=args.div,
        sum0=args.sum0,
        precision=args.precision,
        emit_name=args.emit_name,
        delim=args.delim,
    )


if __name__ == "__main__":
    main()
