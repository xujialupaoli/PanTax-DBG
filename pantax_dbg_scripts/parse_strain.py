#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import sys


def parse_strain_from_tre(tre_file: str, out_path: str):
    items = []

    tre_file = Path(tre_file)
    out_path = Path(out_path)

    with tre_file.open("r", encoding="utf-8", errors="ignore") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or not line.startswith("strain\t"):
                continue

            toks = line.split("\t")
            if len(toks) < 9:
                print(f"[PanTax-DBG][warning] line {line_num}: fewer than 9 columns; skipped", file=sys.stderr)
                continue

            sid = toks[1].strip()
            if not sid and len(toks) >= 3:
                sid = toks[2].split("|")[-1].strip()
            if not sid and len(toks) >= 4:
                sid = toks[3].strip()
            if not sid:
                print(f"[PanTax-DBG][warning] line {line_num}: cannot extract strain id; skipped", file=sys.stderr)
                continue

            try:
                abund = float(toks[8]) / 100.0
            except ValueError:
                print(f"[PanTax-DBG][warning] line {line_num}: invalid abundance {toks[8]!r}; skipped", file=sys.stderr)
                continue

            items.append((sid, abund))

    items.sort(key=lambda x: x[1], reverse=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write("strain_taxid\tabundance\n")
        for sid, abund in items:
            out.write(f"{sid}\t{abund:.10g}\n")

    print(f"[PanTax-DBG] Parsed {len(items)} strain rows -> {out_path}")
    return str(out_path)


def run(input, output):
    return parse_strain_from_tre(input, output)


def main():
    parser = argparse.ArgumentParser(description="Extract strain abundance from a TRE file.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    run(args.input, args.output)


if __name__ == "__main__":
    main()
