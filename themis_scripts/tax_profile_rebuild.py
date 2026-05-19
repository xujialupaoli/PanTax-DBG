#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional


def _read_species_abundance(
    abundance_file: str | Path,
) -> Dict[str, float]:
    """
    Read a two-column TSV:
        species_id <tab> abundance

    PanTax-DBG final species_abundance.txt is expected to have NO header.
    This reader is tolerant: if a header exists, it will be skipped
    when abundance value cannot be parsed as float.
    """
    abundance_file = Path(abundance_file)
    d: Dict[str, float] = {}

    with abundance_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            sid = parts[0].strip()
            v = parts[1].strip()

            try:
                d[sid] = float(v)
            except ValueError:
                # likely header or malformed line
                continue

    return d


def rebuild_tax_profile_with_species_abundance(
    tax_profile: str | Path,
    species_abundance: str | Path,
    out_path: Optional[str | Path] = None,
    drop_root: bool = True,
    drop_strain: bool = True,
    zero_eps: float = 1e-10,
) -> str:
    """
    Rebuild a filtered ganon tax_profile.tre using final species abundance.

    Behavior required by PanTax-DBG integration:
      - default output directory = species_abundance.txt directory
      - default output filename = tax_profile.tre

    The logic follows your original standalone script.
    """
    tax_profile = Path(tax_profile)
    species_abundance = Path(species_abundance)

    if out_path is None:
        out_path = species_abundance.parent / "tax_profile.tre"
    else:
        out_path = Path(out_path)

    if not tax_profile.exists():
        raise FileNotFoundError(f"tax_profile.tre not found: {tax_profile}")
    if not species_abundance.exists():
        raise FileNotFoundError(f"species_abundance.txt not found: {species_abundance}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    sp_abund = _read_species_abundance(species_abundance)

    lines = tax_profile.read_text(encoding="utf-8").splitlines()

    # 1st pass: collect species lines
    species_lines = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        tax_level = parts[0].strip().lower()
        if tax_level != "species":
            continue

        species_id = parts[1].strip()
        id_list = parts[2].strip().split("|") if parts[2].strip() else []
        species_lines.append((species_id, set(id_list)))

    # 2nd pass: compute abundance and write
    with out_path.open("w", encoding="utf-8") as f_out:
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            tax_level_raw = parts[0].strip()
            tax_level = tax_level_raw.lower()

            if drop_root and tax_level == "root":
                continue
            if drop_strain and tax_level == "strain":
                continue

            current_id = parts[1].strip()

            abundance = 0.0
            if tax_level == "species":
                abundance = sp_abund.get(current_id, 0.0)
            else:
                for sp_id, id_set in species_lines:
                    if current_id in id_set:
                        abundance += sp_abund.get(sp_id, 0.0)

            if abs(abundance) < zero_eps:
                continue

            # keep first 4 columns
            if len(parts) >= 4:
                first_four = parts[:4]
            else:
                first_four = parts + [""] * (4 - len(parts))

            f_out.write("\t".join(first_four) + "\t" + str(abundance) + "\n")

    return str(out_path)


def add_header_to_species_abundance(
    species_abundance: str | Path,
    header1: str = "Species_TaxID",
    header2: str = "Relative_Abundance",
) -> str:
    """
    Add an English header to a two-column species_abundance file IN PLACE.

    If the first line already looks like a header, do nothing.
    """
    species_abundance = Path(species_abundance)
    if not species_abundance.exists():
        raise FileNotFoundError(f"species_abundance.txt not found: {species_abundance}")

    lines = species_abundance.read_text(encoding="utf-8").splitlines()
    if not lines:
        species_abundance.write_text(f"{header1}\t{header2}\n", encoding="utf-8")
        return str(species_abundance)

    first = lines[0].split("\t")
    has_header = False
    if len(first) >= 2:
        try:
            float(first[1])
        except ValueError:
            has_header = True

    if has_header:
        return str(species_abundance)

    new_lines = [f"{header1}\t{header2}"] + lines
    species_abundance.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return str(species_abundance)
