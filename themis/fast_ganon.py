#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import time
from typing import Dict

from .utils import ensure_dir, run_cmd
from .paths import get_internal_binary, get_dbg_ganon

from themis_scripts import make_species_rep_from_all
from themis_scripts import species_reassign_em_strict
from themis_scripts import make_virtual_tax_profile_with_strains
from themis_scripts import ganon_species_process
from themis_scripts import parse_strain


GANON_STRAIN_TOPK = 10
GANON_MULTI_PICK = 3
EM_THRESHOLD = 1e-7
EM_MAX_ITER = 0
VIRTUAL_STRAIN_TOPK = 10
VIRTUAL_STRAIN_WEIGHT = "sum"
MAX_RETRIES = 3


def _required_outputs_exist(query_dir: Path) -> bool:
    required = [
        query_dir / "species_abundance.txt",
        query_dir / "strain_abundance.txt",
        query_dir / "predict_spy.ID.abundance",
        query_dir / "tax_profile.tre",
    ]
    return all(p.exists() and p.stat().st_size > 0 for p in required)


def _write_predict_spy(species_abundance: Path, out_abundance: Path, out_ids: Path | None = None) -> None:
    out_abundance.parent.mkdir(parents=True, exist_ok=True)
    ids = []
    with species_abundance.open("r", encoding="utf-8", errors="replace") as fin, \
         out_abundance.open("w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                float(parts[1])
            except ValueError:
                continue
            fout.write(f"{parts[0]}\t{parts[1]}\n")
            ids.append(parts[0])

    if out_ids is not None:
        with out_ids.open("w", encoding="utf-8") as f:
            for sid in ids:
                f.write(sid + "\n")


def run_paired_prehit(
    *,
    db_prefix: str,
    read1: str,
    read2: str,
    query_dir: str | Path,
    threads: int = 8,
    report_type: str = "abundance",
    keep_species_all: bool = False,
    force: bool = False,
) -> Dict[str, str]:
    query_dir = Path(query_dir)
    ensure_dir(query_dir)

    tax_path = Path(f"{db_prefix}.tax")
    ibf_path = Path(f"{db_prefix}.hibf")
    if not tax_path.exists():
        raise FileNotFoundError(f"[PanTax-DBG][error] DB tax file not found: {tax_path}")
    if not ibf_path.exists():
        raise FileNotFoundError(f"[PanTax-DBG][error] DB HIBF file not found: {ibf_path}")

    species_abundance = query_dir / "species_abundance.txt"
    strain_abundance = query_dir / "strain_abundance.txt"
    predict_spy = query_dir / "predict_spy.ID.abundance"
    predict_ids = query_dir / "predict_spy_number.ID"
    tax_profile = query_dir / "tax_profile.tre"

    if (not force) and _required_outputs_exist(query_dir):
        print("[PanTax-DBG] Found existing fast DBG-ganon outputs; skipping pre-screening.")
        return {
            "tax_profile": str(tax_profile),
            "species_abundance": str(species_abundance),
            "strain_abundance": str(strain_abundance),
            "predict_spy": str(predict_spy),
        }

    ganon_classify = get_internal_binary("ganon-classify")
    ganon_frontend = get_dbg_ganon()

    results_prefix = query_dir / "results"
    species_all = query_dir / "results.species.all"
    species_topk_tsv = query_dir / f"results.species_top{GANON_STRAIN_TOPK}_strains.tsv"
    expected_topk = query_dir / "results.species_top10_strains.tsv"

    cmd = [
        ganon_classify,
        "--paired-reads", f"{read1},{read2}",
        "--ibf", f"{db_prefix}.hibf",
        "--tax", f"{db_prefix}.tax",
        "--output-prefix", str(results_prefix),
        "--threads", str(threads),
        "--hibf",
        "--output-species-all",
        "--output-species-topk",
        "--species-topk", str(GANON_STRAIN_TOPK),
        "--multi-pick", str(GANON_MULTI_PICK),
    ]

    last_error = None
    for attempt in range(1, MAX_RETRIES + 1):
        if attempt > 1:
            print(f"[PanTax-DBG][warning] retrying DBG-ganon classify: attempt {attempt}/{MAX_RETRIES}")
            for p in query_dir.glob("results*"):
                if p.is_file():
                    p.unlink(missing_ok=True)
            time.sleep(30)

        try:
            run_cmd(cmd, "[PanTax-DBG][DBG-ganon classify]", echo=True, silence=False)
        except Exception as e:
            last_error = e

        if species_all.exists() and species_all.stat().st_size > 0:
            break
    else:
        raise RuntimeError(f"[PanTax-DBG][error] DBG-ganon classify failed after {MAX_RETRIES} attempts: {last_error}")

    if expected_topk.exists():
        species_topk_tsv = expected_topk
    if not species_topk_tsv.exists():
        candidates = sorted(query_dir.glob("results.species_top*_strains.tsv"))
        if candidates:
            species_topk_tsv = candidates[0]
        else:
            raise FileNotFoundError("[PanTax-DBG][error] missing results.species_top*_strains.tsv from DBG-ganon classify")

    make_species_rep_from_all.run(str(species_all), str(tax_path), str(query_dir / "results.species.rep"))

    em_prefix = query_dir / "results.species.EM"
    em_rep = species_reassign_em_strict.run(str(species_all), str(em_prefix), EM_THRESHOLD, EM_MAX_ITER, str(tax_path))

    species_profile_prefix = query_dir / "species_profile"
    report_cmd = [
        ganon_frontend,
        "report",
        "-i", str(em_rep),
        "--db-prefix", str(db_prefix),
        "--output-prefix", str(species_profile_prefix),
        "--report-type", report_type,
        "-r", "all",
    ]
    run_cmd(report_cmd, "[PanTax-DBG][DBG-ganon report]", echo=True, silence=False)

    species_profile_tre = query_dir / "species_profile.tre"
    if not species_profile_tre.exists():
        raise FileNotFoundError(f"[PanTax-DBG][error] ganon report did not create {species_profile_tre}")

    make_virtual_tax_profile_with_strains.run(
        str(species_profile_tre),
        str(species_topk_tsv),
        str(tax_path),
        str(tax_profile),
        VIRTUAL_STRAIN_TOPK,
        VIRTUAL_STRAIN_WEIGHT,
    )

    ganon_species_process.run(str(tax_profile), str(species_abundance))
    parse_strain.run(str(tax_profile), str(strain_abundance))
    _write_predict_spy(species_abundance, predict_spy, predict_ids)

    if not keep_species_all:
        species_all.unlink(missing_ok=True)

    print("[PanTax-DBG] fast DBG-ganon pre-screening finished.")
    return {
        "tax_profile": str(tax_profile),
        "species_abundance": str(species_abundance),
        "strain_abundance": str(strain_abundance),
        "predict_spy": str(predict_spy),
    }
