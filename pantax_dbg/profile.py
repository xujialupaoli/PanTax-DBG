#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path
import math
import shutil
import argparse
import importlib
import sys

from .utils import ensure_dir
from . import ganon_wrapper, ggcat_wrapper, fast_ganon


def _ts(mod: str):
    return importlib.import_module(f"pantax_dbg_scripts.{mod}")


def _cleanup_intermediate_dirs(*paths):
    """Remove internal work directories after a successful paired-end run."""
    for path in paths:
        p = Path(path)
        if p.exists() and p.is_dir():
            shutil.rmtree(p, ignore_errors=True)

# ==============================================================================
# Helpers
# ==============================================================================

def _filter_and_renormalize_abundance(input_path, output_path, min_val):
    """Filter a two-column abundance table and renormalize retained values."""
    items = []
    total_abund = 0.0
    
    if not Path(input_path).exists():
        return

    with open(input_path, "r", encoding="utf-8") as fin:
        lines = [l.strip() for l in fin if l.strip()]
        
    if not lines:
        with open(output_path, "w", encoding="utf-8") as fout:
            fout.write("speciesID\tabundance\n")
        return

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        
        sid = parts[0]
        try:
            val = float(parts[1])
            if val > min_val:
                items.append((sid, val))
                total_abund += val
        except ValueError:
            continue

    if total_abund > 0:
        norm_factor = 1.0 / total_abund
    else:
        norm_factor = 0.0

    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("speciesID\tabundance\n")
        items.sort(key=lambda x: x[1], reverse=True)
        for sid, val in items:
            new_val = val * norm_factor
            if new_val > 0:
                fout.write(f"{sid}\t{new_val}\n")
    
    print(f"[PanTax-DBG] Filtered species (>{min_val}) and renormalized: {output_path}")


def _count_non_empty_lines_skip_header(path):
    """Count rows whose second column contains a numeric value."""
    count = 0
    if not Path(path).exists():
        return 0
        
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line: 
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    float(parts[1])
                    count += 1
                except ValueError:
                    continue
    return count


# ==============================================================================
# Main Logic
# ==============================================================================

def run(
    reads,
    single,
    db_prefix,
    out_prefix,
    report_type="abundance",
    file_info=None,
    threads=8,
    k=31,
    strain_min_reads=5.0,
    strain_topk=5,
    strain_div=5.0,
    r1_large_n=1000,
    r1_topk_large=3,
    r1_min_abundance=0.0,
    species_min_abundance=0.0,
    ggcat_temp_dir=None,
    ggcat_memory=None,
    ggcat_prefer_memory=False,
    ggcat_disk_optimization_level=None,
    ggcat_intermediate_compression_level=None,
):
    # Public PanTax-DBG release: paired-end short-read workflow only.
    # Internal single-read code paths are intentionally not exposed in this release.
    if len(reads) != 2:
        raise SystemExit(
            "[PanTax-DBG][error] profile expects paired-end short-read input. "
            "Please provide exactly two read files with -r R1 -r R2."
        )

    
    out_prefix = Path(out_prefix).absolute()
    ensure_dir(out_prefix)

    if single:
        return run_single(
            reads=reads,
            db_prefix=db_prefix,
            out_prefix=out_prefix,
            report_type=report_type,
            ref_info=file_info,
            threads=threads,
            k=k,
            species_min_abundance=species_min_abundance,
            ggcat_temp_dir=ggcat_temp_dir,
            ggcat_memory=ggcat_memory,
            ggcat_prefer_memory=ggcat_prefer_memory,
            ggcat_disk_optimization_level=ggcat_disk_optimization_level,
            ggcat_intermediate_compression_level=ggcat_intermediate_compression_level,
        )
    else:
        return run_paired(
            reads=reads,
            db_prefix=db_prefix,
            out_prefix=out_prefix,
            report_type=report_type,
            ref_info=file_info,
            threads=threads,
            k=k,
            strain_min_reads=strain_min_reads,
            strain_topk=strain_topk,
            strain_div=strain_div,
            r1_large_n=r1_large_n,
            r1_topk_large=r1_topk_large,
            r1_min_abundance=r1_min_abundance,
            species_min_abundance=species_min_abundance,
            ggcat_temp_dir=ggcat_temp_dir,
            ggcat_memory=ggcat_memory,
            ggcat_prefer_memory=ggcat_prefer_memory,
            ggcat_disk_optimization_level=ggcat_disk_optimization_level,
            ggcat_intermediate_compression_level=ggcat_intermediate_compression_level,
        )

# =========================
# Paired Mode
# =========================

def run_paired(reads, db_prefix, out_prefix,
               report_type, ref_info, threads, k,
               strain_min_reads=5.0,
               strain_topk=5,
               strain_div=5.0,
               r1_large_n=1000,
               r1_topk_large=3,
               r1_min_abundance=0.0,
               species_min_abundance=0.0,
               ggcat_temp_dir=None,
               ggcat_memory=None,
               ggcat_prefer_memory=False,
               ggcat_disk_optimization_level=None,
               ggcat_intermediate_compression_level=None,
):
    if len(reads) != 2:
        raise SystemExit("[PanTax-DBG] Paired-end mode requires two -r reads (R1 and R2).")
    if not ref_info:
        raise SystemExit("[PanTax-DBG] --ref-info is required.")

    r1, r2 = map(str, reads)

    work_ganon = out_prefix / "query_r1"
    work_db = out_prefix / "database_filter_ccDBG"
    work_q = out_prefix / "res_pantax_dbg"

    ensure_dir(work_ganon)
    ensure_dir(work_db)
    ensure_dir(work_q)

    # -------------------------------------------------------------
    # 1) Fast DBG-ganon pre-screening
    # -------------------------------------------------------------
    ganon_tre = work_ganon / "tax_profile.tre"
    ganon_species = work_ganon / "species_abundance.txt"
    ganon_strain = work_ganon / "strain_abundance.txt"
    ganon_predict_spy = work_ganon / "predict_spy.ID.abundance"

    all_ganon_files_exist = (
        ganon_tre.exists() and
        ganon_species.exists() and
        ganon_strain.exists() and
        ganon_predict_spy.exists()
    )

    if all_ganon_files_exist:
        print("[PanTax-DBG] Found existing fast DBG-ganon outputs. Skipping pre-screening.")
    else:
        print("[PanTax-DBG] Running built-in fast DBG-ganon pre-screening...")
        fast_ganon.run_paired_prehit(
            db_prefix=db_prefix,
            read1=r1,
            read2=r2,
            query_dir=work_ganon,
            threads=threads,
            report_type=report_type,
            keep_species_all=False,
            force=False,
        )

    # Filter first-stage candidates before ccDBG refinement.
    target_species_file = ganon_species

    if r1_min_abundance >= 0:
        filtered_species_file = work_ganon / f"species_abundance_r1_min{r1_min_abundance}.txt"
        
        print(f"[PanTax-DBG] Filtering r1 species > {r1_min_abundance} for ccDBG...")
        _filter_and_renormalize_abundance(str(ganon_species), str(filtered_species_file), r1_min_abundance)
        
        target_species_file = filtered_species_file

    # Select the candidate-genome panel for graph construction.
    line_count = _count_non_empty_lines_skip_header(str(target_species_file))
    print(f"[PanTax-DBG] Effective species count for DBG: {line_count}")

    tmp_id_table = work_q / "tmp_id_table.tsv"

    if line_count > r1_large_n:
        top_tsv = work_db / "ganon_species_strain_top.tsv"
        _ts("make_ganon_pred_symlinks_topk_singleton_filter").run(
            ref_info=str(ref_info),
            ganon_species=str(target_species_file),
            ganon_strain=str(ganon_strain),
            top_k=r1_topk_large,
            singleton_min_abund=7e-7,
            out_tsv=str(top_tsv),
        )
        shutil.copyfile(top_tsv, tmp_id_table)
    elif line_count >= 1:
        sel_tsv = work_db / "ganon_species_strain_selected.tsv"
        _ts("make_ganon_pred_symlinks").run(
            ref_info=str(ref_info),
            ganon_species=str(target_species_file),
            ganon_strain=str(ganon_strain),
            out_tsv=str(sel_tsv),
        )
        shutil.copyfile(sel_tsv, tmp_id_table)
    else:
        print("[PanTax-DBG] Warning: No species remained after filtering. Writing an empty final species profile.")
        final_abundance = out_prefix / "species_abundance.txt"
        with open(final_abundance, "w") as f:
            f.write("speciesID\tabundance\n")
        _cleanup_intermediate_dirs(work_ganon, work_db, work_q, out_prefix / ".temp_files")
        return str(final_abundance)

    # Build the mapping from graph colors to candidate genomes.
    color_map = work_db / "color_mapping.in"

    if color_map.exists():
        print(f"[PanTax-DBG] Found existing color_mapping: {color_map}. Skipping.")
    else:
        _make_color_mapping_from_mapping_tsv(str(tmp_id_table), str(color_map))

    # Build a sample-specific ccDBG for the retained genomes.
    ggcat_db = work_db / "ggcatDB.fasta.lz4"

    if ggcat_db.exists():
        print(f"[PanTax-DBG] Found existing ggcatDB: {ggcat_db}. Skipping build.")
    else:
        ggcat_wrapper.run_build(
            k=k,
            threads=threads,
            color_mapping=str(color_map),
            output=str(ggcat_db),
            temp_dir=ggcat_temp_dir,
            memory=ggcat_memory,
            prefer_memory=ggcat_prefer_memory,
            disk_optimization_level=ggcat_disk_optimization_level,
            intermediate_compression_level=ggcat_intermediate_compression_level,
        )

    # Interleave and preprocess paired reads for graph queries.
    ggcat_reads = work_q / "reads_for_ggcat.fastq"

    if ggcat_reads.exists():
        print(f"[PanTax-DBG] Found existing reads for ggcat: {ggcat_reads}. Skipping fastp.")
    else:
        _prepare_ggcat_reads_paired(r1, r2, str(ggcat_reads), threads=threads)

    # Query reads against the sample-specific ccDBG.
    ggcat_prefix = work_q / "query_ggcatDB"

    check_ggcat_outputs = [
        Path(f"{ggcat_prefix}.species_counts.tsv"),
        Path(f"{ggcat_prefix}.strain_group_abundance.tsv")
    ]
    if all(p.exists() for p in check_ggcat_outputs):
        print(f"[PanTax-DBG] Found existing GGCAT query results. Skipping query.")
    else:
        ggcat_wrapper.run_query(
            db=str(ggcat_db),
            reads=[str(ggcat_reads)],
            k=k,
            threads=threads,
            out_prefix=str(ggcat_prefix),
            single=False,
        )

    # 8) Mix predictions
    temp_mix_abundance = out_prefix / "species_abundance.raw_mix.txt"
    
    _run_threshold_and_mix_for_paired(
        ggcat_prefix=str(ggcat_prefix),
        tmp_id_table=str(tmp_id_table),
        ganon_predict_spy=str(ganon_predict_spy),
        tre_file=str(ganon_tre),
        reads_path=str(ggcat_reads),
        output=str(temp_mix_abundance),
        strain_min_reads=strain_min_reads,
        strain_topk=strain_topk,
        strain_div=strain_div,
    )

    # =========================================================
    # [FILTER 2] Final Result Filtering & Renormalization
    # =========================================================
    
    final_abundance = out_prefix / "species_abundance.txt"

    _filter_and_renormalize_abundance(
        input_path=str(temp_mix_abundance), 
        output_path=str(final_abundance), 
        min_val=species_min_abundance
    )
    
    # 9) Rebuild Tax Profile
    _ts("tax_profile_rebuild").rebuild_tax_profile_with_species_abundance(
        tax_profile=str(ganon_tre),
        species_abundance=str(final_abundance),
        out_path=None,          
        drop_root=True,
        drop_strain=True,
    )
    _ts("tax_profile_rebuild").add_header_to_species_abundance(
        species_abundance=str(final_abundance),
        header1="Species_TaxID",
        header2="Relative_Abundance",
    )

    # =========================================================
    # Strain Level Logic
    # =========================================================
    
    strain_group = f"{ggcat_prefix}.strain_group_abundance.tsv"
    strain_with_id = f"{ggcat_prefix}.strain_group_abundance.with_taxid_name.tsv"

    _ts("add_strain_taxid_and_name_columns").run(
        input=strain_group,
        map=str(tmp_id_table),  
        out=strain_with_id,
    )

    strain_topk_file = f"{ggcat_prefix}.strain_group_abundance.topk.tsv"

    _ts("filter_and_topk_strains_per_species").run(
        input=strain_with_id,
        out=strain_topk_file,
        min_reads=strain_min_reads,
        topk=strain_topk,
        keep_ties=True,
    )

    # Threshold calculation logic...
    species_counts = f"{ggcat_prefix}.species_counts.tsv"
    try:
        remains = 0.0; sc = 0
        if Path(species_counts).exists():
            with open(species_counts, "r") as f:
                for line in f:
                    if line.strip():
                        parts = line.split("\t")
                        if len(parts)>=2: remains += float(parts[1]); sc += 1
        
        total_reads = 1
        if Path(str(ggcat_reads)).exists():
            with open(str(ggcat_reads), "r") as rf:
                total_lines = sum(1 for _ in rf)
            total_reads = max(total_lines // 4, 1)
        
        ratio = remains / total_reads
        expr = (remains / 1_000_000.0) * (ratio ** 2) * (remains / sc) if sc > 0 else 0
        threshold = math.sqrt(expr) if expr > 0 else 0.0
    except Exception:
        threshold = 0.0

    strain_out_prefix = f"{ggcat_prefix}.predict.strain"

    _ts("strain_postprocess_global_abundance").run(
        mix=str(final_abundance), 
        input=strain_topk_file,
        threshold=threshold,
        div=strain_div,
        out_prefix=strain_out_prefix,
        emit_name=True,
    )


# =========================================================
    # Final Output Organization (Species vs. Strain)
    # =========================================================
    
    # 1. 
    raw_strain_abundance = Path(f"{strain_out_prefix}.abundance.tsv")
    final_strain_abundance = out_prefix / "strain_abundance.txt"
    final_strain_tre = out_prefix / "tax_profile_strain.tre"
    
    # 2. 
    if raw_strain_abundance.exists():
        shutil.copyfile(raw_strain_abundance, final_strain_abundance)
        
        # 3. 
        
        
        _ts("tax_profile_strain_update").run(
            tax_profile=str(ganon_tre),
            strain_abundance=str(final_strain_abundance),
            output_tre=str(final_strain_tre)
        )
    else:
        print("[PanTax-DBG] Warning: No strain abundance generated.")

    Path(temp_mix_abundance).unlink(missing_ok=True)

    # The reconstruction step above normally creates the final five-column
    # species tree. Fall back to the original tree only if reconstruction did
    # not produce an output file.
    final_species_tre = out_prefix / "tax_profile.tre"
    if not final_species_tre.exists() and ganon_tre.exists():
        shutil.copyfile(ganon_tre, final_species_tre)

    # Keep only user-facing final outputs after a successful paired-end run.
    _cleanup_intermediate_dirs(work_ganon, work_db, work_q, out_prefix / ".temp_files")

    print("[PanTax-DBG] Paired-end profiling completed.")
    print(f"  - Species profile: {final_abundance}")
    print(f"  - Species tree:    {final_species_tre}")
    print(f"  - Strain profile:  {final_strain_abundance}")
    print(f"  - Strain tree:     {final_strain_tre}")

    return str(final_abundance)
# =========================
# Single Mode
# =========================

def run_single(reads, db_prefix, out_prefix,
               report_type, ref_info, threads, k,
               species_min_abundance=0.0,
               ggcat_temp_dir=None,
               ggcat_memory=None,
               ggcat_prefer_memory=False,
               ggcat_disk_optimization_level=None,
               ggcat_intermediate_compression_level=None):
               
    if len(reads) != 1:
        raise SystemExit("[PanTax-DBG] --single mode requires exactly one -r reads file.")
    if not ref_info:
        raise SystemExit("[PanTax-DBG] --ref-info is required.")

    read = str(reads[0])
    work_ganon = out_prefix / "query_r1"
    work_db = out_prefix / "database_filter_ccDBG"
    work_q = out_prefix / "res_pantax_dbg"

    ensure_dir(work_ganon); ensure_dir(work_db); ensure_dir(work_q)

    # -------------------------------------------------------------
    # 1) Ganon Step (Single): Check existing files
    # -------------------------------------------------------------
    ganon_out_prefix = work_ganon / "results"
    ganon_tre = work_ganon / "tax_profile.tre"
    ganon_species = work_ganon / "species_abundance.txt"
    ganon_strain = work_ganon / "strain_abundance.txt"
    

    all_ganon_files_exist = (
        ganon_tre.exists() and 
        ganon_species.exists() and 
        ganon_strain.exists()
    )

    if all_ganon_files_exist:
        print("[PanTax-DBG] Found existing Ganon results (Single). Skipping classify & report steps.")
    else:
        print("[PanTax-DBG] Running Ganon classification (Single)...")
        ganon_wrapper.run_classify_single(
            db_prefix=db_prefix,
            read=read,
            out_prefix=str(ganon_out_prefix),
            threads=threads,
            report_type=report_type,
        )

        ganon_wrapper.run_report_and_postprocess(
            db_prefix=db_prefix,
            classify_prefix=str(ganon_out_prefix),
            tre_out=str(ganon_tre),
            species_out=str(ganon_species),
            strain_out=str(ganon_strain),
        )

    # Select species and strains from the reference-information table.
    tmp_tsv = work_db / "ganon_species_strain_selected.tsv"
    tmp_id_table = work_q / "tmp_id_table.tsv"

    _ts("make_ganon_pred_symlinks").run(
        ref_info=str(ref_info),
        ganon_species=str(ganon_species), 
        ganon_strain=str(ganon_strain),
        out_tsv=str(tmp_tsv),
    )
    shutil.copyfile(tmp_tsv, tmp_id_table)

    # 4) color_mapping.in
    color_map = work_db / "color_mapping.in"
    if color_map.exists():
        print(f"[PanTax-DBG] Found existing color_mapping: {color_map}. Skipping.")
    else:
        _make_color_mapping_from_mapping_tsv(str(tmp_id_table), str(color_map))

    # 5) ggcat build
    ggcat_db = work_db / "ggcatDB.fasta.lz4"
    if ggcat_db.exists():
        print(f"[PanTax-DBG] Found existing ggcatDB: {ggcat_db}. Skipping build.")
    else:
        ggcat_wrapper.run_build(
            k=k,
            threads=threads,
            color_mapping=str(color_map),
            output=str(ggcat_db),
            temp_dir=ggcat_temp_dir,
            memory=ggcat_memory,
            prefer_memory=ggcat_prefer_memory,
            disk_optimization_level=ggcat_disk_optimization_level,
            intermediate_compression_level=ggcat_intermediate_compression_level,
        )

    # Query single reads against the sample-specific ccDBG.
    ggcat_prefix = work_q / "query_ggcatDB"

    check_ggcat_outputs_single = [
        Path(f"{ggcat_prefix}.species_counts.tsv"),
        
    ]
    if all(p.exists() for p in check_ggcat_outputs_single):
        print(f"[PanTax-DBG] Found existing GGCAT query results (Single). Skipping query.")
    else:
        ggcat_wrapper.run_query(
            db=str(ggcat_db),
            reads=[read],
            k=k,
            threads=threads,
            out_prefix=str(ggcat_prefix),
            single=True,
        )

    # Apply the adaptive threshold and length correction without integration.
    temp_abundance = out_prefix / "species_abundance.raw.txt"
    _run_threshold_and_lc_for_single(
        ggcat_prefix=str(ggcat_prefix),
        tmp_id_table=str(tmp_id_table),
        reads_path=read,
        output=str(temp_abundance),
    )
    
    # [FILTER] Final Result Filtering & Renormalization
    final_abundance = out_prefix / "species_abundance.txt"
    _filter_and_renormalize_abundance(
        input_path=str(temp_abundance),
        output_path=str(final_abundance),
        min_val=species_min_abundance
    )

    # 8.1) rebuild tax_profile
    _ts("tax_profile_rebuild").rebuild_tax_profile_with_species_abundance(
        tax_profile=str(ganon_tre),
        species_abundance=str(final_abundance),
        out_path=None,
        drop_root=True,
        drop_strain=True,
    )
    # 8.2) add header
    _ts("tax_profile_rebuild").add_header_to_species_abundance(
        species_abundance=str(final_abundance),
        header1="Taxonomic_ID",
        header2="Relative_Abundance",
    )
    
    Path(temp_abundance).unlink(missing_ok=True)

    print(f"[PanTax-DBG] Single-end profiling completed: {final_abundance}")
    return str(final_abundance)

# =========================
# Utils
# =========================

def _count_non_empty_lines(path):
    c = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                c += 1
    return c

def _prepare_ggcat_reads_paired(r1, r2, out_fastq, threads=8):
    # Write QC reports to the final output directory rather than the current
    # shell working directory.
    final_output_dir = Path(out_fastq).parent.parent
    fastp_html = final_output_dir / "fastp.html"
    fastp_json = final_output_dir / "fastp.json"

    cmd = [
        "fastp",
        "-i", r1,
        "-I", r2,
        "--thread", str(threads),
        "--html", str(fastp_html),
        "--json", str(fastp_json),
        "--stdout",
    ]
    with open(out_fastq, "w") as out_f:
        try:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            raise SystemExit("[PanTax-DBG][error] fastp not found.")
        except subprocess.CalledProcessError as e:
            raise SystemExit(f"[PanTax-DBG][error] fastp failed with exit code {e.returncode}")

def _make_color_mapping_from_mapping_tsv(mapping_tsv, out_path):
    with open(mapping_tsv, "r", encoding="utf-8") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        header = next(fin, None)
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            new_strain_taxid = parts[3]
            fasta_path = parts[4]
            if new_strain_taxid and fasta_path:
                fout.write(f"{new_strain_taxid}\t{fasta_path}\n")

def _run_threshold_and_mix_for_paired(
    ggcat_prefix,
    tmp_id_table,
    ganon_predict_spy,
    tre_file,
    reads_path,
    output,
    strain_min_reads=5.0,
    strain_topk=5,
    strain_div=5.0,
):
    species_counts = f"{ggcat_prefix}.species_counts.tsv"

    remains = 0.0
    species = 0
    with open(species_counts, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            remains += v
            species += 1

    with open(reads_path, "r", encoding="utf-8") as rf:
        total_lines = sum(1 for _ in rf)
    total_reads = max(total_lines // 4, 1)

    if species == 0:
        print("[PanTax-DBG] No species in ggcat species_counts.tsv (paired).")
        with open(output, "w") as f: f.write("speciesID\tabundance\n")
        return

    ratio = remains / total_reads
    expr = (remains / 1_000_000.0) * (ratio ** 2) * (remains / species)
    threshold = math.sqrt(expr) if expr > 0 else 0.0

    filtered_counts = f"{ggcat_prefix}.species_counts_more{threshold:.6g}.tsv"
    with open(species_counts, "r", encoding="utf-8", errors="ignore") as fin, \
         open(filtered_counts, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                v = float(parts[1])
            except ValueError:
                continue
            if v > threshold:
                fout.write(line)

    abund_more = f"{ggcat_prefix}.species_abundance_more{threshold:.6g}.tsv"
    _ts("length_corrected_abundance").run(
        counts_file=filtered_counts,
        mapping_file=tmp_id_table,
        output_file=abund_more,
    )

    dbg_tmp = f"{abund_more}.tmp"
    dbg_vals = []  # Collect the leading GGCAT abundance estimates.
    with open(abund_more, "r", encoding="utf-8", errors="ignore") as fin, \
         open(dbg_tmp, "w", encoding="utf-8") as fout:
        header = fin.readline()
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            fout.write(f"{parts[0]}\t{parts[5]}\n")
            try:
                dbg_vals.append(float(parts[5]))
            except ValueError:
                continue

    # Obtain the first-stage profile on the same graph-retained species set.
    raw_mix = f"{ggcat_prefix}.raw_filtered_abundance_prediction.tsv"
    ganon_vals = [] 
    _ts("mix_predictions").run(
        dbg_file=dbg_tmp,
        ganon_file=ganon_predict_spy,
        weight=0.0,
        output=raw_mix,
    )
    
    with open(raw_mix, "r") as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 2:
                try: ganon_vals.append(float(p[1]))
                except ValueError: continue

    judge = _ts("judge_strategy")

    # Compare the concentration of the leading abundance estimates.
    top1_ggcat, top2_ggcat = judge.get_top_two_abundance_from_list(dbg_vals)
    top1_ganon, top2_ganon = judge.get_top_two_abundance_from_list(ganon_vals)
    
    strategy = judge.calculate_decision(top1_ganon, top2_ganon, top1_ggcat, top2_ggcat)

    # Scale the integration weight by the unclassified fraction when the
    # first-stage abundance allocation is selected.
    def tre_weight_calc(tre_file):
        w = 0.0; found = False
        try:
            with open(tre_file, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if parts and parts[0] == "unclassified" and len(parts) >= 9:
                        try: val = float(parts[8])
                        except ValueError: val = 0.0
                        w = val / 100.0 / 5.0
                        found = True; break
        except FileNotFoundError: return 0.0
        return max(0.0, min(1.0, w)) if found else 0.0

    tre_w = tre_weight_calc(tre_file)

    if strategy == "GGCAT":
        weight = 1.0
    else:
        weight = tre_w
        
    print(f"[PanTax-DBG] Mixing Strategy: {strategy}, Final Weight: {weight:.6f} (i-logic applied)")

    # Generate the final integrated abundance profile.
    final_tmp = f"{ggcat_prefix}.mix_abundance_prediction.tsv"
    _ts("mix_predictions").run(
        dbg_file=dbg_tmp,
        ganon_file=ganon_predict_spy,
        weight=weight,
        output=final_tmp,
    )

    with open(final_tmp, 'r', encoding='utf-8') as fin, \
        open(output, 'w', encoding='utf-8') as fout:
        fout.write("speciesID\tabundance\n")
        shutil.copyfileobj(fin, fout)

def _run_threshold_and_lc_for_single(
    ggcat_prefix,
    tmp_id_table,
    reads_path,
    output,
):
    species_counts = f"{ggcat_prefix}.species_counts.tsv"

    remains = 0.0
    species = 0
    with open(species_counts, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2: continue
            try:
                v = float(parts[1])
            except ValueError: continue
            remains += v
            species += 1

    with open(reads_path, "r", encoding="utf-8") as rf:
        total_lines = sum(1 for _ in rf)
    total_reads = max(total_lines // 4, 1)

    if species == 0:
         with open(output, "w") as f: f.write("speciesID\tabundance\n")
         return

    ratio = remains / total_reads
    expr = (remains / 1_000_000.0) * (ratio ** 2) * (remains / species)
    threshold = math.sqrt(expr) if expr > 0 else 0.0

    filtered_counts = f"{species_counts}.more{threshold:.6g}.tsv"
    with open(species_counts, "r", encoding="utf-8") as fin, \
         open(filtered_counts, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2: continue
            try:
                v = float(parts[1])
            except ValueError: continue
            if v > threshold:
                fout.write(line)

    abund_more = f"{species_counts}.abundance_more{threshold:.6g}.tsv"
    _ts("length_corrected_abundance").run(
        counts_file=filtered_counts,
        mapping_file=tmp_id_table,
        output_file=abund_more,
    )

    pairs = []
    with open(abund_more, "r", encoding="utf-8", errors="ignore") as fin:
        _ = fin.readline()
        for line in fin:
            if not line.strip(): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                sid = parts[0].strip()
                try:
                    rel = float(parts[5])
                except ValueError: continue
                if sid:
                    pairs.append((sid, rel))

    pairs.sort(key=lambda x: x[1], reverse=True)

    with open(output, "w", encoding="utf-8") as fout:
        fout.write("speciesID\tabundance\n")
        for sid, rel in pairs:
            fout.write(f"{sid}\t{rel}\n")

# =========================
# Module CLI wrapper
# =========================

def cli():
    p = argparse.ArgumentParser("pantax_dbg-profile")
    p.add_argument("-r", "--reads", action="append", required=True, metavar="FASTQ",
                   help="Input paired-end FASTQ file. Supply exactly twice: -r read1.fq.gz -r read2.fq.gz.")
    p.add_argument("--db-prefix", required=True, metavar="PREFIX", help="PanTax-DBG database prefix.")
    p.add_argument("--ref-info", required=True, metavar="FILE", help="Reference metadata TSV.")
    p.add_argument("--out", required=True, metavar="DIR", help="Output directory.")
    p.add_argument("--threads", type=int, default=8, metavar="INT", help="Number of threads. (default: 8)")
    p.add_argument("-k", "--kmer", type=int, default=31, metavar="INT", help="k-mer size. (default: 31)")

    args = p.parse_args()
    if len(args.reads) != 2:
        raise SystemExit(
            "[PanTax-DBG][error] Expected two paired-end short-read FASTQ files. "
            "Please provide: -r read1.fq.gz -r read2.fq.gz"
        )

    run(
        reads=args.reads,
        single=False,
        db_prefix=args.db_prefix,
        out_prefix=args.out,
        report_type="abundance",
        file_info=args.ref_info,
        threads=args.threads,
        k=args.kmer,
    )


if __name__ == "__main__":
    cli()
