# pantax_dbg/ggcat_wrapper.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import tempfile
from pathlib import Path

from .utils import run_cmd, ensure_dir
from .paths import get_dbg_ggcat


def _ggcat_bin():
    return os.environ.get("PANTAX_DBG_GGCAT_BIN") or get_dbg_ggcat()


def _default_tmp_root() -> Path:
    # Prefer user/job TMPDIR when present; otherwise use /tmp.
    return Path(os.environ.get("TMPDIR") or "/tmp")


def run_build(
    k,
    threads,
    color_mapping,
    output,
    temp_dir=None,
    memory=None,
    prefer_memory=False,
    disk_optimization_level=None,
    intermediate_compression_level=None,
):
    """
    Build colored DBG database with bundled DBG-ggcat.

    By default, PanTax-DBG sends DBG-ggcat intermediate files to a unique
    per-run directory under ${TMPDIR:-/tmp}. This reduces network-filesystem I/O
    during benchmarks and large cohort profiling.

    User-facing option:
      --ggcat-temp-dir <DIR>

    Advanced resource knobs are intentionally kept internal for now:
      memory, prefer_memory, disk_optimization_level, intermediate_compression_level
    """
    ensure_dir(Path(output).parent)

    auto_temp_dir = False
    if temp_dir:
        ggcat_temp_dir = Path(temp_dir)
        ensure_dir(ggcat_temp_dir)
    else:
        tmp_root = _default_tmp_root()
        ensure_dir(tmp_root)
        ggcat_temp_dir = Path(tempfile.mkdtemp(prefix="pantax-dbg-ggcat-", dir=str(tmp_root)))
        auto_temp_dir = True

    print(f"[PanTax-DBG][DBG-ggcat build] temporary directory: {ggcat_temp_dir}")

    cmd = [
        _ggcat_bin(), "build",
        "-k", str(k),
        "-j", str(threads),
        "-c",
        "-d", color_mapping,
        "-s", "1",
        "-o", output,
        "-t", str(ggcat_temp_dir),
    ]

    # Internal-only advanced knobs. They are not exposed in the default CLI.
    if memory is not None:
        cmd += ["-m", str(memory)]
    if prefer_memory:
        cmd.append("-p")
    if disk_optimization_level is not None:
        cmd += ["--disk-optimization-level", str(disk_optimization_level)]
    if intermediate_compression_level is not None:
        cmd += ["--intermediate-compression-level", str(intermediate_compression_level)]

    try:
        run_cmd(cmd, "[PanTax-DBG][DBG-ggcat build]", echo=True)
    except Exception:
        if auto_temp_dir:
            print(
                f"[PanTax-DBG][warning] DBG-ggcat build failed; keeping temporary directory for debugging: {ggcat_temp_dir}"
            )
        raise
    else:
        if auto_temp_dir:
            shutil.rmtree(ggcat_temp_dir, ignore_errors=True)
            print(f"[PanTax-DBG][DBG-ggcat build] removed temporary directory: {ggcat_temp_dir}")


def run_query(db, reads, k, threads, out_prefix, single=False):
    ensure_dir(Path(out_prefix).parent)
    cmd = [
        _ggcat_bin(), "query",
        "--colors",
        "-k", str(k),
        "-j", str(threads),
        db,
    ]
    cmd += list(reads)
    cmd += [
        "--colored-query-output-format", "JsonLinesWithNames",
        "-o", out_prefix,
    ]
    if single:
        cmd.append("--single")

    run_cmd(cmd, "[PanTax-DBG][DBG-ggcat query]", echo=True)
