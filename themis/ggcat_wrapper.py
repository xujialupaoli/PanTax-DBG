# themis/ggcat_wrapper.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path
from .utils import run_cmd, ensure_dir
from .paths import get_dbg_ggcat


# GGCAT_BIN = os.environ.get("PANTAX_DBG_GGCAT_BIN") or get_dbg_ggcat()
def _ggcat_bin():
    return os.environ.get("PANTAX_DBG_GGCAT_BIN") or get_dbg_ggcat()

def run_build(k, threads, color_mapping, output):
    """
    
      ggcat build -k K -j threads -c -d color_mapping.in -s 1 -o ggcatDB.fasta.lz4
    """
    ensure_dir(Path(output).parent)
    cmd = [
        _ggcat_bin(), "build",
        "-k", str(k),
        "-j", str(threads),
        "-c",
        "-d", color_mapping,
        "-s", "1",
        "-o", output,
    ]
    run_cmd(cmd, echo=False)


def run_query(db, reads, k, threads, out_prefix, single=False):
    """
    
      ggcat query --colors -k K -j threads db reads... \
        --colored-query-output-format JsonLinesWithNames -o out_prefix [--single]
    """
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

    run_cmd(cmd, echo=False)
