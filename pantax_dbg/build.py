# pantax_dbg/build.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import os
import sys
from pathlib import Path
from .utils import run_cmd
from .paths import get_dbg_ganon


def _ganon_bin():
    return os.environ.get("PANTAX_DBG_GANON_BIN") or get_dbg_ganon()


def build_custom(args):
    cmd = [_ganon_bin(), "build-custom"] + list(args)
    run_cmd(cmd, "[PanTax-DBG][build-custom]")


def main():
    if len(sys.argv) < 2:
        print("Usage: pantax-dbg build-custom [ganon build-custom options...]", file=sys.stderr)
        sys.exit(1)
    cmd = [_ganon_bin(), "build-custom"] + sys.argv[1:]
    run_cmd(cmd, "[PanTax-DBG][build-custom]")


if __name__ == "__main__":
    main()
