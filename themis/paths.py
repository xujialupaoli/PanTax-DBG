from pathlib import Path
import os
import sys


def get_pantax_dbg_libexec() -> Path:
    env_override = os.environ.get("PANTAX_DBG_LIBEXEC")
    if env_override:
        return Path(env_override)

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        return Path(conda_prefix) / "libexec" / "pantax-dbg"

    return Path(sys.prefix) / "libexec" / "pantax-dbg"


def get_internal_binary(name: str) -> str:
    path = get_pantax_dbg_libexec() / name
    if not path.exists():
        raise FileNotFoundError(
            f"Cannot find internal PanTax-DBG binary: {path}\n"
            "PanTax-DBG requires bundled modified backends. "
            "Please check your installation or set PANTAX_DBG_LIBEXEC."
        )
    return str(path)


def get_dbg_ganon() -> str:
    return get_internal_binary("ganon")


def get_dbg_ggcat() -> str:
    return get_internal_binary("dbg-ggcat")
