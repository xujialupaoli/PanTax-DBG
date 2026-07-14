# pantax_dbg/utils.py
import os
import subprocess
import sys
from pathlib import Path

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

# def run_cmd(cmd, log_prefix="[PanTax-DBG] ", *, echo=True, silence=False):
#     if echo and log_prefix:
#         print(log_prefix + " " + " ".join(cmd), file=sys.stderr)
#     stdout = subprocess.DEVNULL if silence else None
    # stderr = subprocess.DEVNULL if silence else None
    # subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)


def _env_with_internal_backends(env=None):
    run_env = dict(os.environ if env is None else env)
    try:
        from .paths import get_pantax_dbg_libexec
        libexec = str(get_pantax_dbg_libexec())
        current = run_env.get("PATH", "")
        run_env["PATH"] = libexec + (os.pathsep + current if current else "")
    except Exception:
        pass
    return run_env


def run_cmd(cmd, log_prefix="[PanTax-DBG]", *, echo=True, silence=False, env=None):
    if echo and log_prefix:
        print(log_prefix + " " + " ".join(map(str, cmd)), file=sys.stderr)

    run_env = _env_with_internal_backends(env)

    if silence:
        result = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=run_env,
        )
        if result.returncode != 0:
            print(result.stderr[-4000:], file=sys.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd)
    else:
        subprocess.run(cmd, check=True, env=run_env)
