# themis/utils.py
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


def run_cmd(cmd, log_prefix="[PanTax-DBG]", *, echo=True, silence=False):
    if echo and log_prefix:
        print(log_prefix + " " + " ".join(map(str, cmd)), file=sys.stderr)

    if silence:
        result = subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode != 0:
            print(result.stderr[-4000:], file=sys.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd)
    else:
        subprocess.run(cmd, check=True)