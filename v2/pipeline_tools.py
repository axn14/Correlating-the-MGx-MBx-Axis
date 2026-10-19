"""Helpers for run_pipeline.sh (stdlib + psutil only).

  python pipeline_tools.py memrun <peak_file> <cmd> [args...]
      Run <cmd>, write the peak total RSS (MB) of it and all child processes
      (e.g. the Jupyter kernel) to <peak_file>; exit with <cmd>'s exit code.

  python pipeline_tools.py scan <notebook.ipynb>
      Print output lines of an executed notebook that look like problems
      (tracebacks, WARNING, SKIPPED, '!!!'). nbconvert's own log does not
      contain cell output, so this is the only place such messages surface.
"""
import json
import re
import subprocess
import sys
import time

PATTERN = re.compile(r"Traceback|WARNING|SKIPPED|SKIP \(|!!!|Error:", re.IGNORECASE)
WARN_CLASS = re.compile(r"\b([A-Z][A-Za-z]*(?:Warning|Error))\b")


def memrun(peak_file, cmd):
    import psutil
    proc = subprocess.Popen(cmd)
    ps = psutil.Process(proc.pid)
    peak = 0.0
    while proc.poll() is None:
        try:
            procs = [ps] + ps.children(recursive=True)
            peak = max(peak, sum(p.memory_info().rss for p in procs if p.is_running()) / 1e6)
        except psutil.Error:
            pass
        time.sleep(0.5)
    with open(peak_file, "w") as f:
        f.write(f"{peak:.0f}\n")
    return proc.returncode


def scan(path):
    """Group problem lines per (cell, warning class) and count them. Nothing is
    truncated: the previous version kept the first 40 distinct lines, and NB09's 37
    non-converged GLASSO permutation fits filled that budget, so 3 of them -- and
    anything later -- never reached the summary."""
    nb = json.load(open(path, encoding="utf-8"))
    groups, order = {}, []
    for c in nb.get("cells", []):
        if c.get("cell_type") != "code":
            continue
        cid = c.get("id") or "?"
        for o in c.get("outputs", []):
            if o.get("output_type") == "error":
                items = [(f"ERROR {o.get('ename')}", f"{o.get('ename')}: {o.get('evalue')}")]
            elif o.get("output_type") == "stream":
                items = []
                for line in "".join(o.get("text", "")).splitlines():
                    s = line.strip()
                    # The second line of every printed warning is a source echo
                    # ("warnings.warn("), which carries no information.
                    if not PATTERN.search(line) or "FutureWarning" in line \
                            or s.startswith("warnings.warn("):
                        continue
                    m = WARN_CLASS.search(line)
                    items.append((m.group(1) if m else s[:120], s))
            else:
                continue
            for key, example in items:
                k = (cid, key)
                if k not in groups:
                    groups[k] = [0, example]
                    order.append(k)
                groups[k][0] += 1
    out = []
    for cid, key in order:
        n, example = groups[(cid, key)]
        out.append(f"[{cid}] {key} x{n}: {example[:180]}" if n > 1
                   else f"[{cid}] {example[:220]}")
    print("\n".join(out) if out else "(no warnings)")


if __name__ == "__main__":
    if sys.argv[1] == "memrun":
        sys.exit(memrun(sys.argv[2], sys.argv[3:]))
    elif sys.argv[1] == "scan":
        scan(sys.argv[2])
    else:
        sys.exit(f"unknown command {sys.argv[1]}")
