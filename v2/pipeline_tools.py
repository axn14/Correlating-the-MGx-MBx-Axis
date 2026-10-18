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
    nb = json.load(open(path, encoding="utf-8"))
    hits = []
    for c in nb.get("cells", []):
        if c.get("cell_type") != "code":
            continue
        for o in c.get("outputs", []):
            if o.get("output_type") == "error":
                hits.append(f"[{c.get('id')}] ERROR {o.get('ename')}: {o.get('evalue')}")
                continue
            text = "".join(o.get("text", "")) if o.get("output_type") == "stream" else ""
            for line in text.splitlines():
                if PATTERN.search(line) and "FutureWarning" not in line:
                    hits.append(f"[{c.get('id')}] {line.strip()[:200]}")
    seen, out = set(), []
    for h in hits:
        if h not in seen:
            seen.add(h)
            out.append(h)
    print("\n".join(out[:40]) if out else "(no warnings)")


if __name__ == "__main__":
    if sys.argv[1] == "memrun":
        sys.exit(memrun(sys.argv[2], sys.argv[3:]))
    elif sys.argv[1] == "scan":
        scan(sys.argv[2])
    else:
        sys.exit(f"unknown command {sys.argv[1]}")
