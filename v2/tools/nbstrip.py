#!/usr/bin/env python3
"""git clean filter for notebooks: what is committed is source only.

Working copies may carry outputs -- from an interactive session in the editor, for
example -- but the committed form never does, so git diffs show source changes and
nothing else. Executed copies of every pipeline run live in Results/runs/<run_id>/.

Canonical form, applied to every notebook:
  - code cells: outputs emptied, execution_count null;
  - per-cell volatile metadata dropped: nbclient's "execution" timings, "ExecuteTime";
  - kernelspec.display_name fixed to "Python 3" (the editor writes the venv name);
  - language_info.version dropped (it records whichever interpreter ran it last);
  - every source stored as nbformat's list of lines.

Invalid JSON exits non-zero, so a half-written notebook cannot be committed silently.

Setup (repo-local; .gitattributes already maps *.ipynb to this filter):
    git config filter.nbstrip.clean "python3 v2/tools/nbstrip.py"
    git config filter.nbstrip.smudge cat
    git config filter.nbstrip.required true

Stdlib only, so it runs under any Python 3 that git invokes.
"""
import json
import sys

VOLATILE_CELL_META = ("execution", "ExecuteTime")


def canonicalise(nb: dict) -> dict:
    for cell in nb.get("cells", []):
        src = cell.get("source", "")
        text = "".join(src) if isinstance(src, list) else src
        cell["source"] = text.splitlines(keepends=True)
        meta = cell.get("metadata")
        if isinstance(meta, dict):
            for key in VOLATILE_CELL_META:
                meta.pop(key, None)
        if cell.get("cell_type") == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
    meta = nb.setdefault("metadata", {})
    if isinstance(meta.get("kernelspec"), dict):
        meta["kernelspec"]["display_name"] = "Python 3"
    if isinstance(meta.get("language_info"), dict):
        meta["language_info"].pop("version", None)
    return nb


def main() -> int:
    raw = sys.stdin.buffer.read().decode("utf-8")
    try:
        nb = json.loads(raw)
    except json.JSONDecodeError as exc:
        sys.stderr.write(f"nbstrip: not valid notebook JSON ({exc}); refusing to stage it\n")
        return 1
    out = json.dumps(canonicalise(nb), indent=1, ensure_ascii=False) + "\n"
    sys.stdout.buffer.write(out.encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
