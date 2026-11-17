"""Apply the v2.0 wiki pages to a local checkout of the GitHub wiki.

    python docs/wiki/sync_wiki.py <path to the .wiki checkout>

Copies the v2 pages, renames the v1 Key-Results page to Key-Results-v1-thesis, and puts a
one-line notice at the top of every page that documents v1.0-thesis. It only edits files;
committing and pushing the wiki is a separate, deliberate step.
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2_PAGES = ["Home.md", "What-changed-in-v2.0.md", "Evidence-classes.md", "Key-results.md"]
V1_PAGES = ["Pipeline-Overview.md", "Methods.md", "Evidence-Streams.md", "Data-Sources.md",
            "Data.md", "Installation.md", "Results.md", "Key-Results-v1-thesis.md"]
NOTICE = ("> This page documents the earlier pipeline (tag `v1.0-thesis`). Its results are "
          "superseded; see [[What changed in v2.0]].\n\n")


def main(wiki: Path) -> None:
    old = wiki / "Key-Results.md"
    if old.exists():
        old.rename(wiki / "Key-Results-v1-thesis.md")
    for p in V2_PAGES:
        shutil.copy(HERE / p, wiki / p)
    for p in V1_PAGES:
        f = wiki / p
        if f.exists():
            text = f.read_text(encoding="utf-8")
            if not text.startswith("> This page documents the earlier pipeline"):
                f.write_text(NOTICE + text, encoding="utf-8")
    print("updated", wiki)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
