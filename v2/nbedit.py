"""Deterministic notebook cell editor used to maintain the v2.1 notebooks.

Why this exists: several v2 notebooks are too large for an interactive
cell-by-cell editor, and the previous approach -- regex/substring patch scripts --
broke repeatedly on whitespace drift and was not idempotent.

This editor does not search for text. Every operation addresses a cell by its
`id` and replaces its source wholesale, so:
  * applying the same edit twice is a no-op (set is idempotent by construction),
  * a missing id is a hard error, never a silent skip,
  * nothing depends on the current content of the cell being replaced.

Operations
----------
set(id, source)                  replace that cell's source
insert_after(id, new_id, type, source)
                                 insert once; if new_id exists, set it instead
delete(id)                       remove that cell (no-op if already absent)

Usage from an edit script:

    from nbedit import Notebook
    nb = Notebook("02_analysis.ipynb")
    nb.set("cell-001", '''...''')
    nb.insert_after("cell-010", "cell-010b", "code", '''...''')
    nb.save()
"""

import json
from pathlib import Path


class Notebook:
    def __init__(self, path):
        self.path = Path(path)
        self.doc = json.loads(self.path.read_text(encoding="utf-8"))
        self.changes = []

    # -- internals ---------------------------------------------------------
    @property
    def cells(self):
        return self.doc["cells"]

    def _index(self, cell_id):
        for i, c in enumerate(self.cells):
            if c.get("id") == cell_id:
                return i
        raise KeyError(
            f"{self.path.name}: no cell with id {cell_id!r}. "
            f"Existing ids: {[c.get('id') for c in self.cells]}"
        )

    def has(self, cell_id):
        return any(c.get("id") == cell_id for c in self.cells)

    @staticmethod
    def _as_source(text):
        """Store source as a list of lines WITH trailing newlines, which is the
        format jupyter itself writes; a single string also loads fine but makes
        git diffs useless."""
        text = text.rstrip("\n")
        lines = text.split("\n")
        return [l + "\n" for l in lines[:-1]] + [lines[-1]]

    # -- operations --------------------------------------------------------
    def set(self, cell_id, source):
        i = self._index(cell_id)
        self.cells[i]["source"] = self._as_source(source)
        if self.cells[i].get("cell_type") == "code":
            self.cells[i]["outputs"] = []
            self.cells[i]["execution_count"] = None
        self.changes.append(f"set {cell_id}")
        return self

    def set_type(self, cell_id, cell_type):
        i = self._index(cell_id)
        c = self.cells[i]
        if c.get("cell_type") == cell_type:
            return self
        c["cell_type"] = cell_type
        if cell_type == "code":
            c.setdefault("outputs", [])
            c.setdefault("execution_count", None)
            c.setdefault("metadata", {})
        else:
            c.pop("outputs", None)
            c.pop("execution_count", None)
        self.changes.append(f"set_type {cell_id} -> {cell_type}")
        return self

    def insert_after(self, cell_id, new_id, cell_type, source):
        if self.has(new_id):                    # idempotent re-run
            return self.set(new_id, source)
        i = self._index(cell_id)
        cell = {"cell_type": cell_type, "id": new_id, "metadata": {},
                "source": self._as_source(source)}
        if cell_type == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
        self.cells.insert(i + 1, cell)
        self.changes.append(f"insert {new_id} after {cell_id}")
        return self

    def delete(self, cell_id):
        if not self.has(cell_id):
            return self
        del self.cells[self._index(cell_id)]
        self.changes.append(f"delete {cell_id}")
        return self

    @classmethod
    def create(cls, path, cells, kernel="python3"):
        """Build a notebook from scratch out of (cell_type, id, source) triples.

        Used where a notebook's entire subject changed and editing cell-by-cell
        would be more error-prone than restating it.
        """
        path = Path(path)
        doc = {
            "cells": [],
            "metadata": {
                "kernelspec": {"display_name": "Python 3", "language": "python",
                               "name": kernel},
                "language_info": {"name": "python"},
            },
            "nbformat": 4,
            "nbformat_minor": 5,
        }
        obj = cls.__new__(cls)
        obj.path = path
        obj.doc = doc
        obj.changes = []
        for cell_type, cell_id, source in cells:
            cell = {"cell_type": cell_type, "id": cell_id, "metadata": {},
                    "source": cls._as_source(source)}
            if cell_type == "code":
                cell["outputs"] = []
                cell["execution_count"] = None
            doc["cells"].append(cell)
            obj.changes.append(f"create {cell_id}")
        return obj

    # -- output ------------------------------------------------------------
    def save(self, verbose=True):
        self.path.write_text(
            json.dumps(self.doc, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        if verbose:
            print(f"{self.path.name}: {len(self.changes)} change(s)")
            for ch in self.changes:
                print(f"  {ch}")
        return self
