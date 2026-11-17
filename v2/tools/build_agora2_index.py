"""Build the AGORA2 species -> {KEGG compound: [EC]} index ONCE, outside the pipeline.

    python tools/build_agora2_index.py            # refuses if the index already exists
    python tools/build_agora2_index.py --force    # rebuild

Reads the 1,439 AGORA2 reconstructions (1.8 GB) under Data/micom/agora2_extracted and
writes utils.AGORA2_INDEX_CACHE plus a sidecar (utils.AGORA2_INDEX_SIDECAR) recording
what it was built from. Notebooks only ever load the result (utils.load_agora2_ec_index);
the scan itself is what overloaded the old analysis server, so it never runs inside a
notebook. checks.py D5 verifies the sidecar against the reconstructions on disk.
"""
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import utils  # noqa: E402


def main() -> int:
    cache, sidecar = utils.AGORA2_INDEX_CACHE, utils.AGORA2_INDEX_SIDECAR
    if cache.exists() and "--force" not in sys.argv:
        print(f"{cache} exists; pass --force to rebuild")
        return 1
    src = utils.find_agora2_data_dir()
    fp = utils.agora2_source_fingerprint(src)
    print(f"scanning {fp['n_files']} reconstructions in {src} ...")
    t0 = time.time()
    index = utils.build_agora2_ec_index(src, cache_path=cache, force_rebuild=True)
    secs = time.time() - t0
    n_links = sum(len(v) for v in index.values())
    sidecar.write_text(json.dumps({
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_dir": str(src.relative_to(utils.DATA_DIR)),
        **fp,
        "n_species_indexed": len(index),
        "n_species_compound_links": n_links,
        "seconds": round(secs, 1),
    }, indent=2), encoding="utf-8")
    print(f"indexed {len(index)} species ({fp['n_files'] - len(index)} without a "
          f"KEGG-annotated producing reaction), {n_links} species-compound links, "
          f"{secs:.0f} s\nwrote {cache}\nwrote {sidecar}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
