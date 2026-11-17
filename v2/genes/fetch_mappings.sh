#!/bin/bash
# Fetch the UniRef90 -> EC and UniRef90 -> KO maps that turn curatedMetagenomicData's
# stratified gene families into the enzyme activities the chapter makes claims about.
#
# Run on Dardel:  bash fetch_mappings.sh
#
# ---------------------------------------------------------------------------
# Why Zenodo and not huttenhower.sph.harvard.edu
#
# HUMAnN's own downloader takes the whole utility_mapping bundle:
#   https://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz
# which is 2,734,293,864 B -- 2.73 GB, measured by HEAD on 2026-11-08 -- to deliver
# the 26 MB file we need. It is also unversioned beyond a date string, served over a
# redirect chain ending at a Globus endpoint, and carries no checksum.
#
# Zenodo record 17100034 ("ChocoPhlAn mapping files", doi:10.5281/zenodo.17100034)
# serves the same files individually, each with an md5, under a DOI that can be cited
# in a methods section and that resolves to a fixed version. 45 MB instead of 2.73 GB,
# and a provenance statement that does not depend on a university web server still
# being up at viva time.
#
# Note the path trap: .../humann_data/full_mapping/<file> is a 404, while
# .../humann_data/full_mapping_v201901b.tar.gz (no subdirectory) is a 200. Several
# tutorials give the first form. Neither matters here, but it is why this script does
# not fall back to that host.
# ---------------------------------------------------------------------------
set -euo pipefail

DEST="${1:-/cfs/klemming/projects/supr/naiss2024-6-169/andnair/mappings}"
BASE="https://zenodo.org/records/17100034/files"

# name -> md5, from the record's own API listing on 2026-11-08. Checked after every
# download: a truncated transfer of a gzip file does not necessarily fail to open, and
# a short map silently drops enzymes, which would read as biological absence.
declare -A WANT=(
  [map_level4ec_uniref90.txt.gz]=3ace48af4aed6a2334959ed64d1f76e8
  [map_ko_uniref90.txt.gz]=44c7c7b6ce27bd4600f2e99de3c0aa43
  [map_ec_name.txt.gz]=d324852a5092d38bd3a7035f3506e8f1
)
declare -A SIZE=(
  [map_level4ec_uniref90.txt.gz]=27591865
  [map_ko_uniref90.txt.gz]=17547123
  [map_ec_name.txt.gz]=75265
)

mkdir -p "$DEST"
cd "$DEST"

fail=0
for f in "${!WANT[@]}"; do
    if [ -f "$f" ]; then
        have=$(md5sum "$f" | awk '{print $1}')
        if [ "$have" = "${WANT[$f]}" ]; then
            printf '  %-32s already present, md5 ok\n' "$f"
            continue
        fi
        echo "  $f present but md5 is $have, expected ${WANT[$f]} -- refetching" >&2
        rm -f "$f"
    fi

    printf '  %-32s fetching %s B ...\n' "$f" "${SIZE[$f]}"
    # --fail so an HTML error page is never written over the target as if it were data.
    curl -fL --retry 3 --retry-delay 5 --max-time 1800 -o "$f.part" "$BASE/$f?download=1"
    mv "$f.part" "$f"

    have=$(md5sum "$f" | awk '{print $1}')
    if [ "$have" != "${WANT[$f]}" ]; then
        echo "  MD5 MISMATCH for $f: got $have, expected ${WANT[$f]}" >&2
        fail=1
    else
        printf '  %-32s md5 ok\n' "$f"
    fi
done

[ "$fail" -eq 0 ] || { echo "at least one file failed verification" >&2; exit 1; }

# Record what was fetched, so the methods section does not have to be reconstructed
# from shell history.
cat > provenance.json <<JSON
{
  "fetched_utc": "$(date -u '+%Y-%m-%dT%H:%M:%SZ')",
  "source": "Zenodo record 17100034",
  "doi": "10.5281/zenodo.17100034",
  "title": "ChocoPhlAn mapping files",
  "rationale": "individual checksummed files under a citable DOI, 45 MB, instead of HUMAnN's 2.73 GB full_mapping_v201901b.tar.gz bundle which carries no checksum",
  "files": {
$(for f in "${!WANT[@]}"; do
    printf '    "%s": {"md5": "%s", "bytes": %s},\n' "$f" "${WANT[$f]}" "$(stat -c %s "$f")"
  done | sed '$ s/,$//')
  }
}
JSON

echo
echo "verified in $DEST:"
ls -la "$DEST"
echo
echo "shape of the EC map (one line per EC, then tab-separated UniRef90 ids):"
# `set -o pipefail` plus `head` closing the pipe makes zcat die of SIGPIPE and the
# script exit 141 having succeeded, so pipefail comes off for the reporting below.
set +o pipefail
zcat map_level4ec_uniref90.txt.gz | head -2 | cut -c1-160
echo "  total EC lines: $(zcat map_level4ec_uniref90.txt.gz | wc -l)"
