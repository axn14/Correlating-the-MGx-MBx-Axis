"""The gene layer, reformulated (Phase D of the 2026-11-11 plan; CONTRACTS 5d class G).

Three things the original `progression_axis.py` analysis could not claim, fixed here.

1. **"Capacity is conserved" was stated on a relative scale.** HUMAnN gene families are
   relative abundances (rows sum to 1), so a flat EC total says nothing about absolute
   capacity per gram of stool -- biomass, bleeding or diarrhoea can move that by orders of
   magnitude while the shares stay put. There is no load or qPCR measurement in Yachida,
   so absolute capacity is NOT estimable, and nothing here pretends otherwise. What can be
   done is to make the relative statement compositionally coherent: each EC total is
   expressed as an **additive log-ratio (ALR) to a fixed reference subcomposition** -- the
   6,336 prevalence-matched random background gene families in `gf_background_long` --
   which is scale-invariant, unlike a raw share.

2. **Controlling for the carrier by partial correlation.** A gene sits in its carrier's
   chromosome, so the two are collinear by construction. The direct question -- did the
   gene's dose *per carrier genome* change? -- is the trend in **log(gene / clade)**, on the
   samples where both are detected. That is the primary test now; the partial correlation
   is kept in `progression_axis.py` as the secondary.

3. **Depth.** `v2/targets/depth_screen.py` found sequencing depth tracks stage weakly
   (rho -0.12, p 0.046), and detection of a gene depends on depth, so log10 depth joins the
   confounders (Age, BMI, Brinkman, Alcohol, Gender) in every test here.

It also writes **class G** for NB09: for each (GTDB species, target), whether HUMAnN reads of
a curated *producing* enzyme for that target are attributed to the species in >= 10% of
Yachida samples and >= 10 samples (thresholds fixed in CONTRACTS 5d before this ran).
Taxonomy: HUMAnN names (ChocoPhlAn/NCBI) and the cohort's GTDB names are matched after
stripping GTDB's alphabetic suffixes from genus AND epithet ("Ruminococcus_B gnavus" ->
"Ruminococcus gnavus"). A HUMAnN species that corresponds to several GTDB species (GTDB
split it) cannot be attributed to any one of them, so class G is **unavailable** there --
the conservative call -- and the row says why.

The finding is renamed accordingly: **taxonomic turnover within a functional group, in
relative metagenomic space** -- not "source replacement" of a conserved capacity.

Outputs (out/):  compositional_ec_alr.tsv, compositional_gene_clade.tsv,
                 class_g.tsv, compositional_provenance.json
Run:  python v2/genes/compositional.py
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import utils  # noqa: E402
import progression_axis as pa  # noqa: E402

OUT = HERE / "out"
DEPTH = HERE.parent / "targets" / "out" / "sample_depth_proxy.tsv"
G_MIN_FRACTION, G_MIN_SAMPLES = 0.10, 10      # CONTRACTS 5d, class G
STAGE_ENZ = HERE / "stage_enzymes.tsv"                       # CONTRACTS 5d-v
STAGE_CLUSTERS = OUT / "stage_clusters.tsv"                  # from extract_stage_clusters.R
STAGE_ROWS = OUT / "gf_stage_rows.tsv.gz"                    # (Dardel)

# EC -> KEGG compounds that EC PRODUCES, from the pre-specified maps in utils, plus the
# two the HUMAnN EC map misses (4.1.1.105, 4.1.1.25; see README "blind spot").
EC_PRODUCTS: dict[str, set[str]] = {}
for _sub, (_prod, _enz, _ec) in utils.DECARBOXYLATION_SUBSTRATE_PRODUCT.items():
    EC_PRODUCTS.setdefault(_ec, set()).add(_prod)
for _kegg, _ecs in {**utils.POLYAMINE_EC, **utils.ARGININE_ROUTE_EC}.items():
    for _ec in _ecs:
        EC_PRODUCTS.setdefault(_ec, set()).add(_kegg)
EC_PRODUCTS.setdefault("4.1.1.105", set()).update({"C00398", "C05332"})  # Trp/aromatic DC
EC_PRODUCTS.setdefault("4.1.1.25", set()).add("C00483")                  # Tyr DC -> tyramine


def norm_name(name: str) -> str:
    """Binomial with GTDB alphabetic suffixes stripped from genus and epithet."""
    s = name.replace("_", " ").strip()
    parts = s.split()
    if len(parts) < 2:
        return s
    genus = re.sub(r" ?[A-Z]$", "", parts[0].replace(" ", ""))
    genus = re.sub(r"(?<=[a-z])[A-Z]$", "", genus)
    epithet = parts[1]
    return f"{genus} {epithet}"


def norm_gtdb(gtdb: str) -> str:
    """'Ruminococcus_B gnavus' / 'Faecalibacterium prausnitzii_G' -> 'Ruminococcus gnavus'."""
    g, _, e = gtdb.partition(" ")
    g = re.sub(r"_[A-Z]+$", "", g)
    e = re.sub(r"_[A-Z]+$", "", e)
    return f"{g} {e}".strip()


def norm_humann(sp: str) -> str:
    """'g__Blautia.s__Ruminococcus_gnavus' -> 'Ruminococcus gnavus'."""
    name = sp.split(".s__")[-1]
    parts = name.split("_")
    return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else name


def cluster_ecs() -> dict[str, set[str]]:
    t = pd.read_csv(pa.TARGETS, sep="\t", dtype=str)
    c = pd.read_csv(pa.CURATED, sep="\t", dtype=str)
    out: dict[str, set[str]] = {}
    for u, ec in zip(t["uniref90"], t["ec"]):
        out.setdefault(u, set()).add(ec.strip())
    for u, ecs in zip(c["uniref90"], c["ecs"].fillna("")):
        for ec in ecs.split(";"):
            if ec.strip():
                out.setdefault(u, set()).add(ec.strip())
    return out


def design(cw: pd.DataFrame, accessions: list[str]):
    """Stage score, confounders + log10 depth, on the primary-axis samples."""
    g = cw.set_index("Accession")
    score = g["Stage.Group"].map(utils.STAGE_GROUP_NUMERIC)
    Z = utils.build_confounder_matrix(cw)
    Z.index = cw["Accession"].values
    depth = pd.read_csv(DEPTH, sep="\t", dtype={"Sample": str}).set_index("Sample")["log10_bytes"]
    Z["log10_depth"] = g["Sample"].astype(str).map(depth)
    keep = [a for a in accessions if pd.notna(score.get(a))
            and Z.loc[a, pa.CONFOUNDERS + ["log10_depth"]].notna().all()]
    return keep, score.loc[keep].astype(float), Z.loc[keep, pa.CONFOUNDERS + ["log10_depth"]].astype(float)


def trend(x: pd.Series, y: pd.Series, Z: pd.DataFrame) -> tuple[float, float, int]:
    ok = x.notna() & np.isfinite(x)
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    r, p = utils.partial_spearman_grid(x[ok].to_numpy()[:, None],
                                       y[ok].to_numpy()[:, None], Z[ok].to_numpy())
    return float(r[0, 0]), float(p[0, 0]), int(ok.sum())


def main() -> None:
    cw = pa.load_crosswalk()
    acc = cw["Accession"].tolist()
    comm, spw, sp_long, _ = pa.load_long(pa.LONG, acc)
    bg_comm, bg_sp, _, _ = pa.load_long(pa.BACKGROUND_LONG, acc)
    keep, y, Z = design(cw, acc)
    ecs_of = cluster_ecs()

    # ---- 1. EC totals as ALR to the background reference --------------------------------
    bg_total = bg_comm.sum(axis=1)
    ec_cols: dict[str, list[str]] = {}
    for u in comm.columns:
        for ec in ecs_of.get(u, ()):
            ec_cols.setdefault(ec, []).append(u)
    rows = []
    for ec, cols in sorted(ec_cols.items()):
        tot = comm[cols].sum(axis=1)
        pos = tot[tot > 0]
        delta = pos.min() / 2 if len(pos) else np.nan
        alr = np.log((tot.loc[keep] + delta) / bg_total.loc[keep])
        r_alr, p_alr, n = trend(alr, y, Z)
        r_raw, p_raw, _ = trend(tot.loc[keep], y, Z.drop(columns="log10_depth"))
        rows.append({"ec": ec, "n_clusters": len(cols), "n": n,
                     "zero_samples": int((tot.loc[keep] == 0).sum()),
                     "rho_alr_depth_adj": r_alr, "p_alr_depth_adj": p_alr,
                     "rho_share_original_covariates": r_raw, "p_share_original_covariates": p_raw})
    ec_df = pd.DataFrame(rows)
    ec_df["q_alr_depth_adj"] = multipletests(ec_df["p_alr_depth_adj"].fillna(1), method="fdr_bh")[1]
    ec_df.to_csv(OUT / "compositional_ec_alr.tsv", sep="\t", index=False)

    # ---- 2. log(gene / clade) per stratified row -----------------------------------------
    sp_long["uniref90"] = sp_long["row"].str.split("|").str[0]
    sp_long["humann_species"] = sp_long["row"].str.split("|").str[1]
    strat_rows = sorted(sp_long["row"].unique())
    prev = (spw.loc[keep] > 0).sum(axis=0)
    eligible = [r for r in strat_rows
                if prev.get(r, 0) >= max(G_MIN_SAMPLES, G_MIN_FRACTION * len(keep))]
    species_needed = {r.split("|")[1] for r in eligible}
    carriers, carrier_report = pa.carrier_abundance(species_needed, cw)
    # carrier_abundance uses the epithet-only GTDB match; add genus-suffix matches it misses
    header = pd.read_csv(pa.SPECIES_TABLE, sep="\t", nrows=0).columns.tolist()
    gtdb_by_norm: dict[str, list[str]] = {}
    for c in header:
        m = re.search(r"s__(.+)$", c)
        if m:
            gtdb_by_norm.setdefault(norm_gtdb(m.group(1).strip()), []).append(c)
    missing = [s for s in species_needed if s not in carriers and norm_humann(s) in gtdb_by_norm]
    if missing:
        cols = sorted({c for s in missing for c in gtdb_by_norm[norm_humann(s)]})
        df = pd.read_csv(pa.SPECIES_TABLE, sep="\t", usecols=["Sample"] + cols, dtype={"Sample": str})
        df["Accession"] = df["Sample"].map(cw.set_index("Sample")["Accession"])
        df = df.dropna(subset=["Accession"]).set_index("Accession")
        for s in missing:
            carriers[s] = df[gtdb_by_norm[norm_humann(s)]].sum(axis=1)

    grows = []
    for r in eligible:
        u, s = r.split("|", 1)
        rec = {"row": r, "uniref90": u, "humann_species": s,
               "ecs": ";".join(sorted(ecs_of.get(u, ()))),
               "n_samples_detected": int(prev.get(r, 0))}
        if s not in carriers:
            rec.update({"carrier_available": False}); grows.append(rec); continue
        gene = spw.loc[keep, r]
        clade = carriers[s].reindex(keep)
        both = (gene > 0) & (clade > 0)
        ratio = pd.Series(np.nan, index=keep)
        ratio[both] = np.log(gene[both] / clade[both])
        r_ratio, p_ratio, n = trend(ratio, y, Z)
        r_gene, p_gene, _ = trend(gene, y, Z)
        r_clade, p_clade, _ = trend(clade, y, Z)
        rec.update({"carrier_available": True, "n_both_detected": n,
                    "rho_gene_per_clade": r_ratio, "p_gene_per_clade": p_ratio,
                    "rho_gene_share": r_gene, "p_gene_share": p_gene,
                    "rho_clade": r_clade, "p_clade": p_clade})
        grows.append(rec)
    gdf = pd.DataFrame(grows)
    ok = gdf.get("p_gene_per_clade", pd.Series(dtype=float)).notna()
    gdf["q_gene_per_clade"] = np.nan
    if ok.any():
        gdf.loc[ok, "q_gene_per_clade"] = multipletests(gdf.loc[ok, "p_gene_per_clade"], method="fdr_bh")[1]
    gdf.to_csv(OUT / "compositional_gene_clade.tsv", sep="\t", index=False)

    # ---- 3. class G --------------------------------------------------------------------
    all_yachida = cw["Accession"].tolist()
    prev_all = (spw.reindex(all_yachida, fill_value=0) > 0).sum(axis=0)
    n_all = len(all_yachida)
    humann_taxa = set(sp_long["humann_species"]) | set(
        bg_sp.columns.str.split("|").str[1]) if len(bg_sp.columns) else set(sp_long["humann_species"])
    # Stage-associated targets (CONTRACTS 5d-v): producing ECs from stage_enzymes.tsv, the
    # EC route only; per-row cohort prevalence counted on Dardel (extract_stage_clusters.R).
    stage_enz = pd.read_csv(STAGE_ENZ, sep="\t", dtype=str)
    stage_targets: set[str] = set()
    stage_support_rows = []
    if STAGE_ROWS.exists() and STAGE_CLUSTERS.exists():
        ec_to_kegg = stage_enz.groupby("ec")["kegg"].apply(set).to_dict()
        scl = pd.read_csv(STAGE_CLUSTERS, sep="\t", dtype=str)
        stage_ecs_of = scl.groupby("uniref90")["ec"].apply(set).to_dict()
        srows = pd.read_csv(STAGE_ROWS, sep="\t", dtype={"n_nz_cohort": int})
        humann_taxa |= set(srows["species"])
        stage_targets = set(stage_enz["kegg"])
        for u, s, n_nz in zip(srows["uniref90"], srows["species"], srows["n_nz_cohort"]):
            for ec in stage_ecs_of.get(u, ()):
                for k in ec_to_kegg.get(ec, ()):
                    stage_support_rows.append((s, k, u, ec, int(n_nz)))
    else:
        print("WARNING: no stage-target extraction in out/; class G covers the a-priori "
              "targets only (CONTRACTS 5d-v)")
    humann_by_norm: dict[str, set[str]] = {}
    for t in humann_taxa:
        if t and t != "unclassified":
            humann_by_norm.setdefault(norm_humann(t), set()).add(t)
    gtdb_species = sorted({norm_gtdb(m.group(1).strip()) and m.group(1).strip()
                           for c in header for m in [re.search(r"s__(.+)$", c)] if m})
    gtdb_count_by_norm: dict[str, int] = {}
    for g in gtdb_species:
        gtdb_count_by_norm[norm_gtdb(g)] = gtdb_count_by_norm.get(norm_gtdb(g), 0) + 1
    targets_with_enzyme = sorted({k for ks in EC_PRODUCTS.values() for k in ks} | stage_targets)

    support: dict[tuple[str, str], list] = {}
    for r in strat_rows:
        u, s = r.split("|", 1)
        n_nz = int(prev_all.get(r, 0))
        if n_nz < G_MIN_SAMPLES or n_nz < G_MIN_FRACTION * n_all:
            continue
        for ec in ecs_of.get(u, ()):
            for k in EC_PRODUCTS.get(ec, ()):
                support.setdefault((norm_humann(s), k), []).append((u, ec, n_nz))
    for s, k, u, ec, n_nz in stage_support_rows:          # same thresholds, 5d-v item 3
        if n_nz < G_MIN_SAMPLES or n_nz < G_MIN_FRACTION * n_all:
            continue
        support.setdefault((norm_humann(s), k), []).append((u, ec, n_nz))

    grow = []
    for g in gtdb_species:
        key = norm_gtdb(g)
        hm = humann_by_norm.get(key, set())
        n_gtdb = gtdb_count_by_norm.get(key, 0)
        if not hm:
            match = "no_humann_counterpart"
        elif n_gtdb > 1:
            match = "ambiguous_gtdb_split"
        elif len(hm) > 1:
            match = "ambiguous_humann"
        else:
            match = "unique"
        for k in targets_with_enzyme:
            hits = support.get((key, k), []) if match == "unique" else []
            grow.append({"Species": g, "kegg": k, "match_type": match,
                         "Source_gene_level_available": match == "unique",
                         "Source_gene_level_support": bool(hits),
                         "clusters": ";".join(sorted({h[0] for h in hits})),
                         "max_samples_detected": max((h[2] for h in hits), default=0)})
    gcls = pd.DataFrame(grow)
    # Base-rate gate (CONTRACTS 5d-v item 2, council D7): the share of uniquely matched
    # species meeting the read threshold. Above utils.EVIDENCE_MAX_BASE_RATE, the call
    # cannot single out a producer, so G is available but does not support; raw kept.
    uniq = gcls[gcls["match_type"] == "unique"]
    base_rate = uniq.groupby("kegg")["Source_gene_level_support"].mean()
    gcls["gene_level_raw_support"] = gcls["Source_gene_level_support"]
    gcls["gene_level_base_rate"] = gcls["kegg"].map(base_rate)
    gcls["Source_gene_level_support"] = (gcls["gene_level_raw_support"]
                                         & (gcls["gene_level_base_rate"]
                                            <= utils.EVIDENCE_MAX_BASE_RATE))
    gcls["gene_level_route"] = np.where(gcls["kegg"].isin(stage_targets),
                                        "stage_agora2_ec", "curated_class")
    match_counts = gcls.drop_duplicates("Species")["match_type"].value_counts().to_dict()
    # A species with no HUMAnN counterpart is unavailable for class G by definition, and
    # NB09's join defaults to unavailable, so those ~57k rows are not written.
    gcls = gcls[gcls["match_type"] != "no_humann_counterpart"]
    gcls.to_csv(OUT / "class_g.tsv", sep="\t", index=False)

    prov = {"script": "v2/genes/compositional.py",
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "n_primary_axis_samples": len(keep), "covariates": list(Z.columns),
            "background_reference_families": int(bg_comm.shape[1]),
            "ec_tests": int(len(ec_df)),
            "gene_clade_rows_eligible": len(eligible),
            "gene_clade_rows_tested": int(ok.sum()),
            "gene_clade_q05": int((gdf["q_gene_per_clade"] < 0.05).sum()),
            "class_g_species": len(gtdb_species),
            "class_g_match_types": match_counts,
            "class_g_supported_pairs": int(gcls["Source_gene_level_support"].sum()),
            "class_g_raw_supported_pairs": int(gcls["gene_level_raw_support"].sum()),
            "class_g_base_rate_max": utils.EVIDENCE_MAX_BASE_RATE,
            "class_g_base_rates": {k: round(float(v), 4) for k, v in base_rate.items()},
            "class_g_stage_targets_covered": sorted(stage_targets),
            "carrier_report": {k: v for k, v in carrier_report.items() if k != "gtdb_columns_used"}}
    (OUT / "compositional_provenance.json").write_text(json.dumps(prov, indent=1, default=str))
    print(json.dumps(prov, indent=1, default=str))
    print(ec_df.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
