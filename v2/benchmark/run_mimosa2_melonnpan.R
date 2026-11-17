#!/usr/bin/env Rscript
# Benchmark step 2/3 — run MelonnPan and MIMOSA2 on the exported Yachida matrices.
#
# >>> VERIFY before first run (environment-specific, cannot be tested from the pipeline):
#   1. Packages:  MelonnPan (biobakery/melonnpan) and MIMOSA2 (this repo, ../MIMOSA2).
#   2. MIMOSA2 reference data: AGORA/KEGG files. Set MIMOSA2_KEGG_PREFIX / MIMOSA2_DATA_PREFIX
#      below (MIMOSA2::download_reference_data() fetches them if absent).
#   3. The MIMOSA2 return-object field used for contributions (var_shares) — the code
#      writes whatever the installed version returns and warns if the field name differs.
#
# Outputs (benchmark/outputs/):
#   melonnpan_weights.tsv        taxa x metabolites EN weights  (predictor importance)
#   melonnpan_trainperf.tsv      per-metabolite training Spearman (well-predicted set)
#   mimosa2_contributions.tsv    compound, taxon, var_share     (mechanistic contribution)

suppressWarnings(suppressMessages({
  library(data.table)
}))

here     <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
IN       <- file.path(here, "inputs")
OUT      <- file.path(here, "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# here = <project root>/v2/benchmark, so two levels up is the project root.
MIMOSA2_REPO        <- normalizePath(file.path(here, "..", "..", "MIMOSA2"), mustWork = FALSE)  # <project root>/MIMOSA2
MIMOSA2_KEGG_PREFIX <- Sys.getenv("MIMOSA2_KEGG_PREFIX", unset = file.path(here, "ref", "KEGG"))
MIMOSA2_DATA_PREFIX <- Sys.getenv("MIMOSA2_DATA_PREFIX", unset = file.path(here, "ref"))

# ======================================================================= MelonnPan
message("== MelonnPan ==")
mp_ok <- requireNamespace("melonnpan", quietly = TRUE)
if (!mp_ok) {
  message("  melonnpan not installed — skipping. Install: ",
          "devtools::install_github('biobakery/melonnpan')")
} else {
  metag <- as.data.frame(fread(file.path(IN, "mp_metag_samples_x_species.tsv")))
  metab <- as.data.frame(fread(file.path(IN, "mp_metab_samples_x_metabolites.tsv")))
  rownames(metag) <- metag$Sample; metag$Sample <- NULL
  rownames(metab) <- metab$Sample; metab$Sample <- NULL
  common <- intersect(rownames(metag), rownames(metab))
  metag <- metag[common, , drop = FALSE]; metab <- metab[common, , drop = FALSE]
  message(sprintf("  samples=%d, taxa=%d, mets=%d", length(common), ncol(metag), ncol(metab)))

  # melonnpan.train fits an elastic net per metabolite and writes weights + training rho.
  mp_out <- file.path(OUT, "melonnpan_run")
  dir.create(mp_out, showWarnings = FALSE)
  fit <- tryCatch(
    melonnpan::melonnpan.train(metab = metab, metag = metag, output = mp_out),
    error = function(e) { message("  melonnpan.train error: ", conditionMessage(e)); NULL })

  # Weights file written by melonnpan.train (taxa x metabolites). Name varies by version.
  wfile <- list.files(mp_out, pattern = "Weight", full.names = TRUE, ignore.case = TRUE)
  if (length(wfile)) {
    file.copy(wfile[1], file.path(OUT, "melonnpan_weights.tsv"), overwrite = TRUE)
    message("  weights -> melonnpan_weights.tsv")
  } else if (!is.null(fit) && !is.null(fit$weights)) {
    fwrite(data.table(fit$weights, keep.rownames = "feature"),
           file.path(OUT, "melonnpan_weights.tsv"), sep = "\t")
  } else message("  >>> VERIFY: no MelonnPan weights file found in ", mp_out)

  pfile <- list.files(mp_out, pattern = "(perf|spearman|metric)", full.names = TRUE, ignore.case = TRUE)
  if (length(pfile)) file.copy(pfile[1], file.path(OUT, "melonnpan_trainperf.tsv"), overwrite = TRUE)
}

# ========================================================================= MIMOSA2
message("== MIMOSA2 ==")
m2_ok <- requireNamespace("mimosa", quietly = TRUE) || requireNamespace("MIMOSA2", quietly = TRUE)
if (!m2_ok && dir.exists(MIMOSA2_REPO)) {
  message("  loading MIMOSA2 from source repo: ", MIMOSA2_REPO)
  m2_ok <- tryCatch({ devtools::load_all(MIMOSA2_REPO, quiet = TRUE); TRUE },
                    error = function(e) { message("  load_all failed: ", conditionMessage(e)); FALSE })
}
if (!m2_ok) {
  message("  MIMOSA2 unavailable — skipping. See repo at ", MIMOSA2_REPO)
} else {
  # Config table per MIMOSA2's format (cf. tests/testthat/configs_AGORA.txt).
  # Species abundances mapped to AGORA genome-scale models; metabolites as KEGG IDs.
  cfg <- data.table(
    V1 = c("file1", "file2", "ref_choices", "file1_type",
           "metType", "kegg_prefix", "data_prefix", "simThreshold", "signifThreshold"),
    V2 = c(file.path(IN, "mimosa2_species_x_samples.tsv"),
           file.path(IN, "mimosa2_mets_x_samples.tsv"),
           "AGORA genomes and models",       # >>> VERIFY: exact ref label in your MIMOSA2
           "Greengenes 13_5 or 13_8 OTUs",   # >>> VERIFY: taxonomy input mode for species names
           "KEGG Compound IDs",
           MIMOSA2_KEGG_PREFIX, MIMOSA2_DATA_PREFIX, "0.95", "0.1"))
  fwrite(cfg, file.path(OUT, "mimosa2_config_used.tsv"), sep = "\t", col.names = FALSE)

  res <- tryCatch(
    MIMOSA2::run_mimosa2(cfg),
    error = function(e) tryCatch(mimosa::run_mimosa2(cfg),
            error = function(e2) { message("  run_mimosa2 error: ", conditionMessage(e2)); NULL }))

  # Extract per-(taxon, compound) contribution (variance share). Field name varies;
  # try the common ones and warn if none match.
  contrib <- NULL
  if (!is.null(res)) {
    for (nm in c("var_shares", "varShares", "contributions", "species_contributions")) {
      if (!is.null(res[[nm]])) { contrib <- as.data.table(res[[nm]]); break }
    }
  }
  if (is.null(contrib)) {
    message("  >>> VERIFY: could not find a contributions field in run_mimosa2() output; ",
            "names = ", paste(names(res), collapse = ", "))
  } else {
    fwrite(contrib, file.path(OUT, "mimosa2_contributions.tsv"), sep = "\t")
    message("  contributions -> mimosa2_contributions.tsv (", nrow(contrib), " rows)")
  }
}

message("Done. Next: python compare_frameworks.py")
