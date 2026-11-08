# Subset curatedMetagenomicData's stratified gene families to the enzyme activities the
# chapter makes claims about, and write everything needed afterwards so the 3.26 GB
# object is never loaded a second time.
#
# Run via extract_gene_families.sbatch. Reads:
#   profiles/gene_families.rda          the matrix, bzip2, 3.26 GB on disk
#   genes/target_clusters.tsv           EC route:      clusters from HUMAnN's EC map
#   genes/curated_clusters.tsv          curated route: clusters holding Swiss-Prot seqs
#
# ---------------------------------------------------------------------------
# Why there are two routes and not one
#
# validate_mapping.py measured HUMAnN's UniRef90 -> EC map against the 363 curated
# Swiss-Prot sequences. Overall sensitivity is 77%, but it is 0% for exactly the three
# activities the chapter leads with:
#
#   EC 4.1.1.105  tryptophan decarboxylase      0 of 2 curated sequences recoverable
#   EC 4.1.1.28   aromatic-L-amino-acid dc      0 of 1
#   EC 4.1.1.25   tyrosine decarboxylase        0 of 4
#
# UniRef90_A7B1V0 (RUMGNA_01526, R. gnavus) and UniRef90_J7SZ64 (CLOSPO_02083,
# C. sporogenes) are both reviewed clusters named "Tryptophan decarboxylase", 20 and 21
# members, and NEITHER appears anywhere in the EC map -- nor in the KO map. Regrouping
# gene families to EC 4.1.1.105 and reporting the result would therefore have measured
# four clusters that do not include either characterised gut bacterial enzyme, and a
# near-zero result would have read as "no bacterial tryptamine synthesis" when it
# actually means "the mapping does not carry these enzymes". That is a false negative of
# the silent kind, so the cluster ids are queried directly as well.
#
# The two routes are reported separately and never summed. EC-route numbers are
# comparable with published HUMAnN EC output; curated-route numbers are the ones that
# speak to characterised enzymes.
# ---------------------------------------------------------------------------
suppressWarnings(suppressMessages({
  ok <- requireNamespace("Matrix", quietly = TRUE)
}))

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else
  "/cfs/klemming/projects/supr/naiss2024-6-169/andnair"
PROF <- file.path(ROOT, "profiles")
GENE <- file.path(ROOT, "genes")
OUT  <- file.path(GENE, "out")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n", sep = "")
t0 <- Sys.time()
lap <- function(what) say("  [%6.1f s] %s", as.numeric(difftime(Sys.time(), t0,
                                                                units = "secs")), what)

# ---------------------------------------------------------------------------
# 1. The target clusters, from both routes
# ---------------------------------------------------------------------------
say("=== target clusters ===")
ec_tab  <- read.delim(file.path(GENE, "target_clusters.tsv"), stringsAsFactors = FALSE)
cur_tab <- read.delim(file.path(GENE, "curated_clusters.tsv"), stringsAsFactors = FALSE)
ec_cl   <- unique(ec_tab$uniref90)
cur_cl  <- unique(cur_tab$uniref90)
all_cl  <- union(ec_cl, cur_cl)
say("  EC route      : %d clusters over %d (ec, cluster) pairs", length(ec_cl), nrow(ec_tab))
say("  curated route : %d clusters holding characterised sequences", length(cur_cl))
say("  curated-only  : %d clusters the EC route would miss", length(setdiff(cur_cl, ec_cl)))
say("  union         : %d clusters to look for", length(all_cl))

# ---------------------------------------------------------------------------
# 2. Load, and assert the object is the study we think it is
#
# The replan documents that ExperimentHub's fetch endpoint keys on rdatapaths.id and
# that a wrong id returns a valid .rda from a DIFFERENT study rather than an error
# (fetching 7254 yields MatsonV_2018.relative_abundance). So the only safe check is
# that the object names itself, and every load in this programme asserts it.
# ---------------------------------------------------------------------------
say("=== loading gene_families.rda ===")
f <- file.path(PROF, "gene_families.rda")
say("  on disk       : %s B", format(file.size(f), big.mark = ","))
env <- new.env()
nm <- load(f, envir = env)
lap("load() returned")
say("  object name   : %s", paste(nm, collapse = ", "))
if (!any(grepl("YachidaS_2019", nm))) {
  stop("object does not name itself YachidaS_2019 -- this is the WRONG STUDY. ",
       "Got: ", paste(nm, collapse = ", "))
}
say("  names itself YachidaS_2019: YES")

gf <- get(nm[1], envir = env)
say("  class         : %s", paste(class(gf), collapse = "/"))
say("  typeof        : %s", typeof(gf))
say("  dim           : %s", paste(dim(gf), collapse = " x "))
say("  object.size   : %s", format(object.size(gf), units = "GB"))
if (inherits(gf, "sparseMatrix")) {
  say("  stored sparse : %d nonzero of %.3g cells (%.2f%%)",
      length(gf@x), prod(dim(gf)), 100 * length(gf@x) / prod(dim(gf)))
}

rn <- rownames(gf)
cn <- colnames(gf)
if (is.null(rn)) stop("no rownames: cannot identify gene families")
say("  rows          : %d", length(rn))
say("  cols          : %d", length(cn))

strat <- grepl("|g__", rn, fixed = TRUE)
say("  stratified    : %d of %d (%.1f%%)", sum(strat), length(rn),
    100 * sum(strat) / length(rn))

# ---------------------------------------------------------------------------
# 3. The row-name index, written once so nothing reloads this object
# ---------------------------------------------------------------------------
say("=== writing the row-name index ===")
gzf <- gzfile(file.path(OUT, "rownames_index.txt.gz"), "w")
writeLines(rn, gzf); close(gzf)
writeLines(cn, file.path(OUT, "colnames.txt"))
lap(sprintf("index written (%s)",
            format(file.size(file.path(OUT, "rownames_index.txt.gz")), big.mark = ",")))

# ---------------------------------------------------------------------------
# 4. Subset. The cluster id is the part before the "|"; the species suffix, where
#    present, is what makes this the attribution step rather than a gene count.
# ---------------------------------------------------------------------------
say("=== matching target clusters ===")
base_id <- sub("\\|.*$", "", rn)
# %in% hashes internally; a per-row closure over millions of rows does not.
hit <- base_id %in% all_cl
lap(sprintf("matched: %d rows", sum(hit)))
say("  distinct clusters found : %d of %d targets",
    length(unique(base_id[hit])), length(all_cl))
say("  of the matched rows, %d are species-stratified", sum(hit & strat))

if (!any(hit)) stop("no target cluster appears in this matrix at all -- check the ",
                    "cluster id spelling against rownames_index.txt.gz")

sub_m <- gf[hit, , drop = FALSE]
rm(gf, env); invisible(gc())
lap("subset taken, full object released")

saveRDS(sub_m, file.path(OUT, "gf_targets.rds"), compress = "gzip")
lap(sprintf("gf_targets.rds (%s)",
            format(file.size(file.path(OUT, "gf_targets.rds")), big.mark = ",")))

# ---------------------------------------------------------------------------
# 5. A per-row summary and a long-format table of nonzeros. Long format is what the
#    downstream join to the metabolite data wants, and keeping only nonzeros makes it
#    small enough to pull to a laptop.
# ---------------------------------------------------------------------------
say("=== summarising ===")
srn  <- rownames(sub_m)
sid  <- sub("\\|.*$", "", srn)
ssp  <- ifelse(grepl("|g__", srn, fixed = TRUE), sub("^.*\\|", "", srn), NA_character_)

dm <- as.matrix(sub_m)
nz <- dm != 0
summ <- data.frame(
  row          = srn,
  uniref90     = sid,
  species      = ssp,
  stratified   = !is.na(ssp),
  n_samples_nz = as.integer(rowSums(nz)),
  prevalence   = rowSums(nz) / ncol(dm),
  mean_nz      = ifelse(rowSums(nz) > 0, rowSums(dm) / pmax(rowSums(nz), 1), 0),
  max          = apply(dm, 1, max),
  stringsAsFactors = FALSE
)
# annotate which route(s) asked for each cluster, so no reader has to guess
summ$in_ec_route      <- summ$uniref90 %in% ec_cl
summ$in_curated_route <- summ$uniref90 %in% cur_cl
ec_lab <- tapply(ec_tab$ec, ec_tab$uniref90, function(x) paste(sort(unique(x)), collapse = ";"))
summ$ecs <- unname(ec_lab[summ$uniref90])
write.table(summ, file.path(OUT, "gf_targets_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "")
lap("gf_targets_summary.tsv")

idx <- which(nz, arr.ind = TRUE)
long <- data.frame(row = srn[idx[, "row"]], sample = cn[idx[, "col"]],
                   abundance = dm[idx], stringsAsFactors = FALSE)
gzl <- gzfile(file.path(OUT, "gf_targets_long.tsv.gz"), "w")
write.table(long, gzl, sep = "\t", quote = FALSE, row.names = FALSE)
close(gzl)
lap(sprintf("gf_targets_long.tsv.gz (%d nonzero entries, %s)", nrow(long),
            format(file.size(file.path(OUT, "gf_targets_long.tsv.gz")), big.mark = ",")))

# ---------------------------------------------------------------------------
# 6. The two questions the chapter cannot currently answer, asked directly
# ---------------------------------------------------------------------------
say("=== the named tryptophan decarboxylases, by cluster id ===")
for (c in c("UniRef90_A7B1V0", "UniRef90_J7SZ64")) {
  rows <- srn[sid == c]
  if (!length(rows)) {
    say("  %s : not present in this cohort's gene families", c)
    next
  }
  say("  %s : %d row(s)", c, length(rows))
  for (r in rows) {
    v <- dm[r, ]
    say("      %-62s nonzero in %3d/%d samples, max %.3g",
        r, sum(v != 0), length(v), max(v))
  }
}

say("=== lysine decarboxylase under Allisonella (EC 4.1.1.18) ===")
ldc <- unique(ec_tab$uniref90[ec_tab$ec == "4.1.1.18"])
ldc <- union(ldc, cur_tab$uniref90[grepl("4.1.1.18", cur_tab$ecs)])
alli <- srn[sid %in% ldc & grepl("g__Allisonella", srn, fixed = TRUE)]
if (length(alli)) {
  say("  %d Allisonella-stratified lysine decarboxylase row(s):", length(alli))
  for (r in alli) {
    v <- dm[r, ]
    say("      %-62s nonzero in %3d/%d samples, max %.3g",
        r, sum(v != 0), length(v), max(v))
  }
} else {
  say("  none. Either Allisonella carries no lysine decarboxylase in these samples,")
  say("  or no Allisonella gene family is stratified at all -- which of the two it is")
  say("  has to be read off rownames_index.txt.gz, not assumed.")
  n_alli <- sum(grepl("g__Allisonella", rn, fixed = TRUE))
  say("  Allisonella-stratified rows anywhere in the matrix: %d", n_alli)
}

# ---------------------------------------------------------------------------
# 7. Provenance
# ---------------------------------------------------------------------------
prov <- sprintf('{
  "built_utc": "%s",
  "source_rda": "%s",
  "source_bytes": %d,
  "object_name": "%s",
  "self_named_assertion": "passed: object name contains YachidaS_2019",
  "class": "%s",
  "dim": [%d, %d],
  "rows_stratified": %d,
  "ec_route_clusters": %d,
  "curated_route_clusters": %d,
  "curated_only_clusters": %d,
  "union_clusters": %d,
  "rows_matched": %d,
  "rows_matched_stratified": %d,
  "clusters_found": %d,
  "nonzero_entries": %d,
  "R_version": "%s",
  "elapsed_s": %.1f
}
', format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), f, file.size(f),
   nm[1], paste(class(sub_m), collapse = "/"), nrow(sub_m), ncol(sub_m), sum(strat),
   length(ec_cl), length(cur_cl), length(setdiff(cur_cl, ec_cl)), length(all_cl),
   sum(hit), sum(hit & strat), length(unique(sid)), nrow(long),
   R.version.string, as.numeric(difftime(Sys.time(), t0, units = "secs")))
writeLines(prov, file.path(OUT, "extract_provenance.json"))

say("=== done in %.1f min ===", as.numeric(difftime(Sys.time(), t0, units = "mins")))
print(list.files(OUT, full.names = FALSE))
