# Class G for the stage-associated targets (CONTRACTS 5d-v): how many cohort samples carry
# reads of each stage-target producing enzyme, attributed to each species.
#
# Run via extract_stage_clusters.sbatch. Reads (all under ROOT):
#   genes/stage_enzymes.tsv               EC list from stage_enzymes.py (committed first)
#   genes/stage_accessions.txt            the 347 Yachida cohort accessions, one per line
#   mappings/map_level4ec_uniref90.txt.gz HUMAnN's EC -> UniRef90 map (Zenodo 17100034)
#   profiles/gene_families.rda            the matrix, 3.26 GB, asserted to be YachidaS_2019
# Writes genes/out/:
#   stage_clusters.tsv        ec, uniref90 (the EC route; there is no curated arm here)
#   gf_stage_rows.tsv.gz      per species-stratified row: n_nz_cohort over the 347 accessions
#   stage_extract_provenance.json
#
# Only per-row prevalence leaves Dardel: that is all class G uses (the >= 10% and >= 10
# sample rule), and these ECs include ubiquitous enzymes whose full long table would be
# large for no purpose.
suppressWarnings(suppressMessages(library(Matrix)))

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else
  "/cfs/klemming/projects/supr/naiss2024-6-169/andnair"
GENE <- file.path(ROOT, "genes")
OUT  <- file.path(GENE, "out")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n", sep = "")
t0 <- Sys.time()
lap <- function(what) say("  [%6.1f s] %s", as.numeric(difftime(Sys.time(), t0, units = "secs")), what)

enz <- read.delim(file.path(GENE, "stage_enzymes.tsv"), stringsAsFactors = FALSE,
                  colClasses = "character")
ecs <- unique(enz$ec)
acc <- readLines(file.path(GENE, "stage_accessions.txt"))
acc <- acc[nzchar(acc)]
say("stage ECs: %d   cohort accessions: %d", length(ecs), length(acc))
if (length(acc) != 347) stop("expected 347 cohort accessions, got ", length(acc))

# --- EC -> clusters -------------------------------------------------------------------
say("=== EC map ===")
con <- gzfile(file.path(ROOT, "mappings", "map_level4ec_uniref90.txt.gz"), "r")
pairs <- list()
repeat {
  ln <- readLines(con, n = 5000)
  if (!length(ln)) break
  ec <- sub("\t.*$", "", ln)
  hit <- which(ec %in% ecs)
  for (i in hit) {
    f <- strsplit(ln[i], "\t", fixed = TRUE)[[1]]
    cl <- f[-1][nzchar(f[-1])]
    pairs[[length(pairs) + 1]] <- data.frame(ec = f[1], uniref90 = cl, stringsAsFactors = FALSE)
  }
}
close(con)
ecmap <- if (length(pairs)) do.call(rbind, pairs) else
  data.frame(ec = character(), uniref90 = character())
write.table(ecmap, file.path(OUT, "stage_clusters.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
missing <- setdiff(ecs, ecmap$ec)
say("  ECs in map: %d of %d; clusters: %d", length(unique(ecmap$ec)), length(ecs),
    length(unique(ecmap$uniref90)))
if (length(missing)) say("  NOT IN THE MAP: %s", paste(missing, collapse = ", "))
lap("EC map read")

# --- load, assert, subset ---------------------------------------------------------------
say("=== loading gene_families.rda ===")
env <- new.env()
nm <- load(file.path(ROOT, "profiles", "gene_families.rda"), envir = env)
if (!any(grepl("YachidaS_2019", nm))) stop("object does not name itself YachidaS_2019: ", nm)
gf <- get(nm[1], envir = env)
lap(sprintf("loaded %s, %s", nm[1], paste(dim(gf), collapse = " x ")))
rn <- rownames(gf)
cn <- colnames(gf)
miss_acc <- setdiff(acc, cn)
if (length(miss_acc)) stop(length(miss_acc), " cohort accessions are not columns")
strat <- grepl("|g__", rn, fixed = TRUE)
keep <- which(strat & sub("\\|.*$", "", rn) %in% unique(ecmap$uniref90))
say("  species-stratified rows in the stage clusters: %d", length(keep))
sub_m <- gf[keep, acc, drop = FALSE]
rm(gf, env); invisible(gc())
lap("subset taken, full object released")

nz <- as.integer(Matrix::rowSums(sub_m != 0))
srn <- rownames(sub_m)
res <- data.frame(row = srn, uniref90 = sub("\\|.*$", "", srn),
                  species = sub("^.*\\|", "", srn), n_nz_cohort = nz,
                  stringsAsFactors = FALSE)
res <- res[res$n_nz_cohort > 0, ]
gz <- gzfile(file.path(OUT, "gf_stage_rows.tsv.gz"), "w")
write.table(res, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
lap(sprintf("gf_stage_rows.tsv.gz: %d rows with any read", nrow(res)))

prov <- sprintf(paste0('{\n "script": "v2/genes/extract_stage_clusters.R",\n',
                       ' "run_utc": "%s",\n "object": "%s",\n "ecs_requested": %d,\n',
                       ' "ecs_in_map": %d,\n "ecs_missing": [%s],\n "clusters": %d,\n',
                       ' "stratified_rows_matched": %d,\n "rows_with_reads": %d,\n',
                       ' "cohort_accessions": %d,\n "seconds": %.0f\n}\n'),
                format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), nm[1], length(ecs),
                length(unique(ecmap$ec)),
                if (length(missing)) paste0('"', missing, '"', collapse = ", ") else "",
                length(unique(ecmap$uniref90)), length(keep), nrow(res), length(acc),
                as.numeric(difftime(Sys.time(), t0, units = "secs")))
writeLines(prov, file.path(OUT, "stage_extract_provenance.json"))
say("done")
