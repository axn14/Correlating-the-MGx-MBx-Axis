#!/usr/bin/env Rscript
# A prevalence-matched background set of gene families, for the question the
# progression-axis trends cannot answer on their own.
#
# Submitted as (job 27324874, COMPLETED in 5m02s, MaxRSS 5.2 GiB):
#
#   sbatch -A <project> -J ki-background -p shared -n 1 -c 1 --mem=16G -t 00:40:00 \
#          --wrap 'ml cray-R/4.4.0; Rscript --vanilla $ROOT/genes/background.R $ROOT'
#
# 16G rather than the 12G the other two jobs use, because subsetting 6,335 rows out of
# a 379-million-nonzero dgCMatrix allocates a second structure while both are live.
#
# THE PROBLEM THIS SOLVES
# 12 of 455 target clusters trend monotonically with CRC stage at BH q < 0.05.
# That is more than BH allows by chance within the target set -- but it says
# nothing about whether trending with stage is RARE. If a tenth of every gene
# family in the gut tracks the stage axis, because the whole community shifts,
# then 12/455 is unremarkable and the enzymes are not special. The only way to
# know is to run the identical test on gene families chosen for having nothing to
# do with these metabolites.
#
# MATCHING, AND WHY ON PREVALENCE
# A rank-correlation's power depends almost entirely on how many samples carry a
# nonzero value, and the target clusters span 0.2% to 98% prevalence. An unmatched
# random draw from 3.58 M community rows would be far sparser than the targets
# (most gene families are near-absent), would therefore have less power, and the
# comparison would flatter the targets for a purely statistical reason. So the
# background is drawn to reproduce the targets' prevalence distribution decile by
# decile, at BACKGROUND_MULT clusters per target, excluding every target cluster
# id in either route.
#
# Output is the same long (row, sample, abundance) shape as gf_targets_long, so
# progression_axis.py runs the same code over it.
suppressWarnings(suppressMessages(library(Matrix)))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "/cfs/klemming/projects/supr/naiss2024-6-169/andnair"
BACKGROUND_MULT <- 5L
SEED <- 20261108L

f <- file.path(root, "profiles", "gene_families.rda")
stopifnot(file.exists(f))
env <- new.env(parent = emptyenv())
t0 <- Sys.time()
nm <- load(f, envir = env)
if (!any(grepl("YachidaS_2019", nm))) {
  stop("object does not name itself YachidaS_2019 -- WRONG STUDY. Got: ",
       paste(nm, collapse = ", "))
}
gf <- get(nm[1], envir = env)
cat(sprintf("loaded %s %d x %d in %.0f s\n", class(gf)[1], nrow(gf), ncol(gf),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

rn <- rownames(gf)
community <- !grepl("|", rn, fixed = TRUE)
cat(sprintf("community rows: %d\n", sum(community)))

# Nonzero count per row, straight off the sparse row index. tabulate() over
# 379,029,747 integers is one pass; building a logical matrix would not fit.
nz_per_row <- tabulate(gf@i + 1L, nbins = nrow(gf))

# Target ids: every cluster either route asked for, excluded from the pool
# whether or not it was observed.
tg <- read.delim(file.path(root, "genes", "target_clusters.tsv"),
                 stringsAsFactors = FALSE)
cu <- read.delim(file.path(root, "genes", "curated_clusters.tsv"),
                 stringsAsFactors = FALSE)
targets <- unique(c(tg$uniref90, cu$uniref90))
cat(sprintf("target cluster ids: %d\n", length(targets)))

is_target <- rn %in% targets
obs_target <- which(community & is_target & nz_per_row > 0)
pool <- which(community & !is_target & nz_per_row > 0)
cat(sprintf("observed target community rows: %d   eligible pool: %d\n",
            length(obs_target), length(pool)))

# Decile bins of the targets' own prevalence, then the same shape from the pool.
tgt_nz <- nz_per_row[obs_target]
brks <- unique(quantile(tgt_nz, probs = seq(0, 1, by = 0.1), type = 1))
brks[1] <- brks[1] - 1L
brks[length(brks)] <- max(nz_per_row) + 1L
tgt_bin <- cut(tgt_nz, breaks = brks, labels = FALSE, right = TRUE)
pool_bin <- cut(nz_per_row[pool], breaks = brks, labels = FALSE, right = TRUE)

set.seed(SEED)
picked <- integer(0)
report <- data.frame()
for (b in sort(unique(tgt_bin))) {
  want <- sum(tgt_bin == b) * BACKGROUND_MULT
  cand <- pool[which(pool_bin == b)]
  take <- if (length(cand) <= want) cand else sample(cand, want)
  picked <- c(picked, take)
  report <- rbind(report, data.frame(bin = b, prev_lo = brks[b], prev_hi = brks[b + 1],
                                     n_target = sum(tgt_bin == b), wanted = want,
                                     available = length(cand), taken = length(take)))
}
print(report)
cat(sprintf("background clusters picked: %d\n", length(picked)))

sub <- gf[picked, , drop = FALSE]
rm(gf, env); invisible(gc())

tri <- as(sub, "TsparseMatrix")
long <- data.frame(row = rownames(sub)[tri@i + 1L],
                   sample = colnames(sub)[tri@j + 1L],
                   abundance = tri@x, stringsAsFactors = FALSE)
long <- long[long$abundance > 0, , drop = FALSE]
cat(sprintf("nonzero triples: %d\n", nrow(long)))

d <- file.path(root, "genes", "out")
dir.create(d, recursive = TRUE, showWarnings = FALSE)
gz <- gzfile(file.path(d, "gf_background_long.tsv.gz"), "w")
write.table(long, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
write.table(report, file.path(d, "background_matching.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(jsonlite_lite <- sprintf(
  '{\n  "built_utc": "%s",\n  "seed": %d,\n  "background_mult": %d,\n  "community_rows": %d,\n  "target_ids": %d,\n  "observed_target_community_rows": %d,\n  "eligible_pool": %d,\n  "background_clusters": %d,\n  "nonzero_triples": %d\n}',
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), SEED, BACKGROUND_MULT,
  sum(community), length(targets), length(obs_target), length(pool),
  length(picked), nrow(long)),
  file.path(d, "background_provenance.json"))
cat("wrote", file.path(d, "gf_background_long.tsv.gz"), "\n")
