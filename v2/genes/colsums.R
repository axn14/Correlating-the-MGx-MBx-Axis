#!/usr/bin/env Rscript
# Per-sample normalization check for gene_families.rda.
#
# Submitted as (job 27324295, COMPLETED in 5m02s, MaxRSS 5.2 GiB):
#
#   sbatch -A <project> -J ki-colsums -p shared -n 1 -c 1 --mem=12G -t 00:30:00 \
#          --wrap 'ml cray-R/4.4.0; Rscript --vanilla $ROOT/genes/colsums.R $ROOT'
#
# The memory figure comes from extract_gene_families.sbatch, which measured the load
# at 5.5 GiB; this job adds one logical vector over the nonzeros, 379 MB.
#
# The figures compare one gene family's abundance between groups of samples, so the
# question that has to be settled before any log2FC is drawn is whether the stored
# values are already compositional. If they are relative abundances, the community-level
# rows -- the ones with no `|g__...` suffix -- sum to 1 in every column and no further
# normalization is needed. If they are RPKs, the column sums vary with sequencing depth
# and a between-group difference could be depth alone.
#
# This is deliberately the only thing the job does. load() of the 3.26 GB object is the
# whole cost (288 s); everything after it is two row masks and two column sums.
suppressWarnings(suppressMessages({
  library(Matrix)
}))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "/cfs/klemming/projects/supr/naiss2024-6-169/andnair"
f <- file.path(root, "profiles", "gene_families.rda")
stopifnot(file.exists(f))

env <- new.env(parent = emptyenv())
t0 <- Sys.time()
nm <- load(f, envir = env)
if (!any(grepl("YachidaS_2019", nm))) {
  stop("object does not name itself YachidaS_2019 -- this is the WRONG STUDY. ",
       "Got: ", paste(nm, collapse = ", "))
}
gf <- get(nm[1], envir = env)
cat(sprintf("loaded %s  %s  %d x %d  in %.0f s\n", nm[1], class(gf)[1],
            nrow(gf), ncol(gf), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

rn <- rownames(gf)
strat <- grepl("|", rn, fixed = TRUE)
comm <- !strat
cat(sprintf("community rows %d   stratified rows %d\n", sum(comm), sum(strat)))

# Neither `colSums(gf[comm, ])` nor a TsparseMatrix conversion fits in memory here.
# Subsetting 4.4 M rows allocates a second matrix; `as(gf, "TsparseMatrix")` materializes
# an explicit column index for all 379,029,747 nonzeros, which is 1.5 GB of integers on
# top of a 6 GB copy of i and x -- about 10.5 GB peak against a 12G request.
#
# A dgCMatrix already stores its nonzeros grouped by column: column j occupies slots
# (p[j]+1):p[j+1] of i and x. So walk the 616 columns and sum each slice in place. The
# only large temporary is the row mask lifted onto the nonzeros, 379 M logicals = 379 MB,
# and peak stays near the 5.5 GiB the load itself needs.
p  <- gf@p
ii <- gf@i + 1L
xx <- gf@x
is_comm <- comm[ii]
rm(ii); invisible(gc())

nc <- ncol(gf)
sum_comm  <- numeric(nc)
sum_strat <- numeric(nc)
nz_comm   <- numeric(nc)
for (j in seq_len(nc)) {
  lo <- p[j] + 1L
  hi <- p[j + 1L]
  if (hi < lo) next
  idx <- lo:hi
  cm <- is_comm[idx]
  v  <- xx[idx]
  sum_comm[j]  <- sum(v[cm])
  sum_strat[j] <- sum(v[!cm])
  nz_comm[j]   <- sum(cm)
}

out <- data.frame(sample = colnames(gf), sum_community = sum_comm,
                  sum_stratified = sum_strat, nonzero_community = nz_comm,
                  stringsAsFactors = FALSE)

d <- file.path(root, "genes", "out")
dir.create(d, recursive = TRUE, showWarnings = FALSE)
write.table(out, file.path(d, "sample_normalization.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\n-- community-level column sums --\n")
print(summary(sum_comm))
cat("\n-- stratified column sums --\n")
print(summary(sum_strat))
cat("\n-- nonzero community gene families per sample --\n")
print(summary(nz_comm))
cat(sprintf("\nverdict: community sums in [%.6f, %.6f]; %s\n",
            min(sum_comm), max(sum_comm),
            if (all(abs(sum_comm - 1) < 1e-4)) "RELATIVE ABUNDANCE (sum to 1), no renormalization needed"
            else "NOT compositional -- renormalize by sum_community before comparing samples"))
cat("wrote", file.path(d, "sample_normalization.tsv"), "\n")
