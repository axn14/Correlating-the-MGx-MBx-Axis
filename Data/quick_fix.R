################################################################################
# MINIMAL STANDALONE FIX
################################################################################

library(tidyverse)
library(mixOmics)

# Fixed preprocessing function
preprocess_omics_data_fixed <- function(data_obj, min_prevalence = 0.1) {
  
  cat("\n=== Preprocessing Data ===\n")
  
  # Convert to matrix - samples as rows, features as columns
  species_mat <- data_obj$species %>%
    column_to_rownames(var = "Sample") %>%
    as.matrix() %>%
    t()
  
  mtb_mat <- data_obj$metabolites %>%
    column_to_rownames(var = "Sample") %>%
    as.matrix() %>%
    t()
  
  cat(paste("Species:", nrow(species_mat), "samples x", ncol(species_mat), "features\n"))
  cat(paste("Metabolites:", nrow(mtb_mat), "samples x", ncol(mtb_mat), "features\n"))
  
  # Filter by prevalence
  species_prev <- colSums(species_mat > 0) / nrow(species_mat)
  mtb_prev <- colSums(mtb_mat > 0) / nrow(mtb_mat)
  
  species_mat <- species_mat[, species_prev >= min_prevalence]
  mtb_mat <- mtb_mat[, mtb_prev >= min_prevalence]
  
  cat(paste("After filtering (", min_prevalence*100, "% prevalence):\n", sep=""))
  cat(paste("  Species:", ncol(species_mat), "\n"))
  cat(paste("  Metabolites:", ncol(mtb_mat), "\n"))
  
  # Log transform
  species_log <- log10(species_mat + 1e-6)
  mtb_log <- log10(mtb_mat + 1e-6)
  
  # Scale
  species_scaled <- scale(species_log, center = TRUE, scale = TRUE)
  mtb_scaled <- scale(mtb_log, center = TRUE, scale = TRUE)
  
  # Clean up
  species_scaled[is.na(species_scaled) | is.infinite(species_scaled)] <- 0
  mtb_scaled[is.na(mtb_scaled) | is.infinite(mtb_scaled)] <- 0
  
  # Convert to plain matrix
  species_scaled <- as.matrix(species_scaled)
  mtb_scaled <- as.matrix(mtb_scaled)
  
  cat("\nFinal matrices ready!\n")
  cat(paste("  Species range: [", round(min(species_scaled),2), ",", 
            round(max(species_scaled),2), "]\n"))
  cat(paste("  Metabolites range: [", round(min(mtb_scaled),2), ",", 
            round(max(mtb_scaled),2), "]\n"))
  
  return(list(
    species = species_scaled,
    metabolites = mtb_scaled,
    metadata = data_obj$metadata,
    mtb_map = data_obj$mtb_map
  ))
}

# Run preprocessing
cat("\n### RUNNING FIXED PREPROCESSING ###\n")
processed_data_fixed <- preprocess_omics_data_fixed(data, min_prevalence = 0.05)

# Quick mixOmics run
cat("\n### RUNNING MIXOMICS ###\n")

X <- processed_data_fixed$species
Y <- processed_data_fixed$metabolites

cat(paste("X:", nrow(X), "x", ncol(X), "\n"))
cat(paste("Y:", nrow(Y), "x", ncol(Y), "\n"))

# Adaptive parameters
ncomp <- 3
keepX <- min(50, floor(ncol(X) * 0.3))
keepY <- min(30, floor(ncol(Y) * 0.5))

cat(paste("Parameters: ncomp=", ncomp, ", keepX=", keepX, ", keepY=", keepY, "\n", sep=""))

# Run sPLS
result_scca <- spls(
  X = X,
  Y = Y,
  ncomp = ncomp,
  mode = "regression",
  keepX = rep(keepX, ncomp),
  keepY = rep(keepY, ncomp)
)

cat("\n✓ mixOmics completed successfully!\n")

# Make plots
pdf("mixomics_plots.pdf", width = 12, height = 8)

plotIndiv(result_scca, comp = c(1,2), ind.names = FALSE,
          title = "sPLS - Sample Space")

plotVar(result_scca, comp = c(1,2), cex = c(2,2),
        title = "sPLS - Variable Correlation")

plotArrow(result_scca, comp = c(1,2), ind.names = FALSE,
          title = "sPLS - Sample Trajectories")

dev.off()

cat("\nPlots saved to: mixomics_plots.pdf\n")

# Extract top features
loadings_X <- selectVar(result_scca, comp = 1)
loadings_Y <- selectVar(result_scca, comp = 1, block = "Y")

cat("\n### TOP SPECIES (Component 1) ###\n")
print(head(loadings_X$value, 10))

cat("\n### TOP METABOLITES (Component 1) ###\n")
print(head(loadings_Y$value, 10))

# Save results
write_csv(as.data.frame(loadings_X$value), "species_loadings_comp1.csv")
write_csv(as.data.frame(loadings_Y$value), "metabolite_loadings_comp1.csv")

cat("\n✓ All done! Check your PDF and CSV files.\n")

# Update global environment
cat("\nUpdating global variables...\n")
assign("processed_data", processed_data_fixed, envir = .GlobalEnv)
assign("scca_result", result_scca, envir = .GlobalEnv)

cat("\n=== COMPLETE ===\n")
