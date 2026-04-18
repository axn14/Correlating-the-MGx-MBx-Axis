
# CORRECTED FIX (primary correction)


library(tidyverse)
library(mixOmics)

# The issue: Your data has SAMPLES in columns, not rows
# We need to handle this differently

preprocess_omics_data_corrected <- function(data_obj, min_prevalence = 0.1) {
  
  cat("\n=== Preprocessing Data (CORRECTED) ===\n")
  
  # Row 1 = Species/Metabolite ID
  # Columns 2-end = Sample_1, Sample_2, etc.
  
  species_df <- data_obj$species
  mtb_df <- data_obj$metabolites
  
  cat(paste("Original species dimensions:", nrow(species_df), "x", ncol(species_df), "\n"))
  cat(paste("Original metabolite dimensions:", nrow(mtb_df), "x", ncol(mtb_df), "\n"))
  
  # Extract feature names
  species_names <- species_df$Sample
  mtb_names <- mtb_df$Sample
  
  # Extract data (remove first column which is feature IDs)
  species_data <- species_df %>% select(-Sample) %>% as.matrix()
  mtb_data <- mtb_df %>% select(-Sample) %>% as.matrix()
  
  # Set row names
  rownames(species_data) <- species_names
  rownames(mtb_data) <- mtb_names
  
  # Get sample names from column names
  sample_names <- colnames(species_data)
  
  cat(paste("\nNumber of samples:", length(sample_names), "\n"))
  cat(paste("Number of species:", nrow(species_data), "\n"))
  cat(paste("Number of metabolites:", nrow(mtb_data), "\n"))
  
  # transpose so samples are rows and features are columns
  species_mat <- t(species_data)
  mtb_mat <- t(mtb_data)
  
  cat(paste("\nAfter transpose:\n"))
  cat(paste("  Species:", nrow(species_mat), "samples x", ncol(species_mat), "features\n"))
  cat(paste("  Metabolites:", nrow(mtb_mat), "samples x", ncol(mtb_mat), "features\n"))
  
  # Verify sample names match
  if (!all(rownames(species_mat) == rownames(mtb_mat))) {
    stop("Sample names don't match between species and metabolites!")
  }
  
  # Filter by prevalence (features present in at least X% of samples)
  species_prev <- colSums(species_mat > 0) / nrow(species_mat)
  mtb_prev <- colSums(mtb_mat > 0) / nrow(mtb_mat)
  
  cat(paste("\nPrevalence filtering (>", min_prevalence*100, "%):\n", sep=""))
  cat(paste("  Species retained:", sum(species_prev >= min_prevalence), "/", ncol(species_mat), "\n"))
  cat(paste("  Metabolites retained:", sum(mtb_prev >= min_prevalence), "/", ncol(mtb_mat), "\n"))
  
  species_mat <- species_mat[, species_prev >= min_prevalence, drop = FALSE]
  mtb_mat <- mtb_mat[, mtb_prev >= min_prevalence, drop = FALSE]
  
  # Log transformation
  cat("\nLog-transforming...\n")
  species_log <- log10(species_mat + 1e-6)
  mtb_log <- log10(mtb_mat + 1e-6)
  
  # Center and scale
  cat("Scaling...\n")
  species_scaled <- scale(species_log, center = TRUE, scale = TRUE)
  mtb_scaled <- scale(mtb_log, center = TRUE, scale = TRUE)
  
  # Clean up any NAs or Inf values
  species_scaled[is.na(species_scaled) | is.infinite(species_scaled)] <- 0
  mtb_scaled[is.na(mtb_scaled) | is.infinite(mtb_scaled)] <- 0
  
  # Convert to plain matrix (remove scale attributes)
  species_scaled <- matrix(species_scaled, 
                           nrow = nrow(species_scaled), 
                           ncol = ncol(species_scaled),
                           dimnames = list(rownames(species_mat), 
                                         colnames(species_mat)[species_prev >= min_prevalence]))
  
  mtb_scaled <- matrix(mtb_scaled, 
                      nrow = nrow(mtb_scaled), 
                      ncol = ncol(mtb_scaled),
                      dimnames = list(rownames(mtb_mat), 
                                    colnames(mtb_mat)[mtb_prev >= min_prevalence]))
  
  cat("\n✓ Preprocessing complete!\n")
  cat(paste("  Final dimensions:\n"))
  cat(paste("    Species:", nrow(species_scaled), "samples x", ncol(species_scaled), "features\n"))
  cat(paste("    Metabolites:", nrow(mtb_scaled), "samples x", ncol(mtb_scaled), "features\n"))
  cat(paste("  Data ranges:\n"))
  cat(paste("    Species: [", round(min(species_scaled),2), ",", round(max(species_scaled),2), "]\n"))
  cat(paste("    Metabolites: [", round(min(mtb_scaled),2), ",", round(max(mtb_scaled),2), "]\n"))
  
  return(list(
    species = species_scaled,
    metabolites = mtb_scaled,
    metadata = data_obj$metadata,
    mtb_map = data_obj$mtb_map
  ))
}


# RUN THE CORRECTED PREPROCESSING


cat("\n########## RUNNING CORRECTED ANALYSIS ##########\n")

# Preprocess with corrected function
processed_data <- preprocess_omics_data_corrected(data, min_prevalence = 0.05)

# Verify dimensions match
cat("\n### Verifying data ###\n")
cat(paste("X (species):", nrow(processed_data$species), "x", ncol(processed_data$species), "\n"))
cat(paste("Y (metabolites):", nrow(processed_data$metabolites), "x", ncol(processed_data$metabolites), "\n"))
cat(paste("Sample names match:", all(rownames(processed_data$species) == rownames(processed_data$metabolites)), "\n"))

# Run mixOmics
cat("\n### Running mixOmics sPLS ###\n")

X <- processed_data$species
Y <- processed_data$metabolites

# Adaptive parameters based on your data size
ncomp <- min(3, floor(nrow(X) / 10))  # Conservative
keepX <- min(50, floor(ncol(X) * 0.2))  # Keep 20% of species features
keepY <- min(30, ncol(Y))  # All metabolites if < 30, otherwise 30

cat(paste("Parameters:\n"))
cat(paste("  ncomp:", ncomp, "\n"))
cat(paste("  keepX:", keepX, "(selecting from", ncol(X), "species)\n"))
cat(paste("  keepY:", keepY, "(selecting from", ncol(Y), "metabolites)\n"))

# Run sPLS
result_scca <- spls(
  X = X,
  Y = Y,
  ncomp = ncomp,
  mode = "regression",
  keepX = rep(keepX, ncomp),
  keepY = rep(keepY, ncomp)
)

cat("\n✓ sPLS completed successfully!\n")

# Generate plots
cat("\nGenerating visualizations...\n")
pdf("mixomics_plots.pdf", width = 14, height = 10)

# Plot 1: Sample plot
plotIndiv(result_scca, 
          comp = c(1, 2),
          ind.names = FALSE,
          title = "sPLS - Sample Space",
          pch = 16,
          cex = 1.5)

# Plot 2: Variable correlation plot
plotVar(result_scca, 
        comp = c(1, 2),
        cex = c(1, 1),
        title = "sPLS - Variable Relationships",
        cutoff = 0.5)

# Plot 3: Arrow plot
plotArrow(result_scca, 
          comp = c(1, 2),
          ind.names = FALSE,
          title = "sPLS - Sample Correlations Between Omics")

# Plot 4: Individual plots for X and Y
par(mfrow = c(1, 2))
plotLoadings(result_scca, comp = 1, contrib = 'max', method = 'mean', 
             title = "Species Loadings - Component 1")
plotLoadings(result_scca, comp = 1, contrib = 'max', method = 'mean', 
             block = 'Y', title = "Metabolite Loadings - Component 1")

dev.off()

cat("✓ Plots saved to: mixomics_plots.pdf\n")

# Extract and save loadings
cat("\n### Extracting important features ###\n")

loadings_species <- selectVar(result_scca, comp = 1)
loadings_metabolites <- selectVar(result_scca, comp = 1, block = "Y")

# Create nice dataframes
species_loadings_df <- data.frame(
  species = rownames(loadings_species$value),
  loading = loadings_species$value$value.var,
  abs_loading = abs(loadings_species$value$value.var)
) %>%
  arrange(desc(abs_loading))

metabolite_loadings_df <- data.frame(
  metabolite = rownames(loadings_metabolites$value),
  loading = loadings_metabolites$value$value.var,
  abs_loading = abs(loadings_metabolites$value$value.var)
) %>%
  arrange(desc(abs_loading))

# Save to CSV
write_csv(species_loadings_df, "species_loadings_component1.csv")
write_csv(metabolite_loadings_df, "metabolite_loadings_component1.csv")

cat("\n### TOP 10 SPECIES (Component 1) ###\n")
print(head(species_loadings_df, 10))

cat("\n### TOP 10 METABOLITES (Component 1) ###\n")
print(head(metabolite_loadings_df, 10))

# Get correlations between selected variables
cat("\n### Extracting microbe-metabolite correlations ###\n")

# Get selected features from component 1
selected_species <- species_loadings_df$species[1:min(20, nrow(species_loadings_df))]
selected_metabolites <- metabolite_loadings_df$metabolite[1:min(20, nrow(metabolite_loadings_df))]

# Calculate correlations
cor_matrix <- cor(X[, selected_species], Y[, selected_metabolites], method = "spearman")

# Convert to long format
cor_df <- as.data.frame(as.table(cor_matrix))
names(cor_df) <- c("species", "metabolite", "correlation")
cor_df <- cor_df %>%
  mutate(abs_cor = abs(correlation)) %>%
  arrange(desc(abs_cor))

write_csv(cor_df, "top_microbe_metabolite_correlations.csv")

cat("\n### TOP 10 MICROBE-METABOLITE CORRELATIONS ###\n")
print(head(cor_df, 10))

# Update global environment
assign("processed_data", processed_data, envir = .GlobalEnv)
assign("scca_result", result_scca, envir = .GlobalEnv)

cat("\n")
cat("################################################################################\n")
cat("# ANALYSIS COMPLETE!\n")
cat("################################################################################\n")
cat("\nFiles generated:\n")
cat("  1. mixomics_plots.pdf - All visualizations\n")
cat("  2. species_loadings_component1.csv - Top species ranked by importance\n")
cat("  3. metabolite_loadings_component1.csv - Top metabolites ranked by importance\n")
cat("  4. top_microbe_metabolite_correlations.csv - Pairwise correlations\n")
cat("\nGlobal variables updated:\n")
cat("  - processed_data (cleaned matrices)\n")
cat("  - scca_result (mixOmics result object)\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Examine the plots in mixomics_plots.pdf\n")
cat("  2. Review top_microbe_metabolite_correlations.csv for key associations\n")
cat("  3. Map significant metabolites to KEGG pathways\n")
cat("  4. Run correlation network analysis for more associations:\n")
cat("     source('microbiome_metabolome_integration.R')\n")
cat("     network_result <- run_correlation_network(processed_data)\n")
cat("\n################################################################################\n")
