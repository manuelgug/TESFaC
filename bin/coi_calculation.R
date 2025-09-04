#!/usr/bin/env Rscript

# ========================================
# 2. OFFSET NAIVE COI CALCULATION (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")

cat("=== STARTING COI CALCULATION ===\n")

# Load cleaned data
data_all_files <- list.files(".", pattern = "genomic_updated_.*\\.csv", full.names = TRUE)
metadata_files <- list.files(".", pattern = "metadata_updated_.*\\.csv", full.names = TRUE)

if (length(data_all_files) == 0) stop("No genomic updated file found")
if (length(metadata_files) == 0) stop("No metadata updated file found")

data_all <- read.csv(data_all_files[1], stringsAsFactors = FALSE)
metadata_updated <- read.csv(metadata_files[1], stringsAsFactors = FALSE, 
                            colClasses = c(NIDA = "character"))

# Calculate offset naive COI directly, without moire
cat("Calculating offset naive COI...\n")
coi_stats <- data_all %>%
  group_by(sampleID, locus) %>%
  summarise(n_alleles = n_distinct(allele), .groups = "drop") %>%
  group_by(sampleID) %>%
  summarise(
    offset_naive_coi = {
      vals <- sort(n_alleles, decreasing = TRUE)
      if (length(vals) >= 2) vals[2] else vals[1] # OFFSET = THE 2ND VALUE
    }
  ) %>% 
  rename(NIDA = sampleID)

# Update metadata with COI
cat("Updating metadata with COI values...\n")
coi_stats$NIDA <- gsub("__.*", "", coi_stats$NIDA)
metadata_updated2 <- merge(metadata_updated, coi_stats, by = c("NIDA"))
metadata_updated2 <- metadata_updated2 %>% arrange(SampleID)

# Save updated metadata (keep genomic data unchanged)
write.csv(data_all, paste0("genomic_updated_", site, ".csv"), row.names = FALSE)
write.csv(metadata_updated2, paste0("metadata_updated_", site, ".csv"), row.names = FALSE)

cat("=== COI CALCULATION COMPLETE ===\n")
cat("COI statistics summary:\n")
print(summary(coi_stats$offset_naive_coi))
cat("\n")
