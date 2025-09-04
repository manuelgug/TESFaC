#!/usr/bin/env Rscript

# ========================================
# 3. FNR CALCULATION (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")
FNR_02_03 <- as.numeric(Sys.getenv("NXF_FNR_02_03", "0.11"))
FNR_below_02 <- as.numeric(Sys.getenv("NXF_FNR_BELOW_02", "0.16"))

cat("=== STARTING FNR CALCULATION ===\n")

# Load data
metadata_files <- list.files(".", pattern = "metadata_updated_.*\\.csv", full.names = TRUE)
genomic_files <- list.files(".", pattern = "genomic_updated_.*\\.csv", full.names = TRUE)

if (length(metadata_files) == 0) stop("No metadata updated file found")
if (length(genomic_files) == 0) stop("No genomic updated file found")

metadata_updated <- read.csv(metadata_files[1], 
                             stringsAsFactors = FALSE, 
                             colClasses = c(NIDA = "character"))

data <- read.csv(genomic_files[1],
                 stringsAsFactors = FALSE,
                 colClasses = c(sampleID = "character")) %>%
  mutate(sampleID = gsub("__.*", "", sampleID))

cat("Calculating strain proportions...\n")

# Handle missing time_point values
metadata_updated$time_point <- ifelse(is.na(metadata_updated$time_point), "D0", metadata_updated$time_point)

# Calculate allele count per locus
allele_count_per_locus <- data %>% 
  group_by(sampleID, locus) %>%
  summarise(n_alleles = length(unique(allele)), .groups = "drop")

# Merge with COI data
merged_df <- left_join(allele_count_per_locus, 
                      metadata_updated[c("NIDA", "offset_naive_coi")], 
                      by = c("sampleID" = "NIDA"))

# Filter rows where COI matches n_alleles
matched_df <- merged_df %>%
  filter(n_alleles == offset_naive_coi)

# Subset data by matching NIDA and locus
subset_labcontrols_genomic <- semi_join(data, matched_df, by = c("sampleID", "locus"))
subset_labcontrols_genomic <- unique(subset_labcontrols_genomic)

# Rank norm.reads.locus descendingly within each NIDA and locus
ranked_df <- subset_labcontrols_genomic %>%
  group_by(sampleID, locus) %>%
  arrange(desc(norm.reads.locus), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Compute the average of norm.reads.locus by rank per NIDA
average_by_rank <- ranked_df %>%
  group_by(sampleID, rank) %>%
  summarise(avg_norm_reads = mean(norm.reads.locus, na.rm = TRUE), .groups = "drop") %>%
  arrange(sampleID, rank)

# Convert to wide format
average_by_rank_wide <- average_by_rank %>%
  pivot_wider(
    names_from = rank,
    values_from = avg_norm_reads,
    names_prefix = "rank_"
  )

cat("Calculating FNR values based on strain proportions...\n")

# Calculate FNRs based on strain proportions
FNRs <- average_by_rank_wide %>%
  rowwise() %>%
  mutate(FNR = list({
    ranks <- c_across(starts_with("rank_"))
    fnrs <- c()
    # FNR for strains with abundance 0.02-0.03
    fnrs <- c(fnrs, rep(FNR_02_03, sum(ranks >= 0.02 & ranks < 0.03, na.rm = TRUE)))
    # FNR for strains below 0.02  
    fnrs <- c(fnrs, rep(FNR_below_02, sum(ranks < 0.02, na.rm = TRUE)))
    fnrs
  })) %>%
  ungroup()

# Merge with COI data
FNRs <- left_join(FNRs, metadata_updated[c("NIDA", "offset_naive_coi")], 
                  by = c("sampleID" = "NIDA"))

# Get final FNR values
FNRs <- FNRs %>%
  distinct(offset_naive_coi, FNR)

# Remove empty vectors or COI = 1 (no FNRs applied)
FNRs <- FNRs %>% 
  filter(lengths(FNR) > 0) %>%
  filter(offset_naive_coi > 1)

# Save FNR data
saveRDS(FNRs, paste0("FNRs_", site, ".RDS"))

cat("=== FNR CALCULATION COMPLETE ===\n")
cat("FNR summary:\n")
if(nrow(FNRs) > 0) {
  print(FNRs)
} else {
  cat("No FNR values calculated - all strain proportions above threshold\n")
}
cat("\n")
