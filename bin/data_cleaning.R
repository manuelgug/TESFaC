#!/usr/bin/env Rscript

# ========================================
# 1. DATA CLEANING (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse parameters from environment variables set by Nextflow
site <- Sys.getenv("NXF_SITE", "")
maf_filter <- as.numeric(Sys.getenv("NXF_MAF_FILTER", "0.01"))
min_allele_read_count <- as.numeric(Sys.getenv("NXF_MIN_ALLELE_READ_COUNT", "10"))
cum_curve_threshold <- as.numeric(Sys.getenv("NXF_CUM_CURVE_THRESHOLD", "0.999"))
OPTIM_AMPSET <- as.logical(Sys.getenv("NXF_OPTIM_AMPSET", "FALSE"))

cat("=== STARTING DATA CLEANING ===\n")

# -------------------- 1) IMPORT DATA --------------------
cat("Importing metadata...\n")

# Find input files
metadata_files <- list.files(".", pattern = "metadata_tes_.*\\.csv", full.names = TRUE)
genomic_site_files <- list.files(".", pattern = "genomic_site_.*\\.csv", full.names = TRUE)
madhito_files <- list.files(".", pattern = "pfPHAST_.*\\.csv", full.names = TRUE)

if (length(metadata_files) == 0) stop("No metadata file found")
if (length(genomic_site_files) == 0) stop("No genomic site file found")

metadata <- read.csv(metadata_files[1], 
                     stringsAsFactors = FALSE, 
                     colClasses = c(NIDA = "character"))

cat("Importing TES data from directories...\n")
tes_dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
tes_dirs <- tes_dirs[grepl("_RESULTS", tes_dirs)]
list_of_dfs <- list()

for (dir in tes_dirs) {
  folder_name <- gsub("_RESULTS.*", "", basename(dir))
  csv_file <- file.path(dir, "allele_data_global_max_0_filtered.csv")
  
  if (file.exists(csv_file)) {
    cat(paste0("Importing: ", dir, "\n"))
    df <- read.csv(csv_file)
    df$run <- folder_name
    list_of_dfs[[folder_name]] <- df
  }
}

merged_dfs <- bind_rows(list_of_dfs)

# -------------------- 2) FORMAT DATA --------------------
cat("Formatting sample IDs...\n")
merged_dfs <- merged_dfs %>%
  mutate(
    sampleID = gsub("N|_S.*$", "", sampleID),
    sampleID = if_else(str_detect(sampleID, "_"), sampleID, paste0(sampleID, ".0")),
    sampleID = gsub("_", ".", sampleID)
  )

metadata <- metadata %>%
  mutate(
    NIDA = gsub("N|_S.*$", "", NIDA),
    NIDA = if_else(!str_detect(NIDA, "\\."), paste0(NIDA, ".0"), NIDA)
  )

merged_dfs <- merged_dfs[merged_dfs$sampleID %in% metadata$NIDA,]
merged_dfs$sampleID <- paste0(merged_dfs$sampleID, "__", merged_dfs$run)
merged_dfs <- merged_dfs %>% select(sampleID, locus, pseudo_cigar, reads, norm.reads.locus, run)

# -------------------- 3) CLEAN TES DATA --------------------
cat("Cleaning TES data...\n")
clean_data <- function(df, maf_filter, min_read_count) {
  good_sampleID <- df %>%
    group_by(sampleID, locus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    filter(reads >= 100) %>%
    group_by(sampleID) %>%
    summarise(n_loci = n_distinct(locus), .groups = "drop") %>%
    filter(n_loci >= 50) %>%
    pull(sampleID)
  
  df <- df[df$sampleID %in% good_sampleID, ]
  
  high_read_samples <- df %>%
    group_by(sampleID) %>%
    summarise(total_reads = sum(reads), .groups = "drop") %>%
    filter(total_reads > 10000) %>%
    pull(sampleID)
  
  df <- df[df$sampleID %in% high_read_samples, ]
  df <- df[grepl("-1A$", df$locus), ]
  df$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", df$pseudo_cigar)
  df$pseudo_cigar[df$pseudo_cigar == "" | is.na(df$pseudo_cigar)] <- "."
  
  df <- df %>%
    group_by(sampleID, locus, pseudo_cigar) %>%
    summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus), .groups = "drop") %>%
    mutate(allele = paste0(locus, "__", pseudo_cigar)) %>%
    filter(!grepl("I=|D=", allele)) %>%
    filter(norm.reads.locus > maf_filter & reads > min_read_count) %>%
    select(-pseudo_cigar)
  
  return(df)
}

merged_dfs_agg <- clean_data(merged_dfs, maf_filter, min_allele_read_count)

# -------------------- 4) KEEP SAMPLES WITH PAIRS --------------------
cat("Filtering samples with pairs...\n")
ids <- data.frame(NIDA = unique(gsub("__.*$", "", merged_dfs_agg$sampleID)))
ids_merged <- merge(ids, metadata, by = "NIDA")
pairs_counts <- ids_merged %>% 
  group_by(PairsID) %>% 
  summarise(n_samples = n(), .groups = "drop") %>%
  filter(n_samples == 2)

metadata_updated <- metadata[metadata$PairsID %in% pairs_counts$PairsID, ]
merged_dfs_agg <- merged_dfs_agg[gsub("__.*$", "", merged_dfs_agg$sampleID) %in% metadata_updated$NIDA, ]

# -------------------- 5) CLEAN SITE DATA --------------------
cat("Cleaning site data...\n")
genomic_site <- read.csv(genomic_site_files[1], 
                        stringsAsFactors = FALSE, 
                        colClasses = c(NIDA2 = "character")) %>%
  rename(sampleID = NIDA2)

genomic_site_agg <- clean_data(genomic_site, maf_filter, min_allele_read_count)
genomic_site_agg <- genomic_site_agg[!sapply(genomic_site_agg$sampleID, function(x) any(grepl(x, merged_dfs_agg$sampleID, fixed = TRUE))), ]

# -------------------- 6) MERGE DATA --------------------
cat("Merging TES and site data...\n")
merged_dfs_agg$data_type <- "tes"
genomic_site_agg$data_type <- "site"
data_all <- rbind(merged_dfs_agg, genomic_site_agg)

# -------------------- 7) SHARED AMPLICONS --------------------
cat("Identifying shared amplicons...\n")
amps_data <- data_all %>%
  group_by(locus) %>%
  summarise(in_n_samples = n_distinct(sampleID), .groups = "drop") %>%
  arrange(desc(in_n_samples))

selected_amps <- amps_data$in_n_samples == max(amps_data$in_n_samples)
data_all <- data_all[data_all$locus %in% amps_data[selected_amps, ]$locus, ]

# -------------------- 8) SELECT AMPS --------------------
cat("Selecting amplicons...\n")
heterozygosity_per_sample <- data_all %>%
  group_by(sampleID, locus) %>%
  mutate(freq = norm.reads.locus / sum(norm.reads.locus, na.rm = TRUE)) %>%
  summarise(heterozygosity = 1 - sum(freq^2, na.rm = TRUE), .groups = "drop")

variability_per_locus <- heterozygosity_per_sample %>%
  group_by(locus) %>%
  summarise(mean_He = mean(heterozygosity, na.rm = TRUE), 
           sd_He = sd(heterozygosity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_He))

if (OPTIM_AMPSET) {
  cat("Using optimized amplicon set...\n")
  variability_per_locus <- variability_per_locus %>%
    mutate(prob_identity = 1 - mean_He)
  
  n_loci <- nrow(variability_per_locus)
  cumulative_He <- sapply(1:n_loci, function(i) 1 - prod(variability_per_locus$prob_identity[1:i]))
  cum_curve <- data.frame(loci_included = 1:n_loci, multilocus_He = cumulative_He)
  amps_var_data <- cbind(variability_per_locus, cum_curve)
  
  top_n_amps <- length(cum_curve$multilocus_He[cum_curve$multilocus_He < cum_curve_threshold])
  top_var_amps <- variability_per_locus$locus[1:top_n_amps]
  data_all <- data_all[data_all$locus %in% top_var_amps, ]
  
} else {
  cat("Using Madhito 20 amplicons...\n")
  if (length(madhito_files) > 0) {
    madhito_amps <- read.csv(madhito_files[1])
    data_all <- data_all[data_all$locus %in% madhito_amps$locus, ]
  }
}

# -------------------- 9) SUMMARY CHECKS --------------------
cat("=== DATA CLEANING SUMMARY ===\n")
cat("Total samples (failure pairs + site):", length(unique(data_all$sampleID)), "\n")
cat("Shared loci:", length(unique(data_all$locus)), "\n")
cat("Failure pairs:", length(unique(data_all[data_all$data_type == "tes", ]$sampleID)) / 2, "\n")

# -------------------- 10) SAVE OUTPUTS --------------------
cat("Saving cleaned data...\n")
write.csv(data_all, paste0("genomic_updated_", site, ".csv"), row.names = FALSE)
write.csv(metadata_updated, paste0("metadata_updated_", site, ".csv"), row.names = FALSE)

cat("=== DATA CLEANING COMPLETE ===\n\n")
