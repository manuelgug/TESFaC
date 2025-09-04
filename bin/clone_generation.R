#!/usr/bin/env Rscript

# ========================================
# 4. GENERATE CLONE DATA (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(data.table)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")
N_CLONES_TARGET <- as.numeric(Sys.getenv("NXF_N_CLONES_TARGET", "1000"))

cat("=== STARTING CLONE GENERATION ===\n")

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

# Join time_point info
data <- left_join(data, metadata_updated[c("NIDA", "time_point", "offset_naive_coi")], 
                  by = c("sampleID" = "NIDA"))

# Handle missing time_point values
data$time_point <- ifelse(is.na(data$time_point), "D0", data$time_point)

cat("Identifying natural clones...\n")

# Identify samples with only 1 allele per locus (true clones)
clones <- data %>%
  group_by(sampleID, locus) %>%
  summarise(n = n_distinct(allele), .groups = "drop") %>%
  group_by(sampleID) %>%
  summarise(all_single = all(n == 1), .groups = "drop") %>%
  filter(all_single) %>%
  pull(sampleID)

clones_genomic <- data %>%
  filter(sampleID %in% clones)

# Handle missing time_point values
clones_genomic$time_point <- ifelse(is.na(clones_genomic$time_point), "D0", clones_genomic$time_point)
clones_genomic <- clones_genomic[clones_genomic$time_point == "D0",] # only D0

cat("Natural clones found:", length(unique(clones_genomic$sampleID)), "\n")

### 2) Create artificial clones from polyclonal infections ----
cat("Generating synthetic clones...\n")

N_CLONES <- N_CLONES_TARGET - length(clones)
data_polyclonal <- data %>% filter(!sampleID %in% clones & time_point == "D0") # only D0
polyclonal_samples <- unique(data_polyclonal$sampleID)
iterations <- ceiling(N_CLONES / length(polyclonal_samples))

# Generate synthetic clones
clone_result <- data.frame()

for (sid in polyclonal_samples) {
  sample_data <- data_polyclonal[data_polyclonal$sampleID == sid, ]
  
  for (i in 1:iterations) {
    # Sample one allele per locus using norm.reads.locus as weights
    sampled <- sample_data %>%
      group_by(locus) %>%
      slice_sample(n = 1, weight_by = norm.reads.locus, replace = TRUE) %>%
      ungroup()
    
    # Add unique clone ID
    sampled$sampleID <- paste0(sid, "__clone_", i)
    
    # Append to result
    clone_result <- bind_rows(clone_result, sampled)
  }
}

all_clones <- bind_rows(clone_result, clones_genomic)
cat("Total clones before filtering:", length(unique(all_clones$sampleID)), "\n")

### 3) Calculate pairwise proportion of shared alleles ----
cat("Calculating pairwise similarity to remove redundant clones...\n")

alleles <- all_clones %>%
  group_by(sampleID) %>%
  summarize(alleles = list(allele), .groups = "drop")

n <- nrow(alleles)
comparison_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(alleles$sampleID, alleles$sampleID))

# Precompute all unique index combinations (i < j)
pair_indices <- combn(n, 2)

# Function to compute Jaccard similarity
compute_shared_prop <- function(i, j) {
  a_i <- alleles$alleles[[i]]
  a_j <- alleles$alleles[[j]]
  shared <- length(intersect(a_i, a_j))
  total <- length(unique(c(a_i, a_j)))
  shared / total
}

# Apply function across all combinations
shared_values <- apply(pair_indices, 2, function(x) compute_shared_prop(x[1], x[2]))

# Fill the upper and lower triangle
for (k in seq_along(shared_values)) {
  i <- pair_indices[1, k]
  j <- pair_indices[2, k]
  comparison_matrix[i, j] <- shared_values[k]
  comparison_matrix[j, i] <- shared_values[k]
}

# Fill diagonal with 1s
diag(comparison_matrix) <- 1

comparison_df <- as.data.frame(comparison_matrix) %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  filter(Var1 < Var2) %>%
  arrange(value)

### 4) Clean monoclonal data: remove redundant clone groups ----
cat("Removing redundant clones...\n")

same_clones <- comparison_df %>% filter(value > 0.99)
groups <- list()

for (i in seq_len(nrow(same_clones))) {
  v1 <- same_clones$Var1[i]
  v2 <- same_clones$Var2[i]
  matched <- FALSE
  
  for (j in seq_along(groups)) {
    if (v1 %in% groups[[j]] || v2 %in% groups[[j]]) {
      groups[[j]] <- unique(c(groups[[j]], v1, v2))
      matched <- TRUE
      break
    }
  }
  
  if (!matched) {
    groups <- append(groups, list(c(v1, v2)))
  }
}

# Keep only one representative per group
clones_to_remove <- unlist(lapply(groups, function(g) g[-1]))
all_clones <- all_clones %>% filter(!sampleID %in% clones_to_remove)

cat("Final unique clones after filtering:", length(unique(all_clones$sampleID)), "\n")

### 5) Export cleaned clone data ----
write.csv(all_clones, paste0("clones_genomic_data_", site, ".csv"), row.names = FALSE)

cat("=== CLONE GENERATION COMPLETE ===\n")
cat("Final clone count:", length(unique(all_clones$sampleID)), "\n\n")
