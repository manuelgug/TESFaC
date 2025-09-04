#!/usr/bin/env Rscript

# ========================================
# 6. GENERATE PAIRS DATA (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")

cat("=== STARTING PAIRS GENERATION ===\n")

# Load data
metadata_files <- list.files(".", pattern = "metadata_updated_.*\\.csv", full.names = TRUE)
mixes_metadata_files <- list.files(".", pattern = "MIXES_METADATA_.*\\.RDS", full.names = TRUE)
mixes_genomic_files <- list.files(".", pattern = "MIXES_GENOMIC_.*\\.RDS", full.names = TRUE)

if (length(metadata_files) == 0) stop("No metadata updated file found")
if (length(mixes_metadata_files) == 0) stop("No mixes metadata file found")
if (length(mixes_genomic_files) == 0) stop("No mixes genomic file found")

metadata_updated <- read.csv(metadata_files[1], 
                             stringsAsFactors = FALSE, 
                             colClasses = c(NIDA = "character"))

metadata_updated <- metadata_updated[!is.na(metadata_updated$PairsID),] # only tes data

MIXES_METADATA <- readRDS(mixes_metadata_files[1])
MIXES_GENOMIC <- readRDS(mixes_genomic_files[1])

#------------------------------------------------
# 1. Determine pairs of mixes based on COIs
#------------------------------------------------
cat("Determining pair combinations based on COI values...\n")

metadata_updated$offset_naive_coi <- round(metadata_updated$offset_naive_coi)

metadata_updated_wide <- metadata_updated %>%
  pivot_wider(
    id_cols = PairsID,
    names_from = time_point,
    values_from = c(NIDA, offset_naive_coi),
    names_glue = "{.value}_{time_point}"
  )

unique_combos <- metadata_updated_wide %>%
  distinct(offset_naive_coi_D0, offset_naive_coi_Dx) %>%
  mutate(
    offset_naive_coi_D0 = paste0("mix", offset_naive_coi_D0),
    offset_naive_coi_Dx = paste0("mix", offset_naive_coi_Dx)
  ) %>%
  arrange(offset_naive_coi_D0, offset_naive_coi_Dx)

#------------------------------------------------
# 2. Create pairs
#------------------------------------------------
cat("Creating pairs from mixture data...\n")

nidas_all <- MIXES_METADATA$NIDA
pairs_df <- expand.grid(NIDA1 = nidas_all, NIDA2 = nidas_all, stringsAsFactors = FALSE)

pairs_df <- pairs_df[rowSums(sapply(1:nrow(unique_combos), function(i) {
  grepl(unique_combos$offset_naive_coi_D0[i], pairs_df$NIDA1) &
    grepl(unique_combos$offset_naive_coi_Dx[i], pairs_df$NIDA2)
})) > 0, ]

pairs_df <- pairs_df %>%
  mutate(PairsID = row_number()) %>%
  pivot_longer(cols = c(NIDA1, NIDA2), names_to = "time_point", values_to = "NIDA") %>%
  mutate(time_point = ifelse(time_point == "NIDA1", "D0", "Dx")) %>%
  select(PairsID, NIDA, time_point)

colnames(MIXES_GENOMIC)[colnames(MIXES_GENOMIC) == "mixID"] <- "NIDA"
merged_dfs <- left_join(pairs_df, MIXES_GENOMIC, by = "NIDA", relationship = "many-to-many")

cat("Total pairs created:", length(unique(merged_dfs$PairsID)), "\n")

saveRDS(merged_dfs, paste0("PAIRS_GENOMIC_", site, ".RDS"))

#------------------------------------------------
# 3. Compare allele content of pairs
#------------------------------------------------
cat("Analyzing allele sharing between pairs...\n")

alleles <- merged_dfs %>%
  group_by(PairsID, NIDA, time_point) %>%
  summarize(alleles = list(allele), .groups = "drop") %>%
  arrange(PairsID, time_point)

alleles_shared_prop <- alleles %>%
  group_by(PairsID) %>%
  summarize(
    NIDA1 = NIDA[1],
    NIDA2 = NIDA[2],
    shared_count = length(intersect(alleles[[1]], alleles[[2]])),
    union_count = length(union(alleles[[1]], alleles[[2]])),
    shared_prop = shared_count / union_count,
    .groups = "drop"
  ) %>%
  mutate(
    NIDA1_trimmed = sub("_.*", "", NIDA1),
    NIDA2_trimmed = sub("_.*", "", NIDA2),
    pair_type = paste0(NIDA1_trimmed, "_", NIDA2_trimmed)
  ) %>%
  select(-NIDA1_trimmed, -NIDA2_trimmed)

#------------------------------------------------
# 4. Label data
#------------------------------------------------
cat("Labeling pairs as recrudescent or new infection...\n")

MIXES_METADATA <- MIXES_METADATA %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  rowwise() %>%
  mutate(
    STRAINS = list(na.omit(c_across(matches("^strain\\d?$"))))
  ) %>%
  ungroup()

PAIRS <- merge(pairs_df, MIXES_METADATA[c("NIDA", "STRAINS")], by = "NIDA") %>%
  arrange(PairsID, time_point)

labels <- PAIRS %>%
  group_by(PairsID) %>%
  summarise(
    labels = ifelse(
      any(unlist(STRAINS[time_point == "Dx"]) %in% unlist(STRAINS[time_point == "D0"])),
      "R", "NI"
    )
  )

PAIRS <- PAIRS %>%
  rowwise() %>%
  mutate(
    nstrains = length(STRAINS),
    sample = paste(STRAINS, collapse = "_")
  )

PAIRS_wide <- PAIRS %>%
  pivot_wider(
    names_from = time_point,
    values_from = c(sample, nstrains),
    names_glue = "{time_point}_{.value}"
  )

PAIRS_metadata <- PAIRS_wide %>%
  group_by(PairsID) %>%
  select(-NIDA, -STRAINS) %>%
  reframe(
    D0_sample = first(na.omit(D0_sample)),
    Dx_sample = first(na.omit(Dx_sample)),
    D0_nstrains = first(na.omit(D0_nstrains)),
    Dx_nstrains = first(na.omit(Dx_nstrains))
  )

PAIRS_metadata <- inner_join(PAIRS_metadata, labels, by = "PairsID")
PAIRS_metadata <- inner_join(PAIRS_metadata, alleles_shared_prop, by = "PairsID")
PAIRS_metadata$pair_type <- paste0(PAIRS_metadata$D0_nstrains, "__", PAIRS_metadata$Dx_nstrains)

saveRDS(PAIRS_metadata, paste0("PAIRS_METADATA_", site, ".RDS"))

#------------------------------------------------
# 5. Pair metadata mini EDA
#------------------------------------------------
cat("Creating pairs summary statistics...\n")

PAIRS_summary <- PAIRS_metadata %>%
  group_by(pair_type) %>%
  summarise(
    n_pairs = n(),
    NI_prop = mean(labels == "NI"),
    R_prop = mean(labels == "R"),
    NI_size = NI_prop * n_pairs,
    R_size = R_prop * n_pairs
  )

cat("Pairs summary by type:\n")
print(PAIRS_summary)

write.csv(PAIRS_summary, paste0("PAIRS_SUMMARY_", site, ".csv"), row.names = FALSE)

cat("=== PAIRS GENERATION COMPLETE ===\n")
cat("Total unique pairs:", length(unique(PAIRS_metadata$PairsID)), "\n\n")
