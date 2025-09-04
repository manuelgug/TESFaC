#!/usr/bin/env Rscript

# ========================================
# 7. CALCULATE FEATURES (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(reshape2)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(Matrix)
  library(progress)
  library(parallel)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")
chunk_size <- as.numeric(Sys.getenv("NXF_CHUNK_SIZE", "1000"))
num_cores <- as.numeric(Sys.getenv("NXF_NUM_CORES", "20"))

cat("=== STARTING FEATURE CALCULATION ===\n")

options(datatable.warn = FALSE)

output_suffixes <- c("TRAINING_DATA", "REAL_DATA")

# Helper: load & standardize both PAIRS_GENOMIC & PAIRS_METADATA
load_pairs <- function(type) {
  if (type == "TRAINING_DATA") {
    cat("Loading training data (synthetic pairs)...\n")
    
    pairs_metadata_files <- list.files(".", pattern = "PAIRS_METADATA_.*\\.RDS", full.names = TRUE)
    pairs_genomic_files <- list.files(".", pattern = "PAIRS_GENOMIC_.*\\.RDS", full.names = TRUE)
    
    if (length(pairs_metadata_files) == 0) stop("No pairs metadata file found")
    if (length(pairs_genomic_files) == 0) stop("No pairs genomic file found")
    
    meta <- readRDS(pairs_metadata_files[1])
    geno <- readRDS(pairs_genomic_files[1])
    
    setDT(meta); meta <- copy(meta)
    setDT(geno); geno <- copy(geno)
    geno <- merge(geno, meta, by = "PairsID", all.x = TRUE)
    
  } else if (type == "REAL_DATA") {
    cat("Loading real data (TES pairs)...\n")
    
    genomic_files <- list.files(".", pattern = "genomic_updated_.*\\.csv", full.names = TRUE)
    metadata_files <- list.files(".", pattern = "metadata_updated_.*\\.csv", full.names = TRUE)
    
    if (length(genomic_files) == 0) stop("No genomic updated file found")
    if (length(metadata_files) == 0) stop("No metadata updated file found")
    
    geno <- fread(genomic_files[1])
    meta <- fread(metadata_files[1], colClasses = c(NIDA = "character")) %>%
      filter(!is.na(time_point))
    
    setDT(geno); geno <- copy(geno)
    setDT(meta); meta <- copy(meta)
    
    # rename + split
    setnames(geno, c("reads", "sampleID"), c("read_counts", "NIDA"))
    geno[, c("NIDA", "run") := tstrsplit(NIDA, "__", fixed = TRUE)]
    
    # merge + derive pair_type
    geno <- merge(meta, geno, by = "NIDA", all = FALSE)
    coi_sum <- geno[
      time_point %in% c("D0","Dx"),
      .(
        offset_naive_coi_D0 = unique(offset_naive_coi[time_point=="D0"]),
        offset_naive_coi_Dx = unique(offset_naive_coi[time_point=="Dx"])
      ),
      by = PairsID
    ]
    coi_sum[, pair_type := paste0(offset_naive_coi_D0, "__", offset_naive_coi_Dx)]
    geno <- merge(geno, coi_sum[, .(PairsID, pair_type)], by = "PairsID", all.x = TRUE)
    
    # wide metadata, then convert to data.table and copy to avoid selfref warnings
    meta <- meta %>%
      select(PairsID, NIDA, time_point) %>%
      pivot_wider(
        names_from  = time_point,
        values_from = NIDA,
        names_prefix = "NIDA"
      ) %>%
      rename(NIDA1 = NIDAD0, NIDA2 = NIDADx) %>%
      left_join(coi_sum, by = "PairsID")
    setDT(meta); meta <- copy(meta)
    
  } else {
    stop("Unknown type: ", type)
  }
  list(META = meta, GENO = geno)
}

# Optimized feature calculation function
calculate_features_optimized <- function(sample1, sample2) {
  # 1) Unique alleles & loci
  alleles1 <- unique(sample1$allele)
  alleles2 <- unique(sample2$allele)
  all_alleles <- union(alleles1, alleles2)
  
  # 2) Allele-level intersection & union via set operations
  inter_cnt   <- length(intersect(alleles1, alleles2))
  union_cnt   <- length(all_alleles)
  jaccard     <- if (union_cnt > 0) inter_cnt / union_cnt else 0
  retention   <- if (length(alleles1) > 0) inter_cnt / length(alleles1) else 0
  allele_gain <- if (union_cnt > 0) length(setdiff(alleles2, alleles1)) / union_cnt else 0
  allele_loss <- if (union_cnt > 0) length(setdiff(alleles1, alleles2)) / union_cnt else 0
  
  # 3) Transition asymmetry
  trans_asym <- if ((allele_gain + allele_loss) > 0)
    (allele_gain - allele_loss) / (allele_gain + allele_loss) else 0
  
  # 4) Prepare locus-grouped allele lists
  split1 <- split(sample1$allele, sample1$locus)
  split2 <- split(sample2$allele, sample2$locus)
  loci   <- union(names(split1), names(split2))
  n_loci <- length(loci)
  
  # 5) Compute discordant loci count
  discordant_loci <- sum(vapply(
    loci,
    function(l) {
      length(intersect(split1[[l]] %||% character(0),
                       split2[[l]] %||% character(0))) == 0
    },
    logical(1)
  ))
  
  # 6) Compute replacement pattern sum
  replacement_pattern_sum <- sum(vapply(
    loci,
    function(l) {
      a1 <- split1[[l]] %||% character(0)
      a2 <- split2[[l]] %||% character(0)
      if      (length(a2) == 0)         1
      else if (all(a2 %in% a1))         0
      else                              length(setdiff(a2, a1)) / length(a2)
    },
    numeric(1)
  ))
  
  # 7) Final locus-level rates
  locus_discordance_rate    <- discordant_loci / n_loci
  replacement_pattern_score <- replacement_pattern_sum / n_loci
  
  # 8) Return all features
  c(
    jaccard_similarity          = jaccard,
    allele_retention_rate       = retention,
    allele_gain                 = allele_gain,
    allele_loss                 = allele_loss,
    locus_discordance_rate      = locus_discordance_rate,
    allele_transition_asymmetry = trans_asym,
    replacement_pattern_score   = replacement_pattern_score
  )
}

# Process both types, silencing self-ref warnings
suppressWarnings({
  for (TYPE in output_suffixes) {
    cat("=== Processing", TYPE, "===\n")
    
    parts           <- load_pairs(TYPE)
    PAIRS_METADATA  <- parts$META
    PAIRS_GENOMIC   <- parts$GENO
    
    # Pre-split genomic by NIDA
    cat("Pre-processing genomic data...\n")
    merged_list <- split(PAIRS_GENOMIC, by = "NIDA")
    merged_list <- lapply(merged_list, function(dt) {
      dt[, c("PairsID","time_point") := NULL]
      unique(dt)
    })
    gc()
    
    unique_pairs <- unique(PAIRS_METADATA$PairsID)
    pairs_to_d0  <- PAIRS_METADATA$NIDA1[match(unique_pairs, PAIRS_METADATA$PairsID)]
    pairs_to_dx  <- PAIRS_METADATA$NIDA2[match(unique_pairs, PAIRS_METADATA$PairsID)]
    
    out_file <- paste0("delta_features_", site, "_", TYPE, ".csv")
    
    # progress bar
    pb <- progress_bar$new(
      format = "  [:bar] :percent ETA: :eta",
      total  = ceiling(length(unique_pairs) / chunk_size),
      width  = 60
    )
    
    first_chunk <- TRUE
    for (i in seq(1, length(unique_pairs), by = chunk_size)) {
      idx <- i:min(i + chunk_size - 1, length(unique_pairs))
      cat(" chunk", i, "–", idx[length(idx)], "\n")
      
      results <- mclapply(idx, function(k) {
        s1 <- merged_list[[ pairs_to_d0[k] ]]
        s2 <- merged_list[[ pairs_to_dx[k] ]]
        feats <- calculate_features_optimized(s1, s2)
        data.frame(t(feats), PairsID = unique_pairs[k], check.names = FALSE)
      }, mc.cores = num_cores)
      
      dt_chunk <- rbindlist(Filter(Negate(is.null), results), fill = TRUE)
      fwrite(dt_chunk, out_file, append = !first_chunk)
      first_chunk <- FALSE
      gc()
      pb$tick()
    }
    
    # Final merge back
    cat("Finalizing feature calculations...\n")
    delta_df <- fread(out_file) %>%
      select(PairsID, everything()) %>%
      left_join(PAIRS_METADATA, by = "PairsID")
    
    if (TYPE == "REAL_DATA") {
      delta_df[, coi_change := offset_naive_coi_Dx - offset_naive_coi_D0]
    } else {
      delta_df[, coi_change := Dx_nstrains - D0_nstrains]
    }
    
    write.csv(delta_df,
              paste0(site, "_", TYPE, ".csv"),
              row.names = FALSE)
    
    cat("=> Feature calculation complete:", out_file, "\n\n")
  }
})

cat("=== FEATURE CALCULATION COMPLETE ===\n\n")
