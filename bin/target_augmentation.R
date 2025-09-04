#!/usr/bin/env Rscript

# ========================================
# 8. TARGET AUGMENTATION (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(data.table)
  library(progress)
  library(parallel)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")
clone_seed <- as.numeric(Sys.getenv("NXF_CLONE_SEED", "420"))
pairs_seed <- as.numeric(Sys.getenv("NXF_PAIRS_SEED", "42069"))
chunk_size <- as.numeric(Sys.getenv("NXF_CHUNK_SIZE", "1000"))
num_cores <- as.numeric(Sys.getenv("NXF_NUM_CORES", "1"))

cat("=== STARTING TARGET AUGMENTATION ===\n")

# Load FNR data
fnr_files <- list.files(".", pattern = "FNRs_.*\\.RDS", full.names = TRUE)

if (length(fnr_files) == 0) {
  cat("No FNR file found. Skipping target augmentation.\n")
  
  # Just copy the training data without augmentation
  training_files <- list.files(".", pattern = ".*_TRAINING_DATA\\.csv", full.names = TRUE)
  if (length(training_files) > 0) {
    training_data <- read.csv(training_files[1])
    write.csv(training_data, paste0(site, "_TRAINING_DATA.csv"), row.names = FALSE)
  }
  
  cat("=== TARGET AUGMENTATION SKIPPED ===\n\n")
  quit(save = "no", status = 0)
}

FNRs <- readRDS(fnr_files[1])

if (nrow(FNRs) > 0) {
  cat("FNR data found. Proceeding with target augmentation...\n")
  
  # Load other required files
  pairs_metadata_files <- list.files(".", pattern = "PAIRS_METADATA_.*\\.RDS", full.names = TRUE)
  pairs_genomic_files <- list.files(".", pattern = "PAIRS_GENOMIC_.*\\.RDS", full.names = TRUE)
  clones_files <- list.files(".", pattern = "clones_genomic_data_.*\\.csv", full.names = TRUE)
  training_files <- list.files(".", pattern = ".*_TRAINING_DATA\\.csv", full.names = TRUE)
  
  if (length(pairs_metadata_files) == 0) stop("No pairs metadata file found")
  if (length(pairs_genomic_files) == 0) stop("No pairs genomic file found")
  if (length(clones_files) == 0) stop("No clones genomic data file found")
  if (length(training_files) == 0) stop("No training data file found")
  
  PAIRS_METADATA <- readRDS(pairs_metadata_files[1])
  PAIRS_GENOMIC <- readRDS(pairs_genomic_files[1])
  clones_genomic <- read.csv(clones_files[1], stringsAsFactors = FALSE)
  TRUTH_CLEAN <- read.csv(training_files[1])
  
  ## EXTRACT PAIRS THAT WILL BE FNRed
  cat("Identifying pairs for FNR application...\n")
  FNR_pairs <- PAIRS_METADATA[PAIRS_METADATA$D0_nstrains %in% FNRs$offset_naive_coi | 
                               PAIRS_METADATA$Dx_nstrains %in% FNRs$offset_naive_coi ,]$PairsID
  FNR_pairs_GENOMIC <- PAIRS_GENOMIC[PAIRS_GENOMIC$PairsID %in% FNR_pairs,]
  FNR_pairs_METADATA <- PAIRS_METADATA[PAIRS_METADATA$PairsID %in% FNR_pairs,]
  
  FNR_pairs_METADATA <- FNR_pairs_METADATA %>%
    select(PairsID, D0_sample, Dx_sample, D0_nstrains, Dx_nstrains) %>%
    pivot_longer(
      cols = c(D0_sample, Dx_sample, D0_nstrains, Dx_nstrains),
      names_to = c("time_point", ".value"),
      names_sep = "_"
    ) %>%
    rename(sample = sample, n_strains = nstrains)
  
  FNR_pairs_METADATA <- FNR_pairs_METADATA %>%
    mutate(sample = str_split(sample, "(?<=\\d)_(?=\\d)"))
  
  # Load clone data
  TRUTH <- clones_genomic
  colnames(TRUTH)[colnames(TRUTH) == "sampleID"] <- "sampleID"
  
  # Get unique FNR values
  unique_fnr_values <- unique(unlist(FNRs$FNR))
  
  set.seed(clone_seed)
  
  cat("Applying FNR to clone data...\n")
  
  # Apply FNR to TRUTH
  TRUTH_filtered_list <- map(unique_fnr_values, function(fnr) {
    TRUTH %>%
      group_by(sampleID) %>%
      group_modify(~ {
        n_to_keep <- ceiling((1 - fnr) * nrow(.x))
        .x[sample(nrow(.x), n_to_keep), , drop = FALSE]
      }) %>%
      ungroup()
  })
  
  # Name each element by its FNR
  names(TRUTH_filtered_list) <- as.character(unique_fnr_values)
  
  # Fill FNR vector with TRUTH for easier handling
  FNRs <- FNRs %>%
    mutate(
      FNR = map2(FNR, offset_naive_coi, ~ {
        len <- length(.x)
        if (len < .y) {
          c(.x, rep("TRUTH", .y - len))
        } else {
          .x
        }
      })
    )
  
  # ADD TRUTH DF TO THE FILTEREDLIST as element 1 FOR EASIER MANAGING
  TRUTH_filtered_list <- c(list(as_tibble(TRUTH)), TRUTH_filtered_list)
  names(TRUTH_filtered_list)[1] <- "TRUTH"
  
  # Initialize storage objects
  final_sample_df <- tibble()
  TRUTH_samples_by_row <- list()
  
  cat("Processing FNR schemes...\n")
  
  # Iterate over each row in FNRs
  for (numba in seq(1, nrow(FNRs), 1)) {
    
    # Get iteration data
    iteration <- FNRs[numba, ]
    
    cat(paste0("PROCESSING: offset_naive_coi = ", iteration$offset_naive_coi, 
               "; FNR = ", paste(unlist(iteration$FNR), collapse = "_"), "\n"))
    
    # Subset FNR_pairs_METADATA based on offset_naive_coi
    coi_subset <- FNR_pairs_METADATA %>%
      group_by(PairsID) %>%
      filter(any(n_strains == iteration$offset_naive_coi)) %>%
      ungroup() %>%
      rename(offset_naive_coi = n_strains) %>%
      full_join(iteration, by = "offset_naive_coi")
    
    # Randomly sample 10% of the pairs for FNR scheme
    n_size <- round(length(unique(coi_subset$PairsID))*0.1, 1)
    
    set.seed(pairs_seed)
    sample_pairsid <- sample(unique(coi_subset$PairsID), n_size)
    
    coi_subset <- coi_subset[coi_subset$PairsID %in% sample_pairsid,]
    
    # Add unique ID to PairsID
    coi_subset$PairsID <- paste0(coi_subset$PairsID, "_FNR", numba)
    
    # Function to process one row
    process_row_optimized <- function(pid, tp, strains, fnrs) {
      # Handle NULL FNR case
      if (is.null(fnrs)) {
        fnrs <- rep("TRUTH", length(strains))
      }
      
      # Create strain-fnr pairs for vectorized processing
      strain_fnr_df <- data.frame(
        strain = strains,
        fnr = fnrs,
        stringsAsFactors = FALSE
      )
      
      # Process all strain-fnr combinations at once
      sampled_data <- strain_fnr_df %>%
        pmap_dfr(function(strain, fnr) {
          TRUTH_filtered_list[[fnr]] %>%
            filter(sampleID == strain)
        })
      
      # Create FNR suffix
      fnr_suffix <- paste(fnrs, collapse = "_")
      
      # Return results
      list(
        list_name = paste0(pid, "_", tp, "__", fnr_suffix),
        list_data = sampled_data,
        df_data = tibble(
          PairsID = paste0(pid, "__", fnr_suffix),
          time_point = tp,
          locus = sampled_data$locus,
          pseudo_cigar = sampled_data$pseudo_cigar,
          allele = sampled_data$allele
        )
      )
    }
    
    # Apply processing function to all rows in this iteration
    all_results <- pmap(coi_subset, function(PairsID, time_point, sample, offset_naive_coi, FNR) {
      process_row_optimized(PairsID, time_point, sample, FNR)
    })
    
    # Append TRUTH sample list
    TRUTH_samples_by_row <- c(
      TRUTH_samples_by_row,
      set_names(
        map(all_results, "list_data"),
        map_chr(all_results, "list_name")
      )
    )
    
    # Append to final_sample_df
    final_sample_df <- bind_rows(
      final_sample_df,
      map_dfr(all_results, "df_data")
    )
  }
  
  # Deduplicate and enrich final_sample_df
  final_sample_df <- final_sample_df %>%
    distinct() %>%
    mutate(NIDA = paste0(PairsID, "__", time_point)) %>%
    rename(FNR_scheme = PairsID) %>%
    mutate(PairsID = sub("__.*", "", FNR_scheme)) %>%
    relocate(PairsID, .before = FNR_scheme)
  
  # Create final_sample_df_METADATA
  final_sample_df_METADATA <- final_sample_df %>%
    select(PairsID, time_point, NIDA) %>%
    distinct() %>%
    pivot_wider(
      names_from = time_point,
      values_from = NIDA,
      names_prefix = "NIDA"
    ) %>%
    rename(
      NIDA1 = NIDAD0,
      NIDA2 = NIDADx
    )
  
  cat("Calculating features for FNR data...\n")
  
  #########3 CALCULATE FEATURES OF FNR DATA
  PAIRS_GENOMIC_FNR <- as.data.table(final_sample_df)
  PAIRS_METADATA_FNR <- as.data.table(final_sample_df_METADATA)
  
  PAIRS_GENOMIC_FNR <- left_join(PAIRS_GENOMIC_FNR, PAIRS_METADATA_FNR, by = "PairsID") 
  
  # Use the same feature calculation function from script 7
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
  
  # Extract unique pairs once
  unique_pairs <- unique(PAIRS_METADATA_FNR$PairsID)
  
  # Create a named vector for mapping PairsID to D0 and Dx for faster access
  pairs_to_d0 <- PAIRS_METADATA_FNR$NIDA1[match(unique_pairs, PAIRS_METADATA_FNR$PairsID)]
  pairs_to_dx <- PAIRS_METADATA_FNR$NIDA2[match(unique_pairs, PAIRS_METADATA_FNR$PairsID)]
  
  # Preprocess merged_dfs_filtered into a list for fast access
  merged_dfs_filtered <- split(PAIRS_GENOMIC_FNR, PAIRS_GENOMIC_FNR$NIDA)
  merged_dfs_filtered <- lapply(merged_dfs_filtered, function(df) {
    df %>% select(-PairsID, -time_point) %>% distinct()
  })
  
  gc()
  
  # Initialize the progress bar
  pb <- progress_bar$new(
    format = "Processing [:bar] :percent :elapsedfull",
    total = ceiling(length(unique_pairs) / chunk_size),
    width = 100
  )
  
  # Temporary file for writing output
  output_file <- paste0(site, "_TRUTH_FNR_TRAINING_DATA.csv")
  
  # Function to process each pair
  process_pair <- function(pair_id) {
    sample1 <- merged_dfs_filtered[[pairs_to_d0[pair_id]]] 
    sample2 <- merged_dfs_filtered[[pairs_to_dx[pair_id]]]
    
    # Calculate features
    feats <- as.data.frame(t(calculate_features_optimized(sample1, sample2)))
    feats$PairsID <- unique_pairs[pair_id]
    
    return(feats)
  }
  
  # Function to write results in chunks to avoid memory bloat
  write_chunk_to_csv <- function(data, file_path, append = FALSE) {
    if (append) {
      fwrite(data, file_path, append = TRUE)
    } else {
      fwrite(data, file_path)
    }
  }
  
  gc()
  
  # Process unique_pairs in chunks
  for (i in seq(1, length(unique_pairs), by = chunk_size)) {
    # Define chunk range
    chunk_pair_ids <- i:min(i + chunk_size - 1, length(unique_pairs))
    
    cat(paste0("Processing PairsID ", unique_pairs[chunk_pair_ids[1]], 
               " to ", unique_pairs[chunk_pair_ids[length(chunk_pair_ids)]], "...\n"))
    
    # Parallel processing of the current chunk
    results <- mclapply(chunk_pair_ids, process_pair, mc.cores = num_cores)
    
    # Filter out NULL results
    metrics_list_chunk <- Filter(Negate(is.null), results)
    
    # Check if there are valid results to write
    if (length(metrics_list_chunk) > 0) {
      # Combine results into a data table
      delta_features_df_chunk <- rbindlist(metrics_list_chunk, fill = TRUE)
      
      # Incrementally write results to CSV
      write_chunk_to_csv(delta_features_df_chunk, output_file, append = (i > 1))
    }
    
    pb$tick()
  }
  
  delta_metrics_df_final_FNR <- read.csv(output_file)
  
  delta_metrics_df_final_FNR <- delta_metrics_df_final_FNR %>%
    select(PairsID, everything())
  
  write.csv(delta_metrics_df_final_FNR, output_file, row.names = FALSE)
  
  cat("Merging FNR data with original training data...\n")
  
  ### APPEND FNR DATA TO ORIGINAL DATA
  FEATURES <- left_join(delta_metrics_df_final_FNR, PAIRS_METADATA_FNR, by = "PairsID")
  
  FEATURES <- FEATURES %>%
    mutate(PairsID_clean = sub("_.*", "", PairsID))
  
  FEATURES$PairsID_clean <- as.numeric(FEATURES$PairsID_clean)
  
  ### MERGE WITH CLEAN DATA'S METADATA
  merged_df <- FEATURES %>%
    left_join(TRUTH_CLEAN[c("D0_sample", "Dx_sample", "D0_nstrains", "Dx_nstrains", "PairsID", "coi_change", "labels")], 
              by = c("PairsID_clean" = "PairsID")) %>%
    select(-PairsID_clean)
  
  merged_df$PairsID <- as.character(merged_df$PairsID)
  
  # Output final features with TRUTH + TRUTH_FNR data
  write.csv(merged_df, paste0(site, "_FNR_TRAINING_DATA.csv"), row.names = FALSE)
  
  ### CONCAT WITH CLEAN DATA
  # Find common columns between the two dataframes
  common_cols <- intersect(names(merged_df), names(TRUTH_CLEAN))
  
  TRUTH_CLEAN$PairsID <- as.character(TRUTH_CLEAN$PairsID)
  
  # Bind the data using only the common columns
  CLEAN_and_FNR_final <- bind_rows(
    merged_df %>% select(all_of(common_cols)),
    TRUTH_CLEAN %>% select(all_of(common_cols))
  )
  
  # Delete columns with NA
  CLEAN_and_FNR_final <- CLEAN_and_FNR_final %>% 
    filter(if_all(everything(), ~ !is.na(.)))
  
  CLEAN_and_FNR_final <- unique(CLEAN_and_FNR_final)
  
  # OVERWRITE ORIGINAL TRAINING DATA WITH THE APPENDED ORIGINAL + FNR DATA
  write.csv(CLEAN_and_FNR_final, paste0(site, "_TRAINING_DATA.csv"), row.names = FALSE)
  
  cat("=== TARGET AUGMENTATION COMPLETE ===\n")
  cat("Training data augmented with FNR samples\n\n")
  
} else {
  cat("All strain proportions from D0 failure samples are > 0.03. No target augmentation needed.\n")
  cat("=== TARGET AUGMENTATION SKIPPED ===\n\n")
  
  # Just copy the training data without augmentation
  training_files <- list.files(".", pattern = ".*_TRAINING_DATA\\.csv", full.names = TRUE)
  if (length(training_files) > 0) {
    training_data <- read.csv(training_files[1])
    write.csv(training_data, paste0(site, "_TRAINING_DATA.csv"), row.names = FALSE)
  }
}
