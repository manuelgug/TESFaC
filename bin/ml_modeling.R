#!/usr/bin/env Rscript

# ========================================
# 9. ML MODELING (Nextflow version)
# ========================================

# Load required libraries
suppressPackageStartupMessages({
  library(caret)    
  library(dplyr)    
  library(tidyr)    
  library(ggplot2)  
  library(broom)
  library(glmnet)
})

# Get parameters from environment
site <- Sys.getenv("NXF_SITE", "")
train_test_split <- as.numeric(Sys.getenv("NXF_TRAIN_TEST_SPLIT", "0.7"))
cv_folds <- as.numeric(Sys.getenv("NXF_CV_FOLDS", "10"))

cat("=== STARTING ML MODELING ===\n")

### 1) IMPORT TRAINING AND REAL DATA ----------
cat("Loading training and real data...\n")

# Find input files
training_files <- list.files(".", pattern = ".*_TRAINING_DATA\\.csv", full.names = TRUE)
real_files <- list.files(".", pattern = ".*_REAL_DATA\\.csv", full.names = TRUE)

if (length(training_files) == 0) stop("No training data file found")
if (length(real_files) == 0) stop("No real data file found")

TRAINING_DATA <- read.csv(training_files[1], row.names = 1) 
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains)

LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

REAL_DATA <- read.csv(real_files[1], 
                      stringsAsFactors = FALSE, 
                      colClasses = c(NIDA1 = "character", NIDA2= "character")) 

# Filter features (exclude non-predictive columns)
excluded_cols <- c("PairsID", "NIDA1", "NIDA2", "pair_type", "IBD_estimate",
                   "offset_naive_coi_D0", "offset_naive_coi_Dx", 
                   "replacement_pattern_score", "locus_discordance_rate")

features_to_use_ml <- colnames(REAL_DATA)[!colnames(REAL_DATA) %in% excluded_cols]

cat("Features selected for modeling:\n")
print(features_to_use_ml)

# Correlation plot
png(paste0(site, "_feature_correlation.png"), 
    width = 2000, height = 1600, res = 300)
corrplot::corrplot(cor(TRAINING_DATA %>% select(features_to_use_ml), use = "complete.obs"), "pie")
dev.off()

### 2) SPLIT DATA ------------
cat("Splitting data into train/test sets...\n")

set.seed(420)

# Stratified sampling by `pair_type`
train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = train_test_split, list = FALSE)
train_data <- TRAINING_DATA[train_indices, ]
test_data <- TRAINING_DATA[-train_indices, ]

## Output TEST for later use in false negative tests
test_data$PairsID <- rownames(test_data)
test_data <- test_data %>% select(PairsID, everything())
saveRDS(test_data, paste0(site, "_test_data.RDS"))

# Extract Metadata and Labels
TRAIN_META <- train_data %>% select(-all_of(features_to_use_ml))
TRAIN <- train_data %>% select(all_of(features_to_use_ml))
TRAIN_labels <- LABELS[train_indices, ]

TEST_META <- test_data %>% select(-all_of(features_to_use_ml))
TEST <- test_data %>% select(all_of(features_to_use_ml))
TEST_labels <- LABELS[-train_indices, ]

# Verify Stratification
cat("Training set stratification:\n")
print(prop.table(table(TRAIN_META$eCOI_pairs)))
print(prop.table(table(TRAIN_labels)))
print(table(TRAIN_labels))

cat("Test set stratification:\n")
print(prop.table(table(TEST_META$eCOI_pairs)))
print(prop.table(table(TEST_labels)))
print(table(TEST_labels))

##### 4) TRAIN MODEL USING GLMNET --------
cat("Training GLMNET model...\n")

# Create training and test data frames
df_train_IBD <- data.frame(TRAIN[features_to_use_ml], label = as.factor(TRAIN_labels))
df_test_IBD <- data.frame(TEST[features_to_use_ml], label = as.factor(TEST_labels))

# Set up cross-validation
ctrl <- trainControl(
  method = "cv",
  number = cv_folds,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")

# Train GLMNET model
fit_IBD <- train(
  label ~ .,
  data = df_train_IBD,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 10
)

cat("Model training complete:\n")
print(fit_IBD)

# Feature importance (GLMNET)
cat("Calculating feature importance...\n")

best_lambda <- fit_IBD$bestTune$lambda
coefs <- coef(fit_IBD$finalModel, s = best_lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
colnames(coefs_df) <- "Estimate"
coefs_df$Variable <- rownames(coefs_df)
coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate) + 1e-8)  # avoid log(0)
coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")

importance <- ggplot(coefs_df, aes(x = reorder(Variable, Estimate), y = Estimate, fill = Color)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  theme_minimal() +
  labs(x = "Features", y = "Coefficient Value", fill = "Coefficient Direction") +
  theme(legend.position = "bottom")

ggsave(paste0("feat_importance_", site, ".png"), 
       importance, bg = "white", dpi = 300, height = 5, width = 8)

# Predict on test set
cat("Evaluating model performance...\n")

decision_thresholds <- seq(0, 1, by = 0.05)
preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]

results <- data.frame(eCOI_pairs = character(),
                      decision_threshold = numeric(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

for (thresh in decision_thresholds) {
  preds <- ifelse(preds_prob >= thresh, "R", "NI")
  
  for (strain_comb in unique(TEST_META$eCOI_pairs)) {
    subset_indices <- TEST_META$eCOI_pairs == strain_comb
    subset_TEST_labels <- TEST_labels[subset_indices]
    
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    subset_preds <- preds[subset_indices]
    cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
                                         decision_threshold = thresh,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
}

# Best thresholds (Youden's J)
cat("Optimizing decision thresholds...\n")

best_decision_thresholds <- results %>%
  mutate(youden_j = sensitivity + specificity - 1) %>%
  group_by(eCOI_pairs) %>%
  slice_max(youden_j) %>%
  slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%
  ungroup()

best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity),
               names_to = "Metric",
               values_to = "Value")

metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Pair Type", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")

# Save results
saveRDS(fit_IBD, paste0("GLMNET_model_", site, ".RDS"))
write.csv(best_decision_thresholds, 
          paste0("training_Results_GLMNET_", site, ".csv"), 
          row.names = FALSE)
ggsave(paste0(site, "_model_results_GLMNET.png"), 
       metrics, bg = "white", dpi = 300, height = 5, width = 7)

# Reshape data into long format for plotting
results_long <- results %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Plot both Sensitivity and Specificity in the same graph
sens_spec_plot <- ggplot(results_long, aes(x = decision_threshold, y = Value, linetype = Metric)) +
  geom_line() +
  geom_vline(data = best_decision_thresholds, aes(xintercept = decision_threshold), 
             color = "red", linetype = "solid") +
  facet_wrap(~eCOI_pairs) +
  labs(title = "", x = "Decision Threshold", y = "Value") +
  theme_minimal()

##### 6) TEST MODEL ON REAL DATA --------
cat("Applying model to real data...\n")

# Add threshold data
best_decision_thresholds <- best_decision_thresholds %>% rename(pair_type = eCOI_pairs)
REAL_DATA <- left_join(REAL_DATA, 
                       best_decision_thresholds[c("pair_type", "decision_threshold")], 
                       by = c("pair_type"))

REAL_DATA$prediction_prob <- NA
REAL_DATA$predictions <- NA

# Loop over each row in REAL_DATA to apply the corresponding decision_threshold
for (i in 1:nrow(REAL_DATA)) {
  
  decision_threshold <- REAL_DATA$decision_threshold[i]
  
  feats <- REAL_DATA %>% select(all_of(features_to_use_ml))
  
  newdata <- feats[i, , drop = FALSE]
  
  prediction_prob <- predict(fit_IBD, newdata = newdata, type = "prob")[, "R"]
  
  # Classify using the current (best) decision_threshold (instead of 0.5)
  prediction_class <- ifelse(prediction_prob >= decision_threshold, "R", "NI")
  
  REAL_DATA$prediction_prob[i] <- prediction_prob
  REAL_DATA$predictions[i] <- prediction_class
}

REAL_DATA <- REAL_DATA %>% 
  select(PairsID, NIDA1, NIDA2, pair_type, all_of(features_to_use_ml), 
         decision_threshold, prediction_prob, predictions) %>% 
  arrange(PairsID)

## OUTPUT RESULTS 
write.csv(REAL_DATA, 
          paste0(site, "_REAL_DATA_PREDICTIONS.csv"), 
          row.names = FALSE)

ggsave(paste0(site, "_REAL_DATA_THRESHOLDS_PLOT.png"), 
       sens_spec_plot, bg = "white", dpi = 300, height = 9, width = 12)

# Final summary
cat("=== ML MODELING COMPLETE ===\n")
cat("Model performance summary:\n")
print(best_decision_thresholds %>% select(pair_type, sensitivity, specificity, decision_threshold))

cat("\nReal data predictions summary:\n")
print(table(REAL_DATA$predictions))

cat("All results saved to current directory\n\n")
