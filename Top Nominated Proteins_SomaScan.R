#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library (readxl)
library(writexl)
library(gplots)
library(factoextra)

#Heatmap and PCA plot of the top nominated 61 proteins from pairwise and multigroup comparisons to classify stroke subtypes (SomaScan Proteomics)

#Proteomics data file upload
excel_sheets("Top_proteins_Combined.xlsx")
Top_data= excel_sheets("Top_proteins_Combined.xlsx") %>% map(~read_xlsx("Top_proteins_Combined.xlsx",.))
Top_data

# Heatmap

#Clustering and Dendograms

#Convert to matrix data frame
dat_matrix = as.matrix.data.frame(Top_data[[2]][,2:101]) 
#Name the rows with protein ids
row.names(dat_matrix) <- Top_data[[2]]$Protein_ID
#Transpose and scale the data to a mean of zero and sd of one
dat_scaled <- scale(t(dat_matrix)) %>% t()

#Transpose the matrix to calculate distance between experiments, row-wise
d1 <- dat_scaled %>% t() %>%
  dist(.,method = "euclidean", diag = FALSE, upper = FALSE)
#Calculate the distance between proteins row-wise 
d2 <- dat_scaled %>%
  dist(.,method = "euclidean", diag = FALSE, upper = FALSE)

#Show the values for d1
round(d1,2)

#Clustering distance between experiments using Ward linkage
c1 <- hclust(d1, method = "ward.D2", members = NULL)
#Clustering distance between proteins using Ward linkage
c2 <- hclust(d2, method = "ward.D2", members = NULL)

#Check clustering by plotting dendrograms
X11()
par(mfrow=c(2,1),cex=0.5) # Make 2 rows, 1 col plot frame and shrink labels
plot(c1); plot(c2) # Plot both cluster dendrograms

#Creating a Heatmap

#Set colours for heatmap, 25 increments
my_palette = colorRampPalette(c("blue","white","red"))(n = 100)

#Plot heatmap with heatmap.2
X11()
par(cex.main=0.75) # Shrink title fonts on plot
dat_scaled %>% 
  #Plot heatmap
  gplots::heatmap.2(.,                          # Tidy, normalised data
                    Colv=FALSE,                 # Experiments clusters in cols
                    Rowv=as.dendrogram(c2),     # Protein clusters in rows
                    revC=TRUE,                  # Flip plot to match pheatmap
                    density.info="histogram",   # Plot histogram of data and colour key
                    trace="none",               # Turn of trace lines from heat map
                    col = my_palette,           # Use my colour scheme
                    cexRow=0.6, cexCol=0.75,
                    scale = "none",
                    key=TRUE)             # Amend row and column label fonts

# Split data into different clusters based on protein clusters
clustered_data <- split(dat_scaled, cutree(as.dendrogram(c2), k = n_clusters))

#Create a heatmap using pheatmap package
dev.off() #Use this if pheatmap function is not working

#Plot pheatmap
dat_scaled %>%
  pheatmap(.,
           fontsize = 7,
           cutree_rows = 2, clutree_cols= 4, cluster_cols = FALSE, col= redblue(20), scale = "none", cexRow=0.6, cexCol=0.75,
           border_color = "black") # Create breaks in heatmap

##################################

#PCA plot using top proteins
R.PCA<- prcomp(Top_data[[1]][,3:63], scale=TRUE)
R.PCA
X11()
pca_plot=fviz_pca_ind(R.PCA, col.ind=Top_data[[1]]$Outcome, title= "PCA plot of combined top proteins", addEllipses = FALSE,
                      label= "none", pointsize= 4)

# Define colors for the groups
colors = c("AIS" = "red", "ICH" = "blue", "TIA" = "purple", "MIM"= "green")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

##############################################################################################################################################

## To calculate the AUC and 95% CI for 61 nominating protein hits to classify AIS, ICH, TIA, and MIM

library(readxl)
library(writexl)
library(pROC)

# Step 1: Load the data
proteomics_file <- "Protein_Classifiers.xlsx"
df <- read_excel(proteomics_file)

# Step 2: Define outcome and protein columns
outcome_vars <- c("AIS", "ICH", "TIA", "MIM")
protein_vars <- colnames(df)[6:ncol(df)]  # Columns 6 onward are proteins

# Step 3: Compute AUC + CI and format
results_list <- list()

for (outcome in outcome_vars) {
  formatted_results <- data.frame(
    Protein_ID = protein_vars,
    stringsAsFactors = FALSE
  )
  
  formatted_auc_ci <- character(length(protein_vars))
  
  for (i in seq_along(protein_vars)) {
    protein <- protein_vars[i]
    roc_obj <- try(roc(df[[outcome]], df[[protein]], quiet = TRUE), silent = TRUE)
    
    if (!inherits(roc_obj, "try-error")) {
      auc_val <- round(as.numeric(auc(roc_obj)), 2)
      ci_vals <- round(ci.auc(roc_obj), 2)
      formatted_auc_ci[i] <- sprintf("%.2f (%.2f–%.2f)", auc_val, ci_vals[1], ci_vals[3])
    } else {
      formatted_auc_ci[i] <- NA
    }
  }
  
  formatted_results[[outcome]] <- formatted_auc_ci
  results_list[[outcome]] <- formatted_results
}

# Step 4: Merge results by Protein_ID
final_results <- results_list[[1]][, "Protein_ID", drop = FALSE]

for (outcome in outcome_vars) {
  final_results <- merge(final_results, results_list[[outcome]], by = "Protein_ID", all = TRUE)
}

# Step 5: Export to Excel
write_xlsx(final_results, "Protein_Classifiers_AUC.xlsx")

##################################################################################################################################################
## Creating protein prediction models from 61 nominating hits using regularized LASSO regression and plotting ROC curves and box plots

#Load packages
library(readxl)
library(hdm)
library(caret)
library(boot)
library(pROC)
library(ggplot2)
library(tidyverse)

# Load data
data <- read_excel("Protein_classifiers.xlsx")
View(data)
y <- as.numeric(data$MIM)
X <- as.matrix(data[, -(1:5)])
predictors <- names(data)[-(1:5)]
predictors

# Fit LASSO logistic regression model
model <- rlassologit(X, y)

# Selected proteins
selected_vars <- predictors[which(model$index)]
cat("Selected proteins in the model:\n")
print(selected_vars)

# Predict probabilities on full data
preds_prob <- as.numeric(predict(model, newdata = X, type = "response"))

# Threshold for classification
threshold <- 0.20
preds_class <- ifelse(preds_prob >= threshold, 1, 0)

# Make sure predicted and true classes are factors with levels "1" and "0" in that order
preds_factor <- factor(preds_class, levels = c(1, 0))
true_factor <- factor(y, levels = c(1, 0))

# Confusion matrix with desired ordering
cm <- confusionMatrix(preds_factor, true_factor, positive = "1")
cm

cat("\nConfusion matrix counts (predicted 1 then 0):\n")
print(cm$table)

# Extract point estimates safely
safe_metric <- function(x) ifelse(is.na(x), NA_real_, x)
metrics_estimate <- c(
  Sensitivity = safe_metric(cm$byClass["Sensitivity"]),
  Specificity = safe_metric(cm$byClass["Specificity"]),
  PPV = safe_metric(cm$byClass["Pos Pred Value"]),
  NPV = safe_metric(cm$byClass["Neg Pred Value"])
)

cat("\nPoint estimates at specified threshold:\n")
print(metrics_estimate)

# Prepare data for bootstrap
boot_data <- data.frame(
  pred_class = factor(preds_class, levels = c(0,1)),
  true_class = factor(y, levels = c(0,1))
)

# Bootstrap function for metrics
boot_metrics <- function(data, indices) {
  d <- data[indices, ]
  cm_boot <- confusionMatrix(d$pred_class, d$true_class, positive = "1")
  c(
    Sensitivity = as.numeric(cm_boot$byClass["Sensitivity"]),
    Specificity = as.numeric(cm_boot$byClass["Specificity"]),
    PPV = as.numeric(cm_boot$byClass["Pos Pred Value"]),
    NPV = as.numeric(cm_boot$byClass["Neg Pred Value"])
  )
}

set.seed(123)
boot_results <- boot(data = boot_data, statistic = boot_metrics, R = 1000)

# Calculate 95% CIs (percentile)
ci_metrics <- t(apply(boot_results$t, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))
colnames(ci_metrics) <- c("2.5%", "97.5%")
rownames(ci_metrics) <- c("Sensitivity", "Specificity", "PPV", "NPV")

cat("\nMetrics with 95% bootstrap confidence intervals:\n")
for (metric in rownames(ci_metrics)) {
  cat(sprintf(
    "%s: %.4f (95%% CI: %.4f - %.4f)\n",
    metric,
    metrics_estimate[metric],
    ci_metrics[metric, "2.5%"],
    ci_metrics[metric, "97.5%"]
  ))
}

# ROC curve and AUC with 95% CI (DeLong)
roc_obj <- roc(y, preds_prob)
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj, method = "delong")

cat(sprintf("\nAUC: %.4f (95%% CI: %.4f - %.4f)\n", auc_val, auc_ci[1], auc_ci[3]))

# Plot ROC
plot(roc_obj, col = "blue", main = "ROC Curve for LASSO Logistic Regression",
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)")
legend("bottomright", legend = c(
  sprintf("AUC = %.2f (95%% CI: %.2f - %.2f)", auc_val, auc_ci[1], auc_ci[3])), bty = "n", cex = 0.75)

## Box plot of selected protein classifiers

# Subset selected proteins
selected_data <- data[, selected_vars]

# Add Subtype column based on mutually exclusive binary columns
data$Subtype <- case_when(
  data$AIS == 1 ~ "AIS",
  data$ICH == 1 ~ "ICH",
  data$TIA == 1 ~ "TIA",
  data$MIM == 1 ~ "MIM",
  TRUE ~ NA_character_
)

# Ensure factor levels for Subtype
data$Subtype <- factor(data$Subtype, levels = c("AIS", "ICH", "TIA", "MIM"))

# Convert to long format
long_data <- selected_data %>%
  mutate(SubjectID = row_number(), Subtype = data$Subtype) %>%
  pivot_longer(cols = all_of(selected_vars), names_to = "Protein", values_to = "Log2_Concentration")

# Order proteins by descending median log2 concentration
protein_order <- long_data %>%
  group_by(Protein) %>%
  summarize(median_conc = median(Log2_Concentration, na.rm = TRUE)) %>%
  arrange(desc(median_conc)) %>%
  pull(Protein)

# Apply protein order
long_data$Protein <- factor(long_data$Protein, levels = protein_order)

# Create boxplot
ggplot(long_data, aes(x = Protein, y = Log2_Concentration, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  labs(
    title = "Box plot of Protein Classifiers for stroke subtypes",
    x = "Protein Classifiers",
    y = "Protein Concentration (Log2)"
  ) +
  scale_fill_manual(
    values = c(
      "AIS" = "red",
      "ICH" = "blue",
      "TIA" = "purple",
      "MIM" = "green"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########################################################################################################################################

## Repeated Nested Cross-validation of protein panels

library(readxl)
library(hdm)
library(pROC)
library(caret)
library(boot)
library(ggplot2)

# Load data
data <- read_excel("Protein_classifiers.xlsx")

# Outcome and predictors
y <- as.numeric(data$TIA)
X <- as.matrix(data[, -(1:5)])
predictors <- colnames(X)

threshold <- 0.10  # classification threshold

# Number of outer folds and repetitions
outer_folds <- 10
repeats <- 5

# Store all predictions and labels across repeats and folds
all_preds <- numeric(0)
all_labels <- numeric(0)

# Store selected variables from each outer fold and repetition
selected_vars_list <- list()

# Compute metrics function
compute_metrics <- function(cm) {
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  ppv <- cm$byClass["Pos Pred Value"]
  npv <- cm$byClass["Neg Pred Value"]
  return(c(sensitivity = sens, specificity = spec, PPV = ppv, NPV = npv))
}

# Bootstrap function to calculate metrics confidence intervals
boot_metrics <- function(data, indices) {
  d <- data[indices, ]
  cm <- confusionMatrix(as.factor(d$pred_class), as.factor(d$true_class), positive = "1")
  as.numeric(compute_metrics(cm))
}

set.seed(1234) # For reproducibility

for (r in 1:repeats) {
  cat(sprintf("\n--- Repeat %d of %d ---\n", r, repeats))
  
  # Create folds anew for each repetition
  folds <- createFolds(y, k = outer_folds, list = TRUE, returnTrain = FALSE)
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    X_train <- X[train_idx, , drop=FALSE]
    y_train <- y[train_idx]
    
    X_test <- X[test_idx, , drop=FALSE]
    y_test <- y[test_idx]
    
    # Fit LASSO logistic model on training data
    model <- rlassologit(X_train, y_train)
    
    # Store selected variables for this fold
    selected_vars_fold <- predictors[which(model$index)]
    selected_vars_list <- c(selected_vars_list, list(selected_vars_fold))
    
    # Predict probabilities on test set
    preds_prob <- predict(model, newdata = X_test, type = "response")
    
    # Accumulate predictions and true labels
    all_preds <- c(all_preds, preds_prob)
    all_labels <- c(all_labels, y_test)
  }
}

# Overall ROC curve and AUC on combined predictions
roc_obj <- roc(all_labels, all_preds)
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj, method = "delong")

# Classification at threshold on combined data
all_class <- ifelse(all_preds >= threshold, 1, 0)
cm_overall <- confusionMatrix(as.factor(all_class), as.factor(all_labels), positive = "1")
overall_metrics <- compute_metrics(cm_overall)

# Bootstrap 95% CI for classification metrics
test_data_boot <- data.frame(
  pred_class = factor(all_class, levels = c(0,1)),
  true_class = factor(all_labels, levels = c(0,1))
)
set.seed(789)
boot_obj <- boot(data = test_data_boot, statistic = boot_metrics, R = 2000)
ci <- t(apply(boot_obj$t, 2, quantile, probs = c(0.025, 0.975)))
colnames(ci) <- c("2.5%", "97.5%")
rownames(ci) <- c("sensitivity", "specificity", "PPV", "NPV")

# Train final model on full dataset
final_model <- rlassologit(X, y)
final_selected_vars <- predictors[which(final_model$index)]

# Output selected proteins summary
cat("\nFinal selected proteins on full dataset:\n")
print(final_selected_vars)

# Optional: summary of selected variables frequency across folds
all_selected_vars <- unlist(selected_vars_list)
freq_table <- sort(table(all_selected_vars), decreasing = TRUE)
cat("\nFrequency of protein selection across folds and repeats:\n")
print(freq_table)

# Output performance metrics with CIs
cat("\nOverall model performance (repeated nested 10-fold CV):\n")
metrics_table <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity", "PPV", "NPV"),
  Estimate = c(as.numeric(auc_val), overall_metrics),
  Lower_95_CI = c(auc_ci[1], ci[,1]),
  Upper_95_CI = c(auc_ci[3], ci[,2])
)
print(metrics_table)

# Plot ROC curve
plot(roc_obj, col = "blue", main = "Repeated Nested 10-fold CV ROC Curve (LASSO Logistic Regression)",
     xlab = "1 - Specificity (False Positive Rate)", ylab = "Sensitivity (True Positive Rate)")
legend("bottomright", legend = c(
  sprintf("AUC = %.2f (95%% CI: %.2f - %.2f)", auc_val, auc_ci[1], auc_ci[3])
), bty = "n", cex = 0.75)


## Visualize stability of protein classifiers across folds and repeats selected using LASSO regression

# Convert the frequency table into a data frame
freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("Protein", "Frequency")

# Sort by descending frequency for plotting
freq_df <- freq_df[order(freq_df$Frequency, decreasing = TRUE), ]

# Optional: limit to top N proteins (e.g., top 10)
top_n <- 8
freq_df_plot <- head(freq_df, top_n)

# Plot
ggplot(freq_df_plot, aes(x = reorder(Protein, Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip for readability
  labs(
    title = "Stability of LASSO-selected Top Protein Classifiers",
    x = "Protein",
    y = "Selection Frequency Across Repeats and Folds"
  ) +
  scale_y_continuous(breaks = seq(0, max(freq_df_plot$Frequency), by = 5)) +  # Set axis ticks every 5
  theme_minimal(base_size = 13)