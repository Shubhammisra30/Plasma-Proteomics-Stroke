#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library(readxl)
library(writexl)
library(gplots)
library(ggrepel)
library(factoextra)
library(variancePartition)
library(ggplot2)
library(Matrix)
library(sva)
library(limma)
library(reshape2)

#Proteomics data file upload
excel_sheets("PEP_data.xlsx")
PEP_data= excel_sheets("PEP_data.xlsx") %>% map(~read_xlsx("PEP_data.xlsx",.))
PEP_data

###########################################
## Log Transformation ##
###########################################

#Select columns and log transform the data
dat_log <- PEP_data[[5]] %>%  
  select(-c(Protein_ID,Outcome, Batch)) %>% 
  log2()

# Read metadata columns
dat_log <- bind_cols(
  PEP_data[[5]] %>% select(Protein_ID, Outcome, Batch),
  dat_log
)

#Save log transformed data in excel file
write_xlsx(dat_log, "PEP_Log2_Transformed_.xlsx")

############################################################################
# Distribution of Unnormalized, log transformed, Batch uncorrected DIA data
############################################################################

# Separate metadata and expression
meta_data <- PEP_data[[7]] %>% select(Protein_ID, Outcome, Batch)
expr_data <- PEP_data[[7]] %>% select(-Protein_ID, -Outcome, -Batch)

# Transpose expression data: proteins in rows, samples in columns
expr_data_t <- as.data.frame(t(expr_data))

summary(as.matrix(expr_data_t))
median(as.matrix(expr_data_t))

# Assign sample names (optional, if Protein_IDs represent samples)
colnames(expr_data_t) <- meta_data$Protein_ID  

# Create colors based on Outcome
colors <- c("red", "green", "blue", "purple", "orange", "yellow", "pink", "cyan")
conditionNames <- unique(meta_data$Outcome)
colorVector <- setNames(colors[1:length(conditionNames)], conditionNames)
statusCol <- colorVector[meta_data$Outcome]

# Set margin sizes
par(mar = c(5, 4, 4, 2) + 0.1)

# Boxplot across all samples, colored by Outcome
med_val <- median(as.matrix(expr_data_t), na.rm = TRUE)
boxplot(expr_data_t,
        xlab = "",
        ylab = "Expression levels",
        las = 2,
        col = statusCol,
        main = "Box plot of all samples by Outcome")

# Add horizontal line at median expression
abline(h = med_val, col = "blue", lwd = 2)

##########################################################
## PCA Plot to check batch effect in unnormalized data##
##########################################################

R.PCA<- prcomp(PEP_data[[7]][,4:7156], scale=TRUE)
R.PCA
pca_plot= fviz_pca_ind(R.PCA, col.ind=PEP_data[[7]]$Batch, title= "PCA Plot showing Batch effect", addEllipses = FALSE, 
                       label= "none", pointsize= 3)
# Define colors for the groups
colors = c("Batch_1"= "red", "Batch_2" = "blue", "Batch_3" = "green")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

###########################################
## Total Area Sum Normalization
###########################################

# Separate metadata and expression
meta_data <- PEP_data[[7]] %>% select(Protein_ID, Outcome, Batch)
expr_data <- PEP_data[[7]] %>% select(-Protein_ID, -Outcome, -Batch)

# Transpose expression data: proteins in rows, samples in columns
expr_data_t <- as.data.frame(t(expr_data))
View(expr_data_t)

# Convert to matrix
expr_mat <- as.matrix(expr_data_t)
View(expr_mat)

# Calculate column sums (total intensity per sample)
total_intensity <- colSums(expr_mat, na.rm = TRUE)

# Scale factors: mean total / each sample’s total
scale_factors <- mean(total_intensity) / total_intensity

# Apply normalization (multiply each sample by its scale factor)
expr_data_norm <- sweep(expr_mat, 2, scale_factors, FUN = "*")

# Transpose back: Samples in rows and Proteins in columns
expr_data_norm_t <- as.data.frame(t(expr_data_norm))

# Back to data frame and reattach metadata
expr_data_norm_t <- as.data.frame(expr_data_norm_t)
expr_data_norm_t$Protein_ID <- meta_data$Protein_ID
expr_data_norm_t$Outcome <- meta_data$Outcome
expr_data_norm_t$Batch <- meta_data$Batch

# Reorder columns for clarity
expr_data_norm_t <- expr_data_norm_t %>%
  select(Protein_ID, Outcome, Batch, everything())
View(expr_data_norm_t)

# Write normalized data to Excel
write_xlsx(expr_data_norm_t, "TAS_Normalized_PEP.xlsx")

##############################################################################
# Distribution of TAS Normalized, log transformed, Batch uncorrected DIA data
##############################################################################

# Separate metadata and expression
meta_data <- PEP_data[[8]] %>% select(Protein_ID, Outcome, Batch)
expr_data <- PEP_data[[8]] %>% select(-Protein_ID, -Outcome, -Batch)

# Transpose expression data: proteins in rows, samples in columns
expr_data_t <- as.data.frame(t(expr_data))
View(expr_data_t)

summary(as.matrix(expr_data_t))
median(as.matrix(expr_data_t))

# Assign sample names (optional, if Protein_IDs represent samples)
colnames(expr_data_t) <- meta_data$Protein_ID  

# Create colors based on Outcome
colors <- c("red", "green", "blue", "purple", "orange", "yellow", "pink", "cyan")
conditionNames <- unique(meta_data$Outcome)
colorVector <- setNames(colors[1:length(conditionNames)], conditionNames)
statusCol <- colorVector[meta_data$Outcome]

# Set margin sizes
par(mar = c(5, 4, 4, 2) + 0.1)

# Boxplot across all samples, colored by Outcome
med_val <- median(as.matrix(expr_data_t), na.rm = TRUE)
boxplot(expr_data_t,
        xlab = "",
        ylab = "Expression levels",
        las = 2,
        col = statusCol,
        main = "Box plot of all samples by Outcome (Mean Norm)")

# Add horizontal line at median expression
abline(h = med_val, col = "blue", lwd = 2)

##############################################################
## PCA Plot to check batch effect in TAS normalized data##
##############################################################

R.PCA<- prcomp(PEP_data[[8]][,4:7156], scale=TRUE)
R.PCA
pca_plot= fviz_pca_ind(R.PCA, col.ind=PEP_data[[8]]$Batch, title= "PCA Plot showing Batch effect after Mean Norm", addEllipses = FALSE, 
                       label= "none", pointsize= 3)
# Define colors for the groups
colors = c("Batch_1"= "red", "Batch_2" = "blue", "Batch_3" = "green")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

############################################################
## Batch correction using ComBat function in sva package ##
###########################################################

# Extract the data matrix (assuming rows = samples, columns = features)
expr_data <- t(PEP_data[[8]][,4:7156])  # Transpose for ComBat: features x samples
View(expr_data)

# Define batch variable (Centre) — must be a factor
batch <- as.factor(PEP_data[[8]]$Batch)

# Optional: model matrix for biological covariates (e.g., preserve variation due to Outcome variable)
modcombat <- model.matrix(~as.factor(Outcome), data=PEP_data[[8]])

# Apply ComBat to remove batch effects
combat_edata <- ComBat(dat=expr_data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# Transpose back to samples x features for PCA
combat_edata_t <- t(combat_edata)

# Run PCA on ComBat-corrected data
R.PCA_combat <- prcomp(combat_edata_t, scale=TRUE)

# Plot PCA — coloring still by centre for visual confirmation
pca_plot_combat <- fviz_pca_ind(
  R.PCA_combat,
  col.ind = PEP_data[[8]]$Batch,
  title = "PCA Plot after TAS Norm and ComBat Batch Correction",
  addEllipses = FALSE,
  label = "none",
  pointsize = 3
)

# Define colors for centres
colors = c("Batch_1"= "red", "Batch_2" = "blue", "Batch_3" = "green")

# Add custom color scale
pca_plot_combat + scale_color_manual(values = colors)

# Save the Batch corrected data using Combat to an excel file
# Convert to data frame (optional but recommended)
combat_df <- as.data.frame(combat_edata_t)

# Add back metadata columns from PEP_data[[8]]
combat_df <- cbind(
  PEP_data[[8]] %>% dplyr::select(Protein_ID, Outcome, Batch),
  combat_df
)

# Save to Excel
write_xlsx(combat_df, "PEP_Combat_Batch_Corrected.xlsx")

################################################################################################################################################

# Differential expression analysis using Pairwise comparisons between Stroke subtypes in the External Validation Cohort

##############
# AIS vs ICH
##############

#Visualizing the Data
dat = PEP_data[[10]]
colnames(dat)
View(dat)

#Creating a Welch's T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply Welch's t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc <- dat_combine %>%
  mutate(
    mean_AIS_case = rowMeans(select(., starts_with("AIS")), na.rm = TRUE),
    mean_ICH_case = rowMeans(select(., starts_with("ICH")), na.rm = TRUE),
    log_fc = mean_AIS_case - mean_ICH_case,
    log_pval = -log10(p_val)
  )
View(dat_fc)

#Save final data with list of final data in excel file
write_xlsx(dat_fc, "Final_AIS_ICH_data.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
X11()
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, mean_AIS_case, AIS1410_1_1:AIS185_40_1, ICH3219_3_1:ICH245_34_1, mean_ICH_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in excel file
write_xlsx(final_data, "Final_diffproteins_AIS_ICH_data.xlsx")

#########################################################################################################################################
# AIS vs TIA
#########################################################################################################################################

# Differential analysis

#Visualizing the Data
dat = PEP_data[[11]]
colnames(dat)
View(dat)

#Creating a Welch's T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply Welch's t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc <- dat_combine %>%
  mutate(
    mean_AIS_case = rowMeans(select(., starts_with("AIS")), na.rm = TRUE),
    mean_TIA_case = rowMeans(select(., starts_with("TIA")), na.rm = TRUE),
    log_fc = mean_AIS_case - mean_TIA_case,
    log_pval = -log10(p_val)
  )
View(dat_fc)

#Save final data with list of final data in excel file
write_xlsx(dat_fc, "Final_AIS_TIA_data.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, mean_AIS_case, mean_TIA_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in excel file
write_xlsx(final_data, "Final_diffproteins_AIS_TIA_data.xlsx")

#########################################################################################################################################
# AIS vs MIM
#########################################################################################################################################

#Visualizing the Data
dat = PEP_data[[12]]
colnames(dat)
View(dat)

#Creating a Welch's T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply Welch's t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc <- dat_combine %>%
  mutate(
    mean_AIS_case = rowMeans(select(., starts_with("AIS")), na.rm = TRUE),
    mean_MIM_case = rowMeans(select(., starts_with("OTH")), na.rm = TRUE),
    log_fc = mean_AIS_case - mean_MIM_case,
    log_pval = -log10(p_val)
  )
View(dat_fc)

#Save final data with list of final data in excel file
write_xlsx(dat_fc, "Final_AIS_MIM_data.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, mean_AIS_case, mean_MIM_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in excel file
write_xlsx(final_data, "Final_diffproteins_AIS_MIM_data.xlsx")

#########################################################################################################################################
# ICH vs TIA
#########################################################################################################################################

#Visualizing the Data
dat = PEP_data[[13]]
colnames(dat)
View(dat)

#Creating a Welch's T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply Welch's t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc <- dat_combine %>%
  mutate(
    mean_ICH_case = rowMeans(select(., starts_with("ICH")), na.rm = TRUE),
    mean_TIA_case = rowMeans(select(., starts_with("TIA")), na.rm = TRUE),
    log_fc = mean_ICH_case - mean_TIA_case,
    log_pval = -log10(p_val)
  )
View(dat_fc)

#Save final data with list of final data in excel file
write_xlsx(dat_fc, "Final_ICH_TIA_data.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, mean_ICH_case, mean_TIA_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in excel file
write_xlsx(final_data, "Final_diffproteins_ICH_TIA_data.xlsx")

#########################################################################################################################################
# ICH vs MIM
#########################################################################################################################################

#Visualizing the Data
dat = PEP_data[[14]]
colnames(dat)
View(dat)

#Creating a Welch's T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply Welch's t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc <- dat_combine %>%
  mutate(
    mean_ICH_case = rowMeans(select(., starts_with("ICH")), na.rm = TRUE),
    mean_MIM_case = rowMeans(select(., starts_with("OTH")), na.rm = TRUE),
    log_fc = mean_ICH_case - mean_MIM_case,
    log_pval = -log10(p_val)
  )
View(dat_fc)

#Save final data with list of final data in excel file
write_xlsx(dat_fc, "Final_ICH_MIM_data.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, mean_ICH_case, mean_MIM_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in excel file
write_xlsx(final_data, "Final_diffproteins_ICH_MIM_data.xlsx")

#########################################################################################################################################
# TIA vs MIM
#########################################################################################################################################

#Visualizing the Data
dat = PEP_data[[15]]
colnames(dat)
View(dat)

#Creating a Welch's T-test function for multiple experiments
t_test <- function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

#Apply Welch's t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc <- dat_combine %>%
  mutate(
    mean_TIA_case = rowMeans(select(., starts_with("TIA")), na.rm = TRUE),
    mean_MIM_case = rowMeans(select(., starts_with("OTH")), na.rm = TRUE),
    log_fc = mean_TIA_case - mean_MIM_case,
    log_pval = -log10(p_val)
  )
View(dat_fc)

#Save final data with list of final data in excel file
write_xlsx(dat_fc, "Final_TIA_MIM_data.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
dat_fc %>% ggplot(aes(log_fc,log_pval)) + geom_point()

#Volcano plot

VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column of NAs
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP
#Add lines as before..
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Change point colors
VP3= VP2 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP3= VP2 + scale_color_manual(values=mycolors)
#Create a new column "proteinlabel" that will contain names of differentially expressed protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()
View(dat_fc)

#Plot the Volcano plot using all layers used so far
X11()
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data<-dat_fc %>%
  #Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  #Ungroup the data
  ungroup() %>% 
  #Select columns of interest
  select(Protein_ID, TIA762_1_2:TIA3168_13_2, OTH2618_8_2:OTH895_21_3,  mean_TIA_case, mean_MIM_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in excel file
write_xlsx(final_data, "Final_diffproteins_TIA_MIM_data.xlsx")

#########################################################################################################################################

# To calculate ANOVA, FDR corrected-ANOVA p-val, and pairwise p-value using Tukey HSD for stroke subtypes

# Load required packages
library(readxl)
library(dplyr)
library(stats)

# Step 1: Load your data
#Proteomics data file upload
excel_sheets("ANOVA_PEP_data.xlsx")
ANOVA_PEP_data= excel_sheets("ANOVA_PEP_data.xlsx") %>% map(~read_xlsx("ANOVA_PEP_data.xlsx",.))
ANOVA_PEP_data
df <- ANOVA_PEP_data[[1]]
View(df)

# Step 2: Prepare data
df$Outcome <- as.factor(df$Outcome)
protein_cols <- colnames(df)[-(1:2)]  # assuming col1=Protein_ID, col2=Outcome

# Step 3: Initialize data frame to store ANOVA p-values
anova_results <- data.frame(Protein = character(),
                            anova_pvalue = numeric(),
                            stringsAsFactors = FALSE)

# Step 4: Loop over each protein for ANOVA only
for (protein in protein_cols) {
  safe_protein <- paste0("`", protein, "`")
  formula <- as.formula(paste(safe_protein, "~ Outcome"))
  
  res <- tryCatch({
    fit <- aov(formula, data = df)
    anova_summary <- summary(fit)
    anova_pvalue <- anova_summary[[1]][["Pr(>F)"]][1]
    
    data.frame(Protein = protein, anova_pvalue = anova_pvalue)
  }, error = function(e) {
    message(paste("Failed for protein:", protein))
    NULL
  })
  
  if (!is.null(res)) {
    anova_results <- bind_rows(anova_results, res)
  }
}

# Step 5: Save final results to CSV
write.csv(anova_results, "PEP_Stroke_Subtypes_ANOVA_pvalues.csv", row.names = FALSE)

#############################################################################################################################################
# Box plots of DEPs validated across SomaScan and DIA-MS in pairwise comparisons

# Load packages
library(readxl)
library(ggplot2)
library(dplyr)
library(purrr)   # for map()

# Step 1: Load Excel file
# AIS vs. ICH
excel_sheets("AIS_ICH_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("AIS_ICH_Validated_DEPs.xlsx") %>% map(~ read_xlsx("AIS_ICH_Validated_DEPs.xlsx", .))
Peptide_data

# AIS vs. TIA
excel_sheets("AIS_TIA_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("AIS_TIA_Validated_DEPs.xlsx") %>% map(~ read_xlsx("AIS_TIA_Validated_DEPs.xlsx", .))
Peptide_data

# AIS vs. MIM
excel_sheets("AIS_MIM_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("AIS_MIM_Validated_DEPs.xlsx") %>% map(~ read_xlsx("AIS_MIM_Validated_DEPs.xlsx", .))
Peptide_data

# ICH vs. TIA
excel_sheets("ICH_TIA_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("ICH_TIA_Validated_DEPs.xlsx") %>% map(~ read_xlsx("ICH_TIA_Validated_DEPs.xlsx", .))
Peptide_data

# ICH vs. MIM
excel_sheets("ICH_MIM_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("ICH_MIM_Validated_DEPs.xlsx") %>% map(~ read_xlsx("ICH_MIM_Validated_DEPs.xlsx", .))
Peptide_data

# TIA vs. MIM
excel_sheets("TIA_MIM_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("TIA_MIM_Validated_DEPs.xlsx") %>% map(~ read_xlsx("TIA_MIM_Validated_DEPs.xlsx", .))
Peptide_data

df <- Peptide_data[[2]]
names(df)

# Step 2: Preserve Group order as in the Excel file
# Use the first occurrence of each group to define levels
group_levels <- unique(df$Group)
df$Group <- factor(df$Group, levels = group_levels)

# Step 3: Identify protein columns (from column 5 onward)
protein_cols <- names(df)[5:ncol(df)]

# Step 4: Loop through each protein and create individual boxplots
for (protein in protein_cols) {
  
  # Create plot with white background
  p <- ggplot(df, aes(x = Group, y = .data[[protein]], fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste("Expression of", protein, "by Group"),
      x = "Group",
      y = "Expression Level"
    )
  
  # Step 5: Save each plot as PNG and EPS (white background)
  ggsave(
    filename = paste0("boxplot_", protein, ".png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = paste0("boxplot_", protein, ".eps"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white",
    device = cairo_ps
  )
  
  # Print progress
  message("Saved boxplot for: ", protein)
}

#############################################################################################################################################
# Box plots of individual proteins to group their multiple isoforms into a single plot (Pairwise Comparisons)

# Load packages
library(readxl)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

# Step 1: Load Excel file
excel_sheets("AIS_TIA_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("AIS_TIA_Validated_DEPs.xlsx") %>% map(~ read_xlsx("AIS_TIA_Validated_DEPs.xlsx", .))

excel_sheets("ICH_TIA_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("ICH_TIA_Validated_DEPs.xlsx") %>% map(~ read_xlsx("ICH_TIA_Validated_DEPs.xlsx", .))

df <- Peptide_data[[2]]
names(df)

# Step 2: Identify protein columns (from column 5 onward)
protein_cols <- names(df)[5:ncol(df)]

# Step 3: Specify proteins of interest
proteins_to_plot <- c("F13A1.GTYIPVPIVSELQSGK", "F13B.KTEEVECLTYGWSLTPK")   

# Filter to ensure these proteins exist in the dataset
proteins_to_plot <- intersect(proteins_to_plot, protein_cols)

# Step 4: Reshape data to long format for selected proteins
df_long <- df %>%
  select(Group, all_of(proteins_to_plot)) %>%
  pivot_longer(cols = all_of(proteins_to_plot),
               names_to = "Protein",
               values_to = "Expression")

# Step 4b: Order proteins by mean expression (highest to lowest)
protein_order <- df_long %>%
  group_by(Protein) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE)) %>%
  arrange(desc(mean_expr)) %>%
  pull(Protein)

df_long$Protein <- factor(df_long$Protein, levels = protein_order)

# Step 5: Create single combined box plot
p <- ggplot(df_long, aes(x = Protein, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7, position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = paste("Expression of Selected Proteins by Group"),
    x = "Protein",
    y = "Expression Level"
  )

# Step 6: Save the combined plot as PNG and EPS (white background)
ggsave(
  filename = "boxplot_selected_proteins.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "boxplot_selected_proteins.eps",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white",
  device = cairo_ps
)

message("Saved combined boxplot for selected proteins, ordered by mean expression.")

####################################################################################################################################
# Box plots of DEPs validated across Soma and DIA-MS in multigroup comparisons

# Load packages
library(readxl)
library(ggplot2)
library(dplyr)
library(purrr)   # for map()

# Step 1: Load Excel file
# Multigroup comparison file upload
excel_sheets("Multigroup_Validated_DEPs.xlsx")
Peptide_data <- excel_sheets("Multigroup_Validated_DEPs.xlsx") %>% map(~ read_xlsx("Multigroup_Validated_DEPs.xlsx", .))
Peptide_data

df <- Peptide_data[[2]]
names(df)

# Step 2: Preserve Group order as in the Excel file
# Use the first occurrence of each group to define levels
group_levels <- unique(df$Group)
df$Group <- factor(df$Group, levels = group_levels)

# Step 3: Identify protein columns (from column 3 onward)
protein_cols <- names(df)[3:ncol(df)]

# Step 4: Loop through each protein and create individual boxplots
for (protein in protein_cols) {
  
  # Create plot with white background
  p <- ggplot(df, aes(x = Group, y = .data[[protein]], fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste("Expression of", protein, "by Group"),
      x = "Group",
      y = "Expression Level"
    )
  
  # Step 5: Save each plot as PNG and EPS (white background)
  ggsave(
    filename = paste0("boxplot_", protein, ".png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = paste0("boxplot_", protein, ".eps"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white",
    device = cairo_ps
  )
  
  # Print progress
  message("Saved boxplot for: ", protein)
}

###########################################################################################################################################

# Heatmap and PCA plot of 40 validated protein peptides using DIA-MS in an external validation cohort

#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library (readxl)
library(writexl)
library(gplots)
library(ggrepel)
library(factoextra)
library(pROC)
library(caret)
library(rsample)
library(vip)
library(Matrix)
library(BiocParallel)

#Proteomics data file upload
excel_sheets("Validated_protein_peptides.xlsx")
Val_data= excel_sheets("Validated_protein_peptides.xlsx") %>% map(~read_xlsx("Validated_protein_peptides.xlsx",.))
Val_data

#Clustering and Dendograms

#Convert to matrix data frame
dat_matrix = as.matrix.data.frame(Val_data[[5]][,2:81]) 
#Name the rows with protein ids
row.names(dat_matrix) <- Val_data[[5]]$Protein_ID
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

#Create a heatmap using pheatmap package
dev.off() #Use this if pheatmap function is not working

##############################################################################################################################################
#PCA plot
##############################################################################################################################################

#PCA plot using 40 validated protein peptides
R.PCA<- prcomp(Val_data[[4]][,3:42], scale=TRUE)
R.PCA
X11()
pca_plot=fviz_pca_ind(R.PCA, col.ind=Val_data[[4]]$Outcome, title= "PCA plot of Validated peptides using DIA-MS", addEllipses = FALSE,
                      label= "none", pointsize= 4)

# Define colors for the groups
colors = c("AIS" = "red", "ICH" = "blue", "TIA" = "purple", "MIM"= "green")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

# Plot PC2 vs. PC3
pca_plot <- fviz_pca_ind(
  R.PCA,
  axes = c(2, 3),  # specify which PCs to plot
  col.ind = Val_data[[6]]$Outcome,
  title = "PCA plot (PC2 vs PC3) of Validated peptides using DIA-MS",
  addEllipses = FALSE,
  label = "none",
  pointsize = 4
)

# Define custom colors
colors <- c("AIS" = "red", "ICH" = "blue", "TIA" = "purple", "MIM" = "green")

# Apply custom colors
pca_plot + scale_color_manual(values = colors)