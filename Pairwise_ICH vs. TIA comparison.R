#Pairwise comparisons of Stroke subtypes

#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library (readxl)
library(writexl)
library(gplots)
library(ggrepel)
library(factoextra)
library(Boruta)
library(randomForest)
library(pROC)
library(caret)
library(rsample)
library(vip)
library(variancePartition)
library(Matrix)
library(BiocParallel)

#ICH vs. TIA comparison (SomaScan Proteomics)

#Proteomics data file upload
excel_sheets("ICH_TIA_data.xlsx")
ICH_TIA_data= excel_sheets("ICH_TIA_data.xlsx") %>% map(~read_xlsx("ICH_TIA_data.xlsx",.))
ICH_TIA_data

#Visualizing the Data
dat = ICH_TIA_data[[1]]
View(dat)

#Creating a T-test function for multiple experiments
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

#Apply t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:41)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:41]))$p.value

#Plot histogram of p-values
dat_pvals %>% 
  ggplot(aes(p_val)) + 
  geom_histogram(binwidth = 0.05, 
                 boundary = 0.5, 
                 fill = "darkblue",
                 colour = "white") +
  xlab("p-value") +
  ylab("Frequency") +
  theme_minimal()

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,42])
View (dat_combine)

#Calculating log-fold change
dat_fc = dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_ICH_case = mean(c(ICH1, ICH2, ICH3, ICH4, ICH5, ICH6, ICH7, ICH8, ICH9, ICH10, ICH11, ICH12, ICH13, ICH14, ICH15, ICH16, ICH17,
                                       ICH18, ICH19, ICH20)),
                mean_TIA_case= mean(c(TIA1, TIA2, TIA3, TIA4, TIA5, TIA6, TIA7, TIA8, TIA9, TIA10, TIA11, TIA12, TIA13, TIA14, TIA15, TIA16,
                                      TIA17, TIA18, TIA19, TIA20)),
                log_fc = mean_ICH_case - mean_TIA_case,
                log_pval = -1*log10(p_val))
View(dat_fc)

#Save final data with list of final data in csv file
write.csv(dat_fc, "Final_ICH_TIA_data.csv")

#Visualize transformed data
#Plot a histogram to look at the distribution of log-fold change
dat_fc %>%
  ggplot(aes(log_fc)) + 
  geom_histogram(binwidth = 0.5,
                 boundary = 0.5,
                 fill = "darkblue",
                 colour = "white") +
  xlab("log2 fold change") +
  ylab("Frequency") +
  theme_minimal()

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
  select(Protein_ID, ICH1:ICH20, TIA1:TIA20, mean_ICH_case, mean_TIA_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in csv file
write.csv(final_data, "Final_diffproteins_ICH_TIA.csv")


#PCA plot using differentially expressed proteins
R.PCA = prcomp(ICH_TIA_data[[2]][,3:37], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=ICH_TIA_data[[2]]$Outcome, title= "PCA plot of differentially expressed proteins- ICH vs. TIA", addEllipses = FALSE,
                       label= "none", pointsize= 4)
# Define colors for the groups
colors = c("ICH" = "red", "TIA" = "blue")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

#Random forest- based classification to test Accuracy of differentially expressed proteins

## Set seed for reproducibility
set.seed(123)

## Define repeated cross validation with 10 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=3)

## Set seed for reproducibility
set.seed(123)

## Split the data into equal proportions so that we use 70% of it for training
index.data <- initial_split(ICH_TIA_data[[2]][,-1], prop = 0.7, strata = "Outcome")
View(index.data)
train.data <- training(index.data)
test.data  <- testing(index.data)
View(train.data)
View(test.data)
train.data$Outcome= factor(train.data$Outcome)
test.data$Outcome= factor(test.data$Outcome)

## Set seed for reproducibility
set.seed(444)

## Train a random forest model
forest <- train(Outcome~., data=train.data, method='rf', trControl=repeat_cv, metric='Accuracy', ntree= 1000, keep.forest=TRUE)
confusionMatrix(forest)

plot(forest)
forest$coefnames
forest$bestTune

#Variable importance
X11()
vip(forest, num_features = 10, geom= c("point"))

## Generate predictions for the test set
test_predictions <- predict(
  ## Random forest object
  object=forest, 
  ## Data to use for predictions; remove Outcome
  newdata=test.data[, -1])

## Print the accuracy
accuracy = mean(test_predictions == test.data$Outcome)*100
cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')
#Predictions
pred = predict(forest, newdata=test.data, type="prob")
#Generate ROC curve
roc_curve = roc(test.data$Outcome, pred[, "ICH"])
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

##################################################################################################################################

#Boruta feature selection-based random forest method

set.seed(222)
Boruta.ICH.TIA = Boruta(Outcome~.-Protein_ID, data = ICH_TIA_data[[3]], doTrace = 2, maxRuns=1000)
print(Boruta.ICH.TIA)
#Take a call on tentative features
Boruta.ICH.TIA.all = TentativeRoughFix(Boruta.ICH.TIA)
print(Boruta.ICH.TIA.all)
Sign.boruta= names(Boruta.ICH.TIA.all$finalDecision[Boruta.ICH.TIA.all$finalDecision %in% c("Confirmed")])
print(Sign.boruta)

#Plot Confirmed Boruta Attributes only
#generateCol is needed by plot.Boruta
generateCol<-function(x,colCode,col,numShadow){
  #Checking arguments
  if(is.null(col) & length(colCode)!=4)
    stop('colCode should have 4 elements.');
  #Generating col
  if(is.null(col)){
    rep(colCode[4],length(x$finalDecision)+numShadow)->cc;
    cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1];
    cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2];
    cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3];
    col=cc;
  }
  return(col);
}

# Modified plot.Boruta
plot.Boruta.sel <- function(
    x,
    pars = NULL,
    colCode = c('green','yellow','red','blue'),
    sort = TRUE,
    whichShadow = c(TRUE, TRUE, TRUE),
    col = NULL, xlab = 'Attributes', ylab = 'Importance', ...) {
  
  #Checking arguments
  if(class(x)!='Boruta')
    stop('This function needs Boruta object as an argument.');
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.');
  
  #Removal of -Infs and conversion to a list
  lz <- lapply(1:ncol(x$ImpHistory), function(i)
    x$ImpHistory[is.finite(x$ImpHistory[,i]),i]);
  colnames(x$ImpHistory)->names(lz);
  
  #Selection of shadow meta-attributes
  numShadow <- sum(whichShadow);
  lz <- lz[c(rep(TRUE,length(x$finalDecision)), whichShadow)];
  
  #Generating color vector
  col <- generateCol(x, colCode, col, numShadow);
  
  #Ordering boxes due to attribute median importance
  if (sort) {
    ii <- order(sapply(lz, stats::median));
    lz <- lz[ii];
    col <- col[ii];
  }
  
  # Select parameters of interest
  if (!is.null(pars)) lz <- lz[names(lz) %in% pars];
  
  #Final plotting
  graphics::boxplot(lz, xlab = xlab, ylab = ylab, col = "green", ...);
  invisible(x);
}

#Boruta plot of confirmed features
X11()
plot.Boruta.sel(Boruta.ICH.TIA.all, cex.axis= 0.5, las= 2, xlab= "",
                pars = c("CLIC5.12475-48", "HSPA9.13492-44",  "S100A10.15318-75", "SQSTM1.15448-47",  "MSRB2.16877-19", "PRSS3.18864-7",  "AMY2A.18917-53", 
                         "CDH12.19388-2", "CD1B.20541-73",  "UBL5.21104-37",  "DUSP21.21191-24",  "BTG1.22088-3", "NGDN.22512-18",  "NPEPL1.23639-93",  
                         "PANK1.25277-51",  "SYT12.25424-234", "PRSS1.3049-61", "FN1.3434-34",  "CDKN1B.3719-2",  "PLG.4151-6", "HRG.4996-66",  "NCK1.5014-49",     
                         "S100A9.5339-49",  "CD27.5412-53", "CLEC3B.5701-81", "TAGLN2.9201-13", "NECTIN1.9300-13"));

#Get only selected features
getSelectedAttributes(Boruta.ICH.TIA.all, withTentative = F)
boruta.df = attStats(Boruta.ICH.TIA.all)
class(boruta.df)
print(boruta.df)

write.csv(boruta.df, "Boruta_results_ICH_TIA.csv")

#PCA plot using Confirmed Boruta Features
R.PCA<- prcomp(ICH_TIA_data[[4]][,3:29], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=ICH_TIA_data[[4]]$Outcome, title= "PCA plot of Confirmed Boruta features- ICH vs. TIA", addEllipses = FALSE, 
                       label= "none", pointsize= 4)
# Define colors for the groups
colors = c("ICH" = "red", "TIA" = "blue")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

#Random forest- based classification to test Accuracy of Boruta Confirmed Features

## Set seed for reproducibility
set.seed(123)

## Define repeated cross validation with 10 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=3)

## Set seed for reproducibility
set.seed(123)

## Split the data into equal proportions so that we use 70% of it for training
index.data <- initial_split(ICH_TIA_data[[4]][,-1], prop = 0.7, strata = "Outcome")
View(index.data)
train.data <- training(index.data)
test.data  <- testing(index.data)
View(train.data)
View(test.data)
train.data$Outcome= factor(train.data$Outcome)
test.data$Outcome= factor(test.data$Outcome)

## Set seed for reproducibility
set.seed(111)

## Train a random forest model
forest <- train(Outcome~., data=train.data, method='rf', trControl=repeat_cv, metric='Accuracy', ntree= 1000, keep.forest=TRUE)
confusionMatrix(forest)

plot(forest)
forest$coefnames
forest$bestTune

#Variable importance
X11()
vip(forest, num_features = 10, geom= c("point"))

## Generate predictions for the test set
test_predictions <- predict(
  ## Random forest object
  object=forest, 
  ## Data to use for predictions; remove Outcome
  newdata=test.data[, -1])

## Print the accuracy
accuracy = mean(test_predictions == test.data$Outcome)*100
cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')
#Predictions
pred = predict(forest, newdata=test.data, type="prob")
#Generate ROC curve
roc_curve = roc(test.data$Outcome, pred[, "ICH"])
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

####################################################################################################################################################################

#Variance Partitioning Analysis (VPA)

#Load the expression dataset
Exp.data= read.csv("Expression_data_ICH_TIA.csv", header = TRUE, row.names = 1)
View(Exp.data)

#Load the meta dataset
Meta.data= read.csv("Meta_data_ICH_TIA.csv", header = TRUE, row.names = 1)
View(Meta.data)

#Convert Expression dataset into a matrix
Exp.data= as.matrix(Exp.data)
class(Exp.data)

#Convert meta data into a data frame
regvars.vp= data.frame(Meta.data)

#Convert Continuous variables to numeric
regvars.vp$Age<-as.numeric(regvars.vp$Age)
regvars.vp$NIHSS0<-as.numeric(regvars.vp$NIHSS0)
regvars.vp$GCS<-as.numeric(regvars.vp$GCS)

#Convert Discrete variables to numeric
regvars.vp$Sex<-factor(regvars.vp$Sex)
regvars.vp$Outcome<-factor(regvars.vp$Outcome)
regvars.vp$Batch<-factor(regvars.vp$Batch)
regvars.vp$Smoking<-factor(regvars.vp$Smoking)
regvars.vp$AFIB<-factor(regvars.vp$AFIB)
regvars.vp$DM<-factor(regvars.vp$DM)
regvars.vp$HTN<-factor(regvars.vp$HTN)
regvars.vp$CAD<-factor(regvars.vp$CAD)

#Formula for indicating which meta data variables to use in the analysis
form = ~ Age + NIHSS0 + GCS + (1|Outcome) + (1|Sex) + (1|DM) + (1|HTN) + (1|CAD) + (1|AFIB) + (1|Smoking) + (1|Batch)  

#Fitting a linear mixed model to extract variance fractions
varPart = fitExtractVarPartModel(Exp.data, form, regvars.vp)
head(varPart)

#Violin plot of contribution of each variable to total variance
#Check if Batch effect is present
plotVarPart(varPart, main="Violin Plot for total variance in proteins between ICH and TIA")

#Remove Batch effect if present
#extract residuals directly without storing intermediate results
residList = fitVarPartModel(Exp.data, ~ (1 | Batch), regvars.vp, fxn = residuals)

#Convert Batch corrected list to matrix
residMatrix <- do.call(rbind, residList)
View(residMatrix)

#Formula for indicating which meta data variables to use in the analysis (Batch removed now)
form = ~ Age + NIHSS0 + GCS + (1|Outcome) + (1|Sex) + (1|DM) + (1|HTN) + (1|CAD) + (1|AFIB) + (1|Smoking)

#Again fit the linear mixed model to extract variance fractions
varPartResid = fitExtractVarPartModel(residMatrix, form, regvars.vp)
head(varPartResid)

#Sort variables (i.e. columns) in a given order
vp = sortCols(varPartResid, last= c("Outcome", "Age", "Sex", "DM", "HTN", "CAD", "AFIB", "Smoking", "NIHSS0", "GCS", "Residuals"))
vp

#Violin plot of contribution of each variable to total variance
plotVarPart(vp, main="Violin Plot for total variance in proteins between ICH and TIA")

#Access first entries for outcome
head(vp$Outcome)

#Sort proteins based on variance explained by Outcome
head(vp[order(vp$Outcome, decreasing = TRUE), ])
vpa.outcome= vp[order(vp$Outcome, decreasing = TRUE), ]
vpa.outcome
write.csv(vpa.outcome, "Protein variance_ICH_TIA.csv")

#Sort proteins by variance explained by Outcome 
OutcomeSortOrder<-order(vp[["Outcome"]],decreasing=TRUE)  
for (i in ls(vp)) { vp[[i]]<-vp[[i]][OutcomeSortOrder]; }
rownames(vp)<-rownames(vp)[OutcomeSortOrder]
X11()
#Plot top proteins associated with outcome
plotPercentBars( vp[1:6,]) + ggtitle( "Proteins with 20% or more variance explained by Outcome (ICH vs. TIA)" )


#Proteins with 20% or more variance explained by stroke subtypes

#PCA plot using VPA proteins
R.PCA<- prcomp(ICH_TIA_data[[5]][,3:8], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=ICH_TIA_data[[5]]$Outcome, title= "PCA plot of proteins with 20% or more variation- ICH vs. TIA", addEllipses = FALSE, label= "none",
                       pointsize= 4)
# Define colors for the groups
colors = c("ICH" = "red", "TIA" = "blue")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

#Random forest- based classification to test Accuracy of proteins with 20% or more variance explained by stroke subtypes

## Set seed for reproducibility
set.seed(123)

## Define repeated cross validation with 10 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=3)

## Set seed for reproducibility
set.seed(123)

## Split the data into equal proportions so that we use 70% of it for training
index.data <- initial_split(ICH_TIA_data[[5]][,-1], prop = 0.7, strata = "Outcome")
View(index.data)
train.data <- training(index.data)
test.data  <- testing(index.data)
View(train.data)
View(test.data)
train.data$Outcome= factor(train.data$Outcome)
test.data$Outcome= factor(test.data$Outcome)

## Set seed for reproducibility
set.seed(333)

## Train a random forest model
forest <- train(Outcome~., data=train.data, method='rf', trControl=repeat_cv, metric='Accuracy', ntree= 1000, keep.forest=TRUE)
confusionMatrix(forest)

print(forest)
plot(forest)
forest$coefnames
forest$finalModel
forest$bestTune

#Variable importance
vip(forest, num_features = 10, geom= c("point"))

## Generate predictions for the test set
test_predictions <- predict(
  ## Random forest object
  object=forest, 
  ## Data to use for predictions; remove Outcome
  newdata=test.data[, -1])

## Print the accuracy
accuracy = mean(test_predictions == test.data$Outcome)*100
cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')
#Predictions
pred = predict(forest, newdata=test.data, type="prob")
#Generate ROC curve
roc_curve = roc(test.data$Outcome, pred[, "ICH"])
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

####################################################################################################################################################################################################################

# Top Proteins

#PCA plot of top proteins
R.PCA<- prcomp(ICH_TIA_data[[7]][,3:6], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=ICH_TIA_data[[7]]$Outcome, title= "PCA plot of top proteins- ICH vs. TIA", addEllipses = FALSE, label= "none",
                       pointsize= 4)
# Define colors for the groups
colors = c("ICH" = "red", "TIA" = "blue")
# Add custom colors
pca_plot + scale_color_manual(values = colors)

#Random forest- based classification to test Accuracy of top proteins

## Set seed for reproducibility
set.seed(123)

## Define repeated cross validation with 10 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=3)

## Set seed for reproducibility
set.seed(123)

## Split the data into equal proportions so that we use 70% of it for training
index.data <- initial_split(ICH_TIA_data[[7]][,-1], prop = 0.7, strata = "Outcome")
View(index.data)
train.data <- training(index.data)
test.data  <- testing(index.data)
View(train.data)
View(test.data)
train.data$Outcome= factor(train.data$Outcome)
test.data$Outcome= factor(test.data$Outcome)

## Set seed for reproducibility
set.seed(123)

## Train a random forest model
forest <- train(Outcome~., data=train.data, method='rf', trControl=repeat_cv, metric='Accuracy', ntree= 1000, keep.forest=TRUE)
confusionMatrix(forest)

print(forest)
plot(forest)
forest$coefnames
forest$finalModel
forest$bestTune

#Variable importance
vip(forest, num_features = 10, geom= c("point"))

## Generate predictions for the test set
test_predictions <- predict(
  ## Random forest object
  object=forest, 
  ## Data to use for predictions; remove Outcome
  newdata=test.data[, -1])

## Print the accuracy
accuracy = mean(test_predictions == test.data$Outcome)*100
cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')
#Predictions
pred = predict(forest, newdata=test.data, type="prob")
#Generate ROC curve
roc_curve = roc(test.data$Outcome, pred[, "ICH"])
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

#Cross-platform Validation
#ICH vs. TIA comparison (BAK-270 Proteomics)

#Proteomics data file upload
excel_sheets("ICH_TIA_Bak.xlsx")
ICH_TIA_data= excel_sheets("ICH_TIA_Bak.xlsx") %>% map(~read_xlsx("ICH_TIA_Bak.xlsx",.))
ICH_TIA_data

#Visualizing the Data
dat = ICH_TIA_data[[1]]
View(dat)

#Creating a T-test function for multiple experiments
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

#Apply t-test function to data using plyr adply
#.margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:21), grp2 = c(22:40)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:21]), as.numeric(dat[1,22:40]))$p.value

#Plot histogram of p-values
dat_pvals %>% 
  ggplot(aes(p_val)) + 
  geom_histogram(binwidth = 0.05, 
                 boundary = 0.5, 
                 fill = "darkblue",
                 colour = "white") +
  xlab("p-value") +
  ylab("Frequency") +
  theme_minimal()

#Bind columns to create transformed data frame
dat_combine = bind_cols(dat, dat_pvals[,41])
View (dat_combine)

#Calculating log-fold change
dat_fc = dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_ICH_case = mean(c(ICH1, ICH2, ICH3, ICH4, ICH5, ICH6, ICH7, ICH8, ICH9, ICH10, ICH11, ICH12, ICH13, ICH14, ICH15, ICH16, ICH17,
                                       ICH18, ICH19, ICH20)),
                mean_TIA_case= mean(c(TIA1, TIA2, TIA3, TIA4, TIA5, TIA6, TIA7, TIA8, TIA9, TIA10, TIA11, TIA12, TIA13, TIA14, TIA15, TIA16,
                                      TIA17, TIA18, TIA19)),
                log_fc = mean_ICH_case - mean_TIA_case,
                log_pval = -1*log10(p_val))
View(dat_fc)

#Save final data with list of final data in csv file
write.csv(dat_fc, "Final_ICH_TIA_Bak.csv")

#Visualize transformed data
#Plot a histogram to look at the distribution of log-fold change
dat_fc %>%
  ggplot(aes(log_fc)) + 
  geom_histogram(binwidth = 0.5,
                 boundary = 0.5,
                 fill = "darkblue",
                 colour = "white") +
  xlab("log2 fold change") +
  ylab("Frequency") +
  theme_minimal()

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
select(Protein_ID, ICH1:ICH20, TIA1:TIA19, mean_ICH_case, mean_TIA_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in csv file
write.csv(final_data, "Final_diffproteins_ICH_TIA_Bak.csv")

##################################################################################################################################

#Boruta feature selection-based random forest method

set.seed(777)
Boruta.ICH.TIA = Boruta(Outcome~.-Protein_ID, data = ICH_TIA_data[[2]], doTrace = 2, maxRuns=1000)
print(Boruta.ICH.TIA)
#Take a call on tentative features
Boruta.ICH.TIA.all = TentativeRoughFix(Boruta.ICH.TIA)
print(Boruta.ICH.TIA.all)
Sign.boruta= names(Boruta.ICH.TIA.all$finalDecision[Boruta.ICH.TIA.all$finalDecision %in% c("Confirmed")])
print(Sign.boruta)

#Plot Confirmed Boruta Attributes only
#generateCol is needed by plot.Boruta
generateCol<-function(x,colCode,col,numShadow){
  #Checking arguments
  if(is.null(col) & length(colCode)!=4)
    stop('colCode should have 4 elements.');
  #Generating col
  if(is.null(col)){
    rep(colCode[4],length(x$finalDecision)+numShadow)->cc;
    cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1];
    cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2];
    cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3];
    col=cc;
  }
  return(col);
}

# Modified plot.Boruta
plot.Boruta.sel <- function(
    x,
    pars = NULL,
    colCode = c('green','yellow','red','blue'),
    sort = TRUE,
    whichShadow = c(TRUE, TRUE, TRUE),
    col = NULL, xlab = 'Attributes', ylab = 'Importance', ...) {
  
  #Checking arguments
  if(class(x)!='Boruta')
    stop('This function needs Boruta object as an argument.');
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.');
  
  #Removal of -Infs and conversion to a list
  lz <- lapply(1:ncol(x$ImpHistory), function(i)
    x$ImpHistory[is.finite(x$ImpHistory[,i]),i]);
  colnames(x$ImpHistory)->names(lz);
  
  #Selection of shadow meta-attributes
  numShadow <- sum(whichShadow);
  lz <- lz[c(rep(TRUE,length(x$finalDecision)), whichShadow)];
  
  #Generating color vector
  col <- generateCol(x, colCode, col, numShadow);
  
  #Ordering boxes due to attribute median importance
  if (sort) {
    ii <- order(sapply(lz, stats::median));
    lz <- lz[ii];
    col <- col[ii];
  }
  
  # Select parameters of interest
  if (!is.null(pars)) lz <- lz[names(lz) %in% pars];
  
  #Final plotting
  graphics::boxplot(lz, xlab = xlab, ylab = ylab, col = "green", ...);
  invisible(x);
}

#Boruta plot of confirmed features
plot.Boruta.sel(Boruta.ICH.TIA.all, cex.axis= 0.5, las= 2, xlab= "",
                pars = c("CRTAC1", "C1RL",  "TNXB", "PARK7",  "PRL"));


#Get only selected features
getSelectedAttributes(Boruta.ICH.TIA.all, withTentative = F)
boruta.df = attStats(Boruta.ICH.TIA.all)
class(boruta.df)
print(boruta.df)

write.csv(boruta.df, "Boruta_results_ICH_TIA_Bak.csv")

#############################################################################################################################################################

#Variance Partitioning Analysis (VPA)

#Load the expression dataset
Exp.data= read.csv("Expression_data_ICH_TIA_Bak.csv", header = TRUE, row.names = 1)
View(Exp.data)

#Load the meta dataset
Meta.data= read.csv("Meta_data_ICH_TIA_Bak.csv", header = TRUE, row.names = 1)
View(Meta.data)

#Convert Expression dataset into a matrix
Exp.data= as.matrix(Exp.data)
class(Exp.data)

#Convert meta data into a data frame
regvars.vp= data.frame(Meta.data)

#Convert Continuous variables to numeric
regvars.vp$Age<-as.numeric(regvars.vp$Age)
regvars.vp$NIHSS0<-as.numeric(regvars.vp$NIHSS0)
regvars.vp$GCS<-as.numeric(regvars.vp$GCS)

#Convert Discrete variables to numeric
regvars.vp$Sex<-factor(regvars.vp$Sex)
regvars.vp$Outcome<-factor(regvars.vp$Outcome)
regvars.vp$Batch<-factor(regvars.vp$Batch)
regvars.vp$Smoking<-factor(regvars.vp$Smoking)
regvars.vp$AFIB<-factor(regvars.vp$AFIB)
regvars.vp$DM<-factor(regvars.vp$DM)
regvars.vp$HTN<-factor(regvars.vp$HTN)
regvars.vp$CAD<-factor(regvars.vp$CAD)

#Formula for indicating which meta data variables to use in the analysis
form = ~ Age + NIHSS0 + GCS + (1|Outcome) + (1|Sex) + (1|DM) + (1|HTN) + (1|CAD) + (1|AFIB) + (1|Smoking) + (1|Batch)  

#Fitting a linear mixed model to extract variance fractions
varPart = fitExtractVarPartModel(Exp.data, form, regvars.vp)
head(varPart)

#Violin plot of contribution of each variable to total variance
#Check if Batch effect is present
plotVarPart(varPart, main="Violin Plot for total variance in proteins between AIS and TIA")

#Remove Batch effect if present
#extract residuals directly without storing intermediate results
residList = fitVarPartModel(Exp.data, ~ (1 | Batch), regvars.vp, fxn = residuals)

#Convert Batch corrected list to matrix
residMatrix <- do.call(rbind, residList)
View(residMatrix)

#Formula for indicating which meta data variables to use in the analysis (Batch removed now)
form = ~ Age + NIHSS0 + GCS + (1|Outcome) + (1|Sex) + (1|DM) + (1|HTN) + (1|CAD) + (1|AFIB) + (1|Smoking)

#Again fit the linear mixed model to extract variance fractions
varPartResid <- fitExtractVarPartModel(residMatrix, form, regvars.vp)
head(varPartResid)

#Sort variables (i.e. columns) in a given order
vp = sortCols(varPartResid, last= c("Outcome", "Age", "Sex", "DM", "HTN", "CAD", "AFIB", "Smoking", "NIHSS0", "GCS", "Residuals"))
vp

#Violin plot of contribution of each variable to total variance
plotVarPart(vp, main="Violin Plot for total variance in proteins between AIS and TIA")

#Access first entries for outcome
head(vp$Outcome)

#Sort proteins based on variance explained by Outcome
head(vp[order(vp$Outcome, decreasing = TRUE), ])
vpa.outcome= vp[order(vp$Outcome, decreasing = TRUE), ]
vpa.outcome
write.csv(vpa.outcome, "Protein variance_Outcome_ICH_TIA_Bak.csv")

#Sort proteins by variance explained by Outcome 
OutcomeSortOrder<-order(vp[["Outcome"]],decreasing=TRUE)  
for (i in ls(vp)) { vp[[i]]<-vp[[i]][OutcomeSortOrder]; }
rownames(vp)<-rownames(vp)[OutcomeSortOrder]
#Plot top proteins associated with outcome
plotPercentBars( vp[1:10,]) + ggtitle( "Top 10 Proteins associated with Outcome (ICH vs. TIA)" )