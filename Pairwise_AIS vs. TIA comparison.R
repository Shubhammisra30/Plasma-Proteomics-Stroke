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

#AIS vs. TIA comparison (SomaScan Proteomics)

#Proteomics data file upload
excel_sheets("IS_TIA_data.xlsx")
IS_TIA_data= excel_sheets("IS_TIA_data.xlsx") %>% map(~read_xlsx("IS_TIA_data.xlsx",.))
IS_TIA_data

#Visualizing the Data
dat = IS_TIA_data[[1]]
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
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:41), grp2 = c(42:61)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:41]), as.numeric(dat[1,42:61]))$p.value

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
dat_combine = bind_cols(dat, dat_pvals[,62])
View (dat_combine)

#Calculating log-fold change
dat_fc = dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_IS_case = mean(c(AIS1, AIS2, AIS3, AIS4, AIS5, AIS6, AIS7, AIS8, AIS9, AIS10, AIS11, AIS12, AIS13, AIS14, AIS15, AIS16,
                                      AIS17, AIS18, AIS19, AIS20, AIS21, AIS22, AIS23, AIS24, AIS25, AIS26, AIS27, AIS28, AIS29, AIS30, AIS31,
                                      AIS32, AIS33, AIS34, AIS35, AIS36, AIS37, AIS38, AIS39, AIS40)),
                mean_TIA_case= mean(c(TIA1, TIA2, TIA3, TIA4, TIA5, TIA6, TIA7, TIA8, TIA9, TIA10, TIA11, TIA12, TIA13, TIA14, TIA15, TIA16,
                                      TIA17, TIA18, TIA19, TIA20)),
                log_fc = mean_IS_case - mean_TIA_case,
                log_pval = -1*log10(p_val))
View(dat_fc)

#Save final data with list of final data in csv file
write.csv(dat_fc, "Final_IS_TIA_data.csv")

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
  select(Protein_ID, AIS1:AIS40, TIA1:TIA20, mean_IS_case, mean_TIA_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in csv file
write.csv(final_data, "Final_diffproteins_IS_TIA.csv")


#PCA plot using differentially expressed proteins
R.PCA<- prcomp(IS_TIA_data[[2]][,3:33], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=IS_TIA_data[[2]]$Outcome, title= "PCA plot of differentially expressed proteins- AIS vs TIA", addEllipses = TRUE, 
                       label= "none", pointsize= 4)
# Define colors for the groups
colors = c("AIS" = "red", "TIA" = "blue")
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
index.data <- initial_split(IS_TIA_data[[2]][,-1], prop = 0.7, strata = "Outcome")
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
roc_curve = roc(test.data$Outcome, pred[, "AIS"])
X11()
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

##################################################################################################################################################

#Boruta feature selection-based random forest method

set.seed(333)
Boruta.IS.TIA = Boruta(Outcome~.-Protein_ID, data = IS_TIA_data[[3]], doTrace = 2, maxRuns=1000)
print(Boruta.IS.TIA)
#Take a call on tentative features
Boruta.IS.TIA.all = TentativeRoughFix(Boruta.IS.TIA)
print(Boruta.IS.TIA.all)
Sign.boruta= names(Boruta.IS.TIA.all$finalDecision[Boruta.IS.TIA.all$finalDecision %in% c("Confirmed")])
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
plot.Boruta.sel(Boruta.IS.TIA.all, cex.axis= 0.5, las= 2, xlab= "",
                pars = c("DUSP4.10035-6", "TPMT.11218-84",  "P4HA2.11348-132",  "ZNF560.12543-76",  "ELF5.13457-33",  "HSPA9.13492-44", "MSRB2.16877-19", 
                         "GTF2F2.20927-43", "ITGA6.21979-8",  "ADGRB2.22588-35",  "SKP2.22777-46",  "LMO4.23284-5", "BICD1.25118-45", "SCCPDH.25436-44",  
                         "NEUROG3.25938-13", "PLG.4151-6",  "SIGLEC1.4464-10",  "CTRB2.5648-28",  "PTPRS.6049-64",  "DDX39B.9742-59", "MTAP.9910-9"));

#Get only selected features
getSelectedAttributes(Boruta.IS.TIA.all, withTentative = F)
boruta.df = attStats(Boruta.IS.TIA.all)
class(boruta.df)
print(boruta.df)

write.csv(boruta.df, "Boruta_results_IS_TIA.csv")

#PCA plot using Confirmed Boruta Features
R.PCA<- prcomp(IS_TIA_data[[4]][,3:23], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=IS_TIA_data[[4]]$Outcome, title= "PCA plot of Confirmed Boruta features- AIS vs TIA", addEllipses = TRUE, 
                       label= "none", pointsize= 4)
# Define colors for the groups
colors = c("AIS" = "red", "TIA" = "blue")
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
index.data <- initial_split(IS_TIA_data[[4]][,-1], prop = 0.7, strata = "Outcome")
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
pred
#Generate ROC curve
roc_curve = roc(test.data$Outcome, pred[, "AIS"])
X11()
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

######################################################################################################################################################

#Variance Partitioning Analysis (VPA)

#Load the expression dataset
Exp.data= read.csv("Expression_data_IS_TIA.csv", header = TRUE, row.names = 1)
View(Exp.data)

#Load the meta dataset
Meta.data= read.csv("Meta_data_IS_TIA.csv", header = TRUE, row.names = 1)
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
varPartResid = fitExtractVarPartModel(residMatrix, form, regvars.vp)
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
write.csv(vpa.outcome, "Protein variance_Outcome_IS_TIA.csv")

#Sort proteins by variance explained by Outcome 
OutcomeSortOrder<-order(vp[["Outcome"]],decreasing=TRUE)  
for (i in ls(vp)) { vp[[i]]<-vp[[i]][OutcomeSortOrder]; }
rownames(vp)<-rownames(vp)[OutcomeSortOrder]
#Plot top proteins associated with outcome
X11()
plotPercentBars( vp[1:2,]) + ggtitle( "Proteins associated with 20% or more variance explained by Outcome (AIS vs. TIA)" )


#Proteins with 20% or more variance explained by stroke subtypes

#PCA plot for VPA- IS vs TIA
R.PCA<- prcomp(IS_TIA_data[[5]][,3:4], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=IS_TIA_data[[5]]$Outcome, title= "PCA plot of proteins with 20% or more variation- AIS vs. TIA", addEllipses = FALSE, label= "none",
                       pointsize= 4)
# Define colors for the groups
colors = c("AIS" = "red", "TIA" = "blue")
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
index.data <- initial_split(IS_TIA_data[[5]][,-1], prop = 0.7, strata = "Outcome")
View(index.data)
train.data <- training(index.data)
test.data  <- testing(index.data)
View(train.data)
View(test.data)
train.data$Outcome= factor(train.data$Outcome)
test.data$Outcome= factor(test.data$Outcome)

## Set seed for reproducibility
set.seed(666)

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
roc_curve = roc(test.data$Outcome, pred[, "AIS"])
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

#############################################################################################################################################################

#Top proteins

#PCA plot of top proteins
R.PCA<- prcomp(IS_TIA_data[[5]][,3:6], scale=TRUE)
R.PCA
X11()
pca_plot= fviz_pca_ind(R.PCA, col.ind=IS_TIA_data[[5]]$Outcome, title= "PCA plot of top proteins- AIS vs TIA", addEllipses = TRUE, label= "none",
                       pointsize= 4)
# Define colors for the groups
colors = c("AIS" = "red", "TIA" = "blue")
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
index.data <- initial_split(IS_TIA_data[[5]][,-1], prop = 0.7, strata = "Outcome")
View(index.data)
train.data <- training(index.data)
test.data  <- testing(index.data)
View(train.data)
View(test.data)
train.data$Outcome= factor(train.data$Outcome)
test.data$Outcome= factor(test.data$Outcome)

## Set seed for reproducibility
set.seed(222)

## Train a random forest model
forest <- train(Outcome~., data=train.data, method='rf', trControl=repeat_cv, metric='Accuracy', ntree= 1000, keep.forest=TRUE)
confusionMatrix(forest)

plot(forest)
print(forest)
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
pred
#Generate ROC curve
roc_curve = roc(test.data$Outcome, pred[, "AIS"])
X11()
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
auc <- auc(roc_curve)
auc
auc_ci= ci.auc(roc_curve)
auc_ci

#############################################################################################################################################################
#############################################################################################################################################################

#Cross-platform Validation
#AIS vs. TIA comparison (BAK-270 Proteomics)

#Proteomics data file upload
excel_sheets("IS_TIA_Bak.xlsx")
IS_TIA_data= excel_sheets("IS_TIA_Bak.xlsx") %>% map(~read_xlsx("IS_TIA_Bak.xlsx",.))
IS_TIA_data

#Visualizing the Data
dat = IS_TIA_data[[1]]
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
dat_pvals = plyr::adply(dat,.margins = 1, .fun = t_test, grp1 = c(2:41), grp2 = c(42:60)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat[1,2:41]), as.numeric(dat[1,42:60]))$p.value

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
dat_combine = bind_cols(dat, dat_pvals[,61])
View (dat_combine)

#Calculating log-fold change
dat_fc = dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_IS_case = mean(c(AIS1, AIS2, AIS3, AIS4, AIS5, AIS6, AIS7, AIS8, AIS9, AIS10, AIS11, AIS12, AIS13, AIS14, AIS15, AIS16,
                                      AIS17, AIS18, AIS19, AIS20, AIS21, AIS22, AIS23, AIS24, AIS25, AIS26, AIS27, AIS28, AIS29, AIS30, AIS31,
                                      AIS32, AIS33, AIS34, AIS35, AIS36, AIS37, AIS38, AIS39, AIS40)),
                mean_TIA_case= mean(c(TIA1, TIA2, TIA3, TIA4, TIA5, TIA6, TIA7, TIA8, TIA9, TIA10, TIA11, TIA12, TIA13, TIA14, TIA15, TIA16,
                                      TIA17, TIA18, TIA19)),
                log_fc = mean_IS_case - mean_TIA_case,
                log_pval = -1*log10(p_val))
View(dat_fc)

#Save final data with list of final data in csv file
write.csv(dat_fc, "Final_IS_TIA_Bak.csv")

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
  select(Protein_ID, AIS1:AIS40, TIA1:TIA19, mean_IS_case, mean_TIA_case, log_fc, log_pval, p_val)
View(final_data)

#Save final data with list of significant proteins in csv file
write.csv(final_data, "Final_diffproteins_IS_TIA_Bak.csv")

##################################################################################################################################################

#Boruta feature selection-based random forest method

set.seed(555)
Boruta.IS.TIA = Boruta(Outcome~.-Protein_ID, data = IS_TIA_data[[2]], doTrace = 2, maxRuns=1000)
print(Boruta.IS.TIA)
#Take a call on tentative features
Boruta.IS.TIA.all = TentativeRoughFix(Boruta.IS.TIA)
print(Boruta.IS.TIA.all)
Sign.boruta= names(Boruta.IS.TIA.all$finalDecision[Boruta.IS.TIA.all$finalDecision %in% c("Confirmed")])
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
plot.Boruta.sel(Boruta.IS.TIA.all, cex.axis= 0.5, las= 2, xlab= "",
                pars = c("APOB",  "AZGP1", "LBP", "SELP", "C1QA", "ATRN" ));

#Get only selected features
getSelectedAttributes(Boruta.IS.TIA.all, withTentative = F)
boruta.df = attStats(Boruta.IS.TIA.all)
class(boruta.df)
print(boruta.df)

write.csv(boruta.df, "Boruta_results_IS_TIA_Bak.csv")

#############################################################################################################################################################

#Variance Partitioning Analysis (VPA)

#Load the expression dataset
Exp.data= read.csv("Expression_data_IS_TIA_Bak.csv", header = TRUE, row.names = 1)
View(Exp.data)

#Load the meta dataset
Meta.data= read.csv("Meta_data_IS_TIA_Bak.csv", header = TRUE, row.names = 1)
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
vp <- sortCols(varPartResid, last= c("Outcome", "Age", "Sex", "DM", "HTN", "CAD", "AFIB", "Smoking", "NIHSS0", "GCS", "Residuals"))
vp

#Violin plot of contribution of each variable to total variance
plotVarPart(vp, main="Violin Plot for total variance in proteins between AIS and TIA")

#Access first entries for outcome
head(vp$Outcome)

#Sort proteins based on variance explained by Outcome
head(vp[order(vp$Outcome, decreasing = TRUE), ])
vpa.outcome= vp[order(vp$Outcome, decreasing = TRUE), ]
vpa.outcome
write.csv(vpa.outcome, "Protein variance_Outcome_IS_TIA_Bak.csv")

#Sort proteins by variance explained by Outcome 
OutcomeSortOrder<-order(vp[["Outcome"]],decreasing=TRUE)  
for (i in ls(vp)) { vp[[i]]<-vp[[i]][OutcomeSortOrder]; }
rownames(vp)<-rownames(vp)[OutcomeSortOrder]
#Plot top proteins associated with outcome
plotPercentBars( vp[1:10,]) + ggtitle( "Top 10 Proteins associated with Outcome (IS vs. TIA)" )