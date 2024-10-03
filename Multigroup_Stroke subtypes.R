#Multigroup comparisons of Stroke subtypes using sPLS-DA approach

#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library (readxl)
library(writexl)
library(gplots)
library(ggrepel)
library(randomForest)
library(pROC)
library(caret)
library(rsample)
library(mixOmics)

#Comparisons between AIS, ICH, TIA, and MIM (SomaScan Proteomics)

#Proteomics data file upload
excel_sheets("All_groups.xlsx")
All_data= excel_sheets("All_groups.xlsx") %>% map(~read_xlsx("All_groups.xlsx",.))
All_data

#Assigning variables
set.seed(123) # for reproducibility, remove for normal use

X= All_data[[4]][,3:7309]
Y= All_data[[4]]$Outcome
View(X)
View(Y)
View(All_data[[4]])
dim(X); length(Y)
Y=as.factor(Y)
summary(Y)

## INITIAL ANALYSIS ##
# Preliminary (unsupervised) Analysis with PCA
pca.data = pca(X, ncomp = 10, center = TRUE, scale = TRUE) # run pca method on data
plot(pca.data)  # barplot of the eigenvalues (explained variance per component)
plotIndiv(pca.data, group = All_data[[4]]$Outcome, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'Stroke subtypes: PCA plot') # onto the PCA subspace

## INITIAL PLS-DA MODEL ##
All.data = splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data , comp = 1:2, 
          group = All_data[[4]]$Outcome, ind.names = FALSE,  # colour points by class
          ellipse = FALSE,
          legend = TRUE, title = "Stroke subtypes: PLS-DA plot")
# use the centroids.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(All.data, comp.predicted=2, dist = "centroids.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data, comp = 1:2,
          group = All_data[[4]]$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: PLS-DA prediction background plot")

## TUNING sPLS-DA ##
# Selecting the number of components:
# undergo performance evaluation in order to tune the number of components to use
set.seed(111)
perf.splsda.data = perf(All.data, validation = "Mfold", 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        progressBar = TRUE, auc = TRUE, cpus = 4) # include AUC values
# plot the outcome of performance evaluation across all ten components
X11()
plot(perf.splsda.data, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.data$choice.ncomp # what is the optimal value of components according to perf()
perf.splsda.data$auc
perf.splsda.data$error.rate
perf.splsda.data$error.rate.class

# Selecting the number of variables:
# grid of possible keepX values that will be tested for each component
set.seed(123)
list.keepX = c(1:9,  seq(10, 30, 10))
list.keepX
# undergo the tuning process to determine the optimal number of variables
tune.splsda.data = tune.splsda(X, Y, ncomp = 5, # calculate for first 5 components
                               validation = 'Mfold',
                               folds = 10, nrepeat = 50, # use repeated cross-validation
                               dist = 'centroids.dist',
                               measure = "BER", # use balanced error rate of dist measure; 'measure' must be either 'BER', 'overall' or 'AUC'
                               test.keepX = list.keepX, progressBar = TRUE,
                               cpus = 4) # allow for parallelisation to decrease runtime

plot(tune.splsda.data, col = color.jet(5)) # plot output of variable number tuning
tune.splsda.data$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
tune.splsda.data$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp = tune.splsda.data$choice.ncomp$ncomp
optimal.ncomp
optimal.keepX = tune.splsda.data$choice.keepX[1:optimal.ncomp]
optimal.keepX

## FINAL MODEL
# form final model with optimised values for component and variable count
final.splsda = splsda(X, Y, 
                      ncomp = optimal.ncomp, 
                      keepX = optimal.keepX)

# Sample Visualization
# Sample Plots
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = All_data[[4]]$Outcome, ind.names = FALSE, # colour by class label
          ellipse = FALSE, legend = TRUE,
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot (Comp 1,2)')
plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = All_data[[4]]$Outcome, ind.names = FALSE, # colour by class label
          ellipse = FALSE, legend = TRUE,
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot (Comp 1,3)')

background = background.predict(final.splsda, comp.predicted=2, dist = "centroids.dist")

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = All_data[[4]]$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,2)")
plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = All_data[[4]]$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,3)")

# Variable Visualization
# Plot Loadings
X11()
plotLoadings(final.splsda, comp = 1, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)
X11()
plotLoadings(final.splsda, comp = 2, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)
X11()
plotLoadings(final.splsda, comp = 3, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)

# Variable Plots:
# form new perf() object which utilises the final model
perf.splsda.data = perf(final.splsda, 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        validation = "Mfold", dist = "centroids.dist",
                        progressBar = TRUE)
perf.splsda.data$error.rate.class

# plot the stability of each feature for the first two components, 'h' type refers to histogram
par(mfrow=c(1,2)) #Sets the layout of plots to be a single row and two columns
#For Component 1
stable.comp1= perf.splsda.data$features$stable$comp1
barplot(stable.comp1, xlab= "Features across CV folds", ylab= "Stability Frequency",
        main = "Feature stability: Comp 1", las=1)
#For Component 2
stable.comp2= perf.splsda.data$features$stable$comp2
barplot(stable.comp2, xlab= "Features across CV folds", ylab= "Stability Frequency",
        main = "Feature stability: Comp 2", las=1)
par(mfrow=c(1,1)) #Sets the layout of plots to be a single row and a single column

# Extract the names of variables selected (for each component) when performing sPLS-DA
# First extract the name of selected var
select.name.comp1= selectVar(final.splsda, comp = 1)$name
select.name.comp2= selectVar(final.splsda, comp = 2)$name
select.name.comp3= selectVar(final.splsda, comp = 3)$name
# Then extract the stability values from perf:
stability.comp1= perf.splsda.data$features$stable$comp1[select.name.comp1]
stability.comp2= perf.splsda.data$features$stable$comp1[select.name.comp2]
stability.comp3= perf.splsda.data$features$stable$comp1[select.name.comp3]
# Display just the head of the stability of the selected var
head(cbind(selectVar(final.splsda, comp=1)$value, stability.comp1))
head(cbind(selectVar(final.splsda, comp=2)$value, stability.comp2))
head(cbind(selectVar(final.splsda, comp=3)$value, stability.comp3))

#Save Variables in Comp 1, Comp 2, and Comp 3 in a .csv file
Proteins_comp1= cbind(selectVar(final.splsda, comp=1)$value, stability.comp1)
write.csv(Proteins_comp1, "Proteins in Comp 1_All groups.csv")
Proteins_comp2= cbind(selectVar(final.splsda, comp=2)$value, stability.comp2)
write.csv(Proteins_comp2, "Proteins in Comp 2_All groups.csv")
Proteins_comp3= cbind(selectVar(final.splsda, comp=3)$value, stability.comp3)
write.csv(Proteins_comp3, "Proteins in Comp 3_All groups.csv")

# correlation circle plot:
plotVar(final.splsda, comp = c(1,2), cex = 3) # generate correlation circle plot

## PERFORMANCE PLOTS (ROC)
auc.splsda = auroc(final.splsda, roc.comp = 1, print = TRUE) # AUROC for the first component
auc.splsda = auroc(final.splsda, roc.comp = 2, print = TRUE) # AUROC for all two components
auc.splsda = auroc(final.splsda, roc.comp = 3, print = TRUE) # AUROC for all three components

#########################################################################################################################################################
#########################################################################################################################################################

#Sensitivity Analysis- Multigroup comparisons between three stroke subtypes (AIS, TIA, and MIM) excluding ICH
#SomaScan Proteomics

#Proteomics data file upload
excel_sheets("Three_groups.xlsx")
All_data= excel_sheets("Three_groups.xlsx") %>% map(~read_xlsx("Three_groups.xlsx",.))
All_data

#Assigning variables
set.seed(123) # for reproducibility, remove for normal use

X= All_data[[2]][,3:7309]
Y= All_data[[2]]$Outcome
View(X)
View(Y)
View(All_data[[2]])
dim(X); length(Y)
Y=as.factor(Y)
summary(Y)

## INITIAL ANALYSIS ##
# Preliminary (unsupervised) Analysis with PCA
pca.data = pca(X, ncomp = 10, center = TRUE, scale = TRUE) # run pca method on data
plot(pca.data)  # barplot of the eigenvalues (explained variance per component)
plotIndiv(pca.data, group = All_data[[2]]$Outcome, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'Stroke subtypes: PCA plot') # onto the PCA subspace

## INITIAL PLS-DA MODEL ##
All.data = splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data , comp = 1:2, 
          group = All_data[[2]]$Outcome, ind.names = FALSE,  # colour points by class
          ellipse = FALSE,
          legend = TRUE, X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2', title = "Stroke subtypes: PLS-DA plot")
# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(All.data, comp.predicted=2, dist = "max.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data, comp = 1:2,
          group = All_data[[2]]$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2', title = "Stroke subtypes: PLS-DA prediction background plot")

## TUNING sPLS-DA ##
# Selecting the number of components:
# undergo performance evaluation in order to tune the number of components to use
set.seed(111)
perf.splsda.data = perf(All.data, validation = "Mfold", 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        progressBar = TRUE, auc = TRUE, cpus = 4) # include AUC values
# plot the outcome of performance evaluation across all ten components
X11()
plot(perf.splsda.data, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.data$choice.ncomp # what is the optimal value of components according to perf()
perf.splsda.data$auc
perf.splsda.data$error.rate
perf.splsda.data$error.rate.class

# Selecting the number of variables:
# grid of possible keepX values that will be tested for each component
set.seed(123)
list.keepX = c(1:9,  seq(10, 20, 5))
list.keepX
# undergo the tuning process to determine the optimal number of variables
tune.splsda.data = tune.splsda(X, Y, ncomp = 5, # calculate for first 5 components
                               validation = 'Mfold',
                               folds = 10, nrepeat = 10, # use repeated cross-validation
                               dist = 'max.dist',
                               measure = "BER", # use balanced error rate of dist measure; 'measure' must be either 'BER', 'overall' or 'AUC'
                               test.keepX = list.keepX, progressBar = TRUE,
                               cpus = 4) # allow for parallelisation to decrease runtime

plot(tune.splsda.data, col = color.jet(5)) # plot output of variable number tuning
tune.splsda.data$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
tune.splsda.data$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp = tune.splsda.data$choice.ncomp$ncomp
optimal.ncomp
optimal.keepX = tune.splsda.data$choice.keepX[1:optimal.ncomp]
optimal.keepX

## FINAL MODEL
# form final model with optimised values for component and variable count
final.splsda = splsda(X, Y, 
                      ncomp = optimal.ncomp, 
                      keepX = optimal.keepX)

# Sample Visualization
# Sample Plots
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = All_data[[2]]$Outcome, ind.names = FALSE, # colour by class label
          ellipse = FALSE, legend = TRUE,
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot (Comp 1,2)')
background = background.predict(final.splsda, comp.predicted=2, dist = "max.dist")
X11()
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = All_data[[2]]$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,2)")

# Variable Visualization
# Plot Loadings
X11()
plotLoadings(final.splsda, comp = 1, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)
X11()
plotLoadings(final.splsda, comp = 2, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)

# Variable Plots:
# form new perf() object which utilises the final model
perf.splsda.data = perf(final.splsda, 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        validation = "Mfold", dist = "max.dist",
                        progressBar = TRUE)
perf.splsda.data$error.rate.class

# plot the stability of each feature for the first two components, 'h' type refers to histogram
par(mfrow=c(1,2)) #Sets the layout of plots to be a single row and two columns
#For Component 1
stable.comp1= perf.splsda.data$features$stable$comp1
barplot(stable.comp1, xlab= "Features across CV folds", ylab= "Stability Frequency",
        main = "Feature stability: Comp 1", las=1)
#For Component 2
stable.comp2= perf.splsda.data$features$stable$comp2
barplot(stable.comp2, xlab= "Features across CV folds", ylab= "Stability Frequency",
        main = "Feature stability: Comp 2", las=1)
par(mfrow=c(1,1)) #Sets the layout of plots to be a single row and a single column

# Extract the names of variables selected (for each component) when performing sPLS-DA
# First extract the name of selected var
select.name.comp1= selectVar(final.splsda, comp = 1)$name
select.name.comp2= selectVar(final.splsda, comp = 2)$name
# Then extract the stability values from perf:
stability.comp1= perf.splsda.data$features$stable$comp1[select.name.comp1]
stability.comp2= perf.splsda.data$features$stable$comp1[select.name.comp2]
# Just the head of the stability of the selected var
head(cbind(selectVar(final.splsda, comp=1)$value, stability.comp1))
head(cbind(selectVar(final.splsda, comp=2)$value, stability.comp2))

#Save Variables in Comp 1 and Comp 2 in a .csv file
Proteins_comp1= cbind(selectVar(final.splsda, comp=1)$value, stability.comp1)
write.csv(Proteins_comp1, "Proteins in Comp 1_Three groups.csv")
Proteins_comp2= cbind(selectVar(final.splsda, comp=2)$value, stability.comp2)
write.csv(Proteins_comp2, "Proteins in Comp 2_Three groups.csv")

# correlation circle plot:
plotVar(final.splsda, comp = c(1,2), cex = 3) # generate correlation circle plot

## PERFORMANCE PLOTS (ROC)
auc.splsda = auroc(final.splsda, roc.comp = 1, print = TRUE) # AUROC for the first component
auc.splsda = auroc(final.splsda, roc.comp = 2, print = TRUE) # AUROC for all two components

#########################################################################################################################################################
#########################################################################################################################################################

#Cross-platform Validation
#Multigroup comparisons between AIS, ICH, TIA, and MIM (BAK-270 Proteomics)

#File upload for BAK-270 validation

Bak= read.csv("Bak_Data.csv", header = TRUE, row.names = 1)
View(Bak)
Outcome= read.csv("Outcome_Data.csv", header = TRUE, row.names = 1)
View(Outcome)

set.seed(123) # for reproducibility, remove for normal use

X= Bak
Y= Outcome$Outcome

View(X)
View(Y)
View(Outcome$Outcome)
dim(X); length(Y)
Y=as.factor(Y)
summary(Y)

## INITIAL ANALYSIS ##
# Preliminary (unsupervised) Analysis with PCA
pca.data = pca(X, ncomp = 10, center = TRUE, scale = TRUE) # run pca method on data
plot(pca.data)  # barplot of the eigenvalues (explained variance per component)
plotIndiv(pca.data, group = Outcome$Outcome, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'Stroke subtypes: PCA plot') # onto the PCA subspace

## INITIAL PLS-DA MODEL ##
All.data = splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data , comp = 1:2, 
          group = Outcome$Outcome, ind.names = FALSE,  # colour points by class
          ellipse = FALSE,
          legend = TRUE, title = "Stroke subtypes: PLS-DA plot")
# use the centroids.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(All.data, comp.predicted=2, dist = "centroids.dist")
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data, comp = 1:2,
          group = Outcome$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: PLS-DA prediction background plot")

## TUNING sPLS-DA ##
# Selecting the number of components:
# undergo performance evaluation in order to tune the number of components to use
set.seed(111)
perf.splsda.data = perf(All.data, validation = "Mfold", 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        progressBar = TRUE, auc = TRUE, cpus = 4) # include AUC values
# plot the outcome of performance evaluation across all ten components
X11()
plot(perf.splsda.data, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.data$choice.ncomp # what is the optimal value of components according to perf()
perf.splsda.data$auc
perf.splsda.data$error.rate
perf.splsda.data$error.rate.class

# Selecting the number of variables:
# grid of possible keepX values that will be tested for each component
set.seed(123)
list.keepX = c(10:14,  seq(15, 20, 5))
list.keepX
# undergo the tuning process to determine the optimal number of variables
tune.splsda.data = tune.splsda(X, Y, ncomp = 5, # calculate for first 5 components
                               validation = 'Mfold',
                               folds = 10, nrepeat = 50, # use repeated cross-validation
                               dist = 'centroids.dist',
                               measure = "BER", # use balanced error rate of dist measure; 'measure' must be either 'BER', 'overall' or 'AUC'
                               test.keepX = list.keepX, progressBar = TRUE,
                               cpus = 4) # allow for parallelisation to decrease runtime

plot(tune.splsda.data, col = color.jet(5)) # plot output of variable number tuning
tune.splsda.data$choice.ncomp$ncomp # what is the optimal value of components according to tune.splsda()
tune.splsda.data$choice.keepX # what are the optimal values of variables according to tune.splsda()

optimal.ncomp = tune.splsda.data$choice.ncomp$ncomp
optimal.ncomp
optimal.keepX = tune.splsda.data$choice.keepX[1:optimal.ncomp]
optimal.keepX

## FINAL MODEL
# form final model with optimised values for component and variable count
final.splsda = splsda(X, Y, 
                      ncomp = optimal.ncomp, 
                      keepX = optimal.keepX)

# Sample Visualization
# Sample Plots
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour by class label
          ellipse = FALSE, legend = TRUE,
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot (Comp 1,2)')
plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour by class label
          ellipse = FALSE, legend = TRUE,
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot (Comp 1,3)')
plotIndiv(final.splsda, comp = c(1,4), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour by class label
          ellipse = FALSE, legend = TRUE,
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot (Comp 1,4)')

background = background.predict(final.splsda, comp.predicted=2, dist = "centroids.dist")

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,2)")
plotIndiv(final.splsda, comp = c(1,3), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,3)")
plotIndiv(final.splsda, comp = c(1,4), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,4)")

# Variable Visualization
# Plot Loadings
X11()
plotLoadings(final.splsda, comp = 1, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)
X11()
plotLoadings(final.splsda, comp = 2, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)
X11()
plotLoadings(final.splsda, comp = 3, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)
X11()
plotLoadings(final.splsda, comp = 4, method = "mean", contrib = "max", size.name = 0.6, size.legend = 0.9, legend = TRUE, ndisplay = 10)

# Variable Plots:
# form new perf() object which utilises the final model
perf.splsda.data = perf(final.splsda, 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        validation = "Mfold", dist = "centroids.dist",
                        progressBar = TRUE)
perf.splsda.data$error.rate.class

# plot the stability of each feature for the first two components, 'h' type refers to histogram
par(mfrow=c(1,2)) #Sets the layout of plots to be a single row and two columns
#For Component 1
stable.comp1= perf.splsda.data$features$stable$comp1
barplot(stable.comp1, xlab= "Features across CV folds", ylab= "Stability Frequency",
        main = "Feature stability: Comp 1", las=1)
#For Component 2
stable.comp2= perf.splsda.data$features$stable$comp2
barplot(stable.comp2, xlab= "Features across CV folds", ylab= "Stability Frequency",
        main = "Feature stability: Comp 2", las=1)
par(mfrow=c(1,1)) #Sets the layout of plots to be a single row and a single column

# Extract the names of variables selected (for each component) when performing sPLS-DA
# First extract the name of selected var
select.name.comp1= selectVar(final.splsda, comp = 1)$name
select.name.comp2= selectVar(final.splsda, comp = 2)$name
select.name.comp3= selectVar(final.splsda, comp = 3)$name
select.name.comp4= selectVar(final.splsda, comp = 4)$name
# Then extract the stability values from perf:
stability.comp1= perf.splsda.data$features$stable$comp1[select.name.comp1]
stability.comp2= perf.splsda.data$features$stable$comp1[select.name.comp2]
stability.comp3= perf.splsda.data$features$stable$comp1[select.name.comp3]
stability.comp4= perf.splsda.data$features$stable$comp1[select.name.comp4]
# Display just the head of the stability of the selected var
head(cbind(selectVar(final.splsda, comp=1)$value, stability.comp1))
head(cbind(selectVar(final.splsda, comp=2)$value, stability.comp2))
head(cbind(selectVar(final.splsda, comp=3)$value, stability.comp3))
head(cbind(selectVar(final.splsda, comp=4)$value, stability.comp4))

#Save Variables in Comp 1, Comp 2, Comp 3, and Comp 4 in a .csv file
Proteins_comp1= cbind(selectVar(final.splsda, comp=1)$value, stability.comp1)
write.csv(Proteins_comp1, "Proteins in Comp 1_All groups.csv")
Proteins_comp2= cbind(selectVar(final.splsda, comp=2)$value, stability.comp2)
write.csv(Proteins_comp2, "Proteins in Comp 2_All groups.csv")
Proteins_comp3= cbind(selectVar(final.splsda, comp=3)$value, stability.comp3)
write.csv(Proteins_comp3, "Proteins in Comp 3_All groups.csv")
Proteins_comp4= cbind(selectVar(final.splsda, comp=4)$value, stability.comp4)
write.csv(Proteins_comp4, "Proteins in Comp 4_All groups.csv")

# correlation circle plot:
plotVar(final.splsda, comp = c(1,2), cex = 3) # generate correlation circle plot

## PERFORMANCE PLOTS (ROC)
auc.splsda = auroc(final.splsda, roc.comp = 1, print = TRUE) # AUROC for the first component
auc.splsda = auroc(final.splsda, roc.comp = 2, print = TRUE) # AUROC for all two components
auc.splsda = auroc(final.splsda, roc.comp = 3, print = TRUE) # AUROC for all three components
auc.splsda = auroc(final.splsda, roc.comp = 4, print = TRUE) # AUROC for all four components