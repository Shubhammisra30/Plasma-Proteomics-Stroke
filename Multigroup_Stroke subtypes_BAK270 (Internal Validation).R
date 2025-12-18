#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library(mixOmics)

#sPLS-DA Multigroup comparison between stroke subtypes (BAK270 Proteomics, Internal Validation)

Bak= read.csv("Log_Bak_data.csv", header = TRUE, row.names = 1)
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
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = "Stroke subtypes: PLS-DA plot with confidence ellipses")
# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(All.data, comp.predicted=2, dist = "centroids.dist") # Choose one of the three following distances: 'max.dist', 'centroids.dist' or 'mahalanobis.dist'
# plot the samples projected onto the first two components of the PLS-DA subspace
X11()
plotIndiv(All.data, comp = 1:2,
          group = Outcome$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: PLS-DA prediction background plot")

## TUNING sPLS-DA ##
# Selecting the number of components:
# undergo performance evaluation in order to tune the number of components to use
set.seed(222)
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
                               dist = 'max.dist', # Choose one of the three following distances: 'max.dist', 'centroids.dist' or 'mahalanobis.dist'
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
plotIndiv(final.splsda, comp = c(1,4), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          star=TRUE, title = 'Stroke subtypes: sPLS-DA plot with confidence ellipses (Comp 1,4)')
background = background.predict(final.splsda, comp.predicted=2, dist = "max.dist")
plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = Outcome$Outcome, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = "Stroke subtypes: sPLS-DA prediction background plot (Comp 1,2)")

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

# Cluster Image Map (CIM):
# set the styling of the legend to be homogeneous with previous plots
legend=list(legend = levels(Y), # set of classes
            col = unique(color.mixo(Y)), # set of colours
            title = "Stroke subtypes: Cluster Image Map", # legend title
            cex = 0.7) # legend size
# generate the CIM, using the legend and colouring rows by each sample's class
# with the elements (genes) for the optimized ncomponents (2 in this example) and with the  elements in each separated Component:
X11()
cim <- cim(final.splsda, row.sideColors = color.mixo(Y), cluster= "column", 
           legend = legend, margins= c(5,5))     # depicts the expression levels of each protein (selected for component construction) for every sample.

# Variable Plots:
# form new perf() object which utilises the final model
perf.splsda.data = perf(final.splsda, 
                        folds = 10, nrepeat = 50, # use repeated cross-validation
                        validation = "Mfold", dist = "max.dist",  # use max.dist measure
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

# How to extract the names of variables selected (for each component) when performing sPLS-DA?
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
# Just the head of the stability of the selected var
head(cbind(selectVar(final.splsda, comp=1)$value, stability.comp1))
head(cbind(selectVar(final.splsda, comp=2)$value, stability.comp2))
head(cbind(selectVar(final.splsda, comp=3)$value, stability.comp3))
head(cbind(selectVar(final.splsda, comp=4)$value, stability.comp4))

#Save Variables in Comp 1, Comp 2, and Comp 3 in a .csv file
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