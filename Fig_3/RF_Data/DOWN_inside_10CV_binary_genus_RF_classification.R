## Binary classification using DOWN sampling within 10-fold CV on AGP or HMP data
## Using Genera as features.

##Set working directory
setwd("PATH/TO/WORKING/DIR/")

library("randomForest")
library("plyr")
#install.packages("rfUtilities")
library("rfUtilities") 
#install.packages("caret")
library("caret") 
library(vegan)
library(pROC)

##Read the genus table for AGP or HMP
## EDIT file name - using taxa table produced from OTU table rarefied at depth = 10000
otu <- read.table("AGP_genus_rarefied.txt",sep="\t", row.names=1, comment.char="", header=TRUE, check.names = FALSE, stringsAsFactors = F)
dim(otu)
otu <- t(otu)
dim(otu)

# Read the mapping file
map <- read.table("subset_map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")
dim(map)

## Note, for HMP data, replace all occurrences of 'race' with 'ethnicity' in this code 
table(map$race) #for HMP, use map$ethnicity

## Not needed for HMP data
## For AGP data, let's relabel everything as "black", "hispanic", "asian" and "white"
#the assignment function levels can be used to change the levels
#the order will remain the same.
levels(map$race) <- c('black', 'asian','white','hispanic')
table(map$race)

##### Preprocessing #####
# Random Forest doesn't is sensitive to noise. 
# Hence, remove rare OTUs

#Filter bugs that don't appear in enough number of samples
## Only keep bugs that appear in greater than prevalence size (prevalence-size = smallest-group/2 = 5)
select <- colSums(otu > 0) > min(table(map$race))/2
otu.rare.removed <- otu[,select]
dim(otu.rare.removed)

#Normalize to relative abundance
otu.orig <- otu
otu <- otu.rare.removed
otu.t <- t(otu)
otu.t.norm <- sweep(otu.t, 2, colSums(otu.t) , '/')

##################### Transformation ################################
### possible transformation
# certain transformations that help make features more comparable across samples 
# can be useful when using relative abundance data.
# arcsine sqrt transformation is frequently used for transforming microbiome data.
otu.t.norm.asin <- asin(sqrt(otu.t.norm))
otu <- otu.t.norm.asin
#################################################################################

##### Prepare data-frame for running RF models ##########
# First transform the otu table to have samples as rows and features as columns
otu <- t(otu)
dim(otu)

#Make sure the mapping file and otu table have samples in the same order
common.ids <- intersect(rownames(map),rownames(otu))
length(common.ids)

otu <- otu[common.ids,]
dim(otu)

map <- map[common.ids,,drop=F] # drop=F from preventing R from converting it to a vector when there is only one column
dim(map)

#merge them into one data frame
data <- data.frame(otu, check.names = FALSE)
data$race = map$race
dim(data)

############ End of preprocessing ###################
################ Binarization of the multi-class###############
##For HMP, replace data$race with data$ethnicity
data$black_vs_all <- ifelse(data$race == "black", "black", "not_black") 
data$black_vs_all <- as.factor(data$black_vs_all)
data$asian_vs_all <- ifelse(data$race == "asian", "asian", "not_asian")
data$asian_vs_all <- as.factor(data$asian_vs_all)
data$hispanic_vs_all <- ifelse(data$race == "hispanic", "hispanic", "not_hispanic")
data$hispanic_vs_all <- as.factor(data$hispanic_vs_all)
data$white_vs_all <- ifelse(data$race == "white", "white", "not_white")
data$white_vs_all <- as.factor(data$white_vs_all)

dim(data)
## Remember that last 5 columns of data are race-related. 

summary(data$race)
summary(data$black_vs_all)
summary(data$asian_vs_all)
summary(data$hispanic_vs_all)
summary(data$white_vs_all)
levels(data$white_vs_all)
# [1] "not_white" "white" 
## This is default factor levels. Let's make the reference level "white"
data$white_vs_all <- relevel(data$white_vs_all, ref="white")
levels(data$white_vs_all) #[1] "white"     "not_white"
summary(data$white_vs_all)

####################### Binary classification using DOWN subsampling inside 10-fold CV ###################

##Set trainControl
ctrl <- trainControl(method = "cv",
                     number = 10,
                     summaryFunction=twoClassSummary,#twoClassSummary computes sensitivity, specificity and the area under the ROC curve
                     savePredictions = "final", 
                     classProbs = TRUE,
                     sampling = "down"
                     )

## Now run each OVA binary RF classifier. 

## ova = black_vs_all
OVA <- data$black_vs_all

## 10-fold CV using DOWN sampling

set.seed(42)
down_inside_B <- train( x = data[,1:(ncol(data)-5)],#features,
                         y = OVA,#response
                         method = "rf",
                         metric="ROC",
                         trControl = ctrl)

#Get the confusion matrix for the heldout (test) samples during CV.
cm_train_B <- confusionMatrix.train(down_inside_B,norm="none")

##Compute classification metrics
acc_train_B <- accuracy(cm_train_B$table)
metrics_B <- c(Sens = acc_train_B$sensitivity, 
               Spec = acc_train_B$specificity, 
               Acc = acc_train_B$PCC/100,
               prec = precision(cm_train_B$table),
               auroc = auc(down_inside_B$pred$obs,
                           down_inside_B$pred$black)[1],
               f1 = F_meas(cm_train_B$table))

######## Repeat for asian ######
## ova = asian_vs_all
OVA <- data$asian_vs_all

## 10-fold CV using DOWN sampling

set.seed(42)
down_inside_A <- train(x = data[,1:(ncol(data)-5)],#features,
                        y = OVA,#response
                        method = "rf",
                        metric="ROC",
                        trControl = ctrl)

#Get the confusion matrix for the heldout (test) samples during CV.
cm_train_A <- confusionMatrix.train(down_inside_A,norm="none")

##Compute classification metrics
acc_train_A <- accuracy(cm_train_A$table)
metrics_A <- c(Sens = acc_train_A$sensitivity, 
               Spec = acc_train_A$specificity, 
               Acc = acc_train_A$PCC/100,
               prec = precision(cm_train_A$table),
               auroc = auc(down_inside_A$pred$obs,
                           down_inside_A$pred$asian)[1],
               f1 = F_meas(cm_train_A$table))

######## Repeat for hispanic ######
OVA <- data$hispanic_vs_all

## 10-fold CV using DOWN sampling
set.seed(42)
down_inside_H <- train(x = data[,1:(ncol(data)-5)],#features,
                        y = OVA,#response
                        method = "rf",
                        metric="ROC",
                        trControl = ctrl)

#Get the confusion matrix for the heldout (test) samples during CV.
cm_train_H <- confusionMatrix.train(down_inside_H,norm="none")

##Compute classification metrics
acc_train_H <- accuracy(cm_train_H$table)
metrics_H <- c(Sens = acc_train_H$sensitivity, 
               Spec = acc_train_H$specificity, 
               Acc = acc_train_H$PCC/100,
               prec = precision(cm_train_H$table),
               auroc = auc(down_inside_H$pred$obs,
                           down_inside_H$pred$hispanic)[1],
               f1 = F_meas(cm_train_H$table))

######## Repeat for white ######

OVA <- data$white_vs_all

## 10-fold CV using DOWN sampling
set.seed(42)
down_inside_W <- train(x = data[,1:(ncol(data)-5)],#features,
                        y = OVA,#response
                        method = "rf",
                        metric="ROC",
                        trControl = ctrl)

#Get the confusion matrix for the heldout (test) samples during CV.
cm_train_W <- confusionMatrix.train(down_inside_W,norm="none")

##Compute classification metrics
acc_train_W <- accuracy(cm_train_W$table)
metrics_W <- c(Sens = acc_train_W$sensitivity, 
               Spec = acc_train_W$specificity, 
               Acc = acc_train_W$PCC/100,
               prec = precision(cm_train_W$table),
               auroc = auc(down_inside_W$pred$obs,
                           down_inside_W$pred$white)[1],
               f1 = F_meas(cm_train_W$table))


## Combine all metrics into one dataframe
metrics_df <- data.frame(metrics_A,metrics_B, metrics_H, metrics_W)
metrics_df <- as.data.frame(t(metrics_df))
rownames(metrics_df) <- c("asian_vs_all","black_vs_all","hispanic_vs_all","white_vs_all")
write.table(metrics_df,file="Results/DOWN-metrics.txt", sep="\t", col.names = NA)

############### Plot ROC curve for the k-fold CV on full data #########
pdf("Plots/Rarefy_DOWN_ROC_ethnicity_inside_genus.pdf")
par(col.lab = "white")
plot.roc(down_inside_A$pred$obs,
         down_inside_A$pred$asian, col="#3490DE",legacy.axes = T)
plot.roc(down_inside_B$pred$obs,
         down_inside_B$pred$black, col="#FFDE7D", add=T,legacy.axes = T)
plot.roc(down_inside_H$pred$obs,
         down_inside_H$pred$hispanic, col="#62D2A2", add=T,legacy.axes = T)
plot.roc(down_inside_W$pred$obs,
         down_inside_W$pred$white, col="#F73859", add=T,legacy.axes = T)
legend("bottomright",legend = c("asian","black","hispanic","white"),col = c("#3490DE","#FFDE7D","#62D2A2","#F73859"), lty=c(1,1,1), box.lty=0, cex=0.8)
par(col.lab = "black")
title(xlab = "False positive rate", ylab = "True positive rate")
dev.off()
