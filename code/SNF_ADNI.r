### Running SNF on ADNI data
### Run this after determining k, alpha, and t using SNF_ADNI_params.r

#--- Include libraries ---#

if (!require("SNFtool")) install.packages("SNFtool")
library(SNFtool)

if (!require("Hmisc")) install.packages("Hmisc")
library(Hmisc)

if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

if (!require("magrittr")) install.packages("magrittr")
library(magrittr)

if (!require("cluster")) install.packages("cluster")
library(cluster)

if (!require("MASS")) install.packages("MASS")
library(MASS)

# Install ADNIMERGE (ADNI data)
library(ADNIMERGE)

# From Grace Jacobs: source code for function to determine data integration and clustering across resampling using SNF and spectral clustering
# https://github.com/gracejacobs/SNF-NDD-Clustering/blob/master/code/1_Running_SNF_analysis/robust_core-clustering_function.R
source("robust_core-clustering_function.R")

#--- Preparing data ---#

# Clinical Dementia Rating
d1 <- cdr
# Raw Alzheimer's Disease Assessment Scale-Cognitive data
d2 <- adas
# Montreal Cognitive Assessment
d3 <- moca
# Mini Mental State Exam
d4 <- mmse
# Cross-Sectional FreeSurfer (6.0)
d5 <- ucsffsx6

# Select columns of interest
ds1 <- subset(d1, select=c(RID, VISCODE, CDRSB))
ds1 <- ds1 %>% drop_na(CDRSB)
ds2 <- subset(d2, select=c(RID, VISCODE, TOTAL13))
ds2 <- ds2 %>% drop_na(TOTAL13)
ds3 <- subset(d3, select=c(RID, VISCODE, MOCA))
ds3 <- ds3 %>% drop_na(MOCA)
ds4 <- subset(d4, select=c(RID, VISCODE, MMSCORE))
ds4 <- ds4 %>% drop_na(MMSCORE)
ds5 <- d5 %>% dplyr::select(ends_with('CV') | ends_with('SA') | ends_with('TA'), RID, VISCODE)

# Merge and reduce dfs
mdf<- Reduce(function(x, y) merge(x, y, by = c("RID", "VISCODE"), all.x=TRUE), list(ds1, ds2, ds3, ds4, ds5))

# Remove all NA
nona_df <- na.omit(mdf)
# Remove duplicates by RID, keeping first value
nona_df <- nona_df[order(nona_df$VISCODE ,decreasing=TRUE),]
nona_df <- nona_df%>% distinct(RID, .keep_all = TRUE)
# Remove any participants that have NA for diagnosis
ddx <- dxsum %>% filter(COLPROT == "ADNI3"| COLPROT == "ADNI1")
ddx <- subset(ddx, select = c(DIAGNOSIS, RID, VISCODE))
ddx1<- Reduce(function(x, y) merge(x, y, by = c("RID", "VISCODE"), all.x=TRUE), list(ddx, nona_df))
ddx1 <- ddx1 %>% drop_na(DIAGNOSIS)
# Removing NAs
ddx1_na <- na.omit(ddx1)

# Participant RID and VISCODE for selection of demo info in setup_cluster_analysis
RIDVIS <- subset(ddx1_na, select=c(RID, VISCODE))
              
# Save RIDVIS
directory <- ("output/1_SNF_analysis/")
write.csv(RIDVIS, file=file.path(directory, paste("ridvis.csv", sep="")))

# Separate data
participant <- ddx1_na$RID
behavioural <- subset(ddx1_na, select=c(CDRSB, TOTAL13, MOCA, MMSCORE))
corticalVolume <- ddx1_na %>% dplyr::select(ends_with('CV'))
surfaceArea <- ddx1_na %>% dplyr::select(ends_with('SA'))
thicknessAvg <- ddx1_na %>% dplyr::select(ends_with('TA'))
# Format
corticalVolume <- as.data.frame(corticalVolume)
surfaceArea <- as.data.frame(surfaceArea)

#--- Option 1: Divide by ICV and remove ICV ---#
#Option 1 is used for analyses presented in paper - comparison of both options showed similar patterns.
              
icv <- corticalVolume$ST10CV
corticalVolume = subset(corticalVolume, select = -c(ST10CV))
CV_icv <- sweep(corticalVolume, 2, icv, FUN = '/')
SA_icv <- sweep(surfaceArea, 2, icv, FUN = '/')

# Write CSV
write.csv(participant, file=file.path(directory, paste("IDs.csv", sep="")))
write.csv(behavioural, file=file.path(directory, paste("Behavioural.csv", sep="")))
write.csv(CV_icv, file=file.path(directory, paste("CV.csv", sep="")))
write.csv(SA_icv, file=file.path(directory, paste("SA.csv", sep="")))
write.csv(thicknessAvg, file=file.path(directory, paste("TA.csv", sep="")))

#--- Option 2: Regressing out age and sex ---#

# Convert matrix to tibble and name the columns
matrix_data <- cbind(participant, behavioural)
matrix_data <- cbind(matrix_data, CV_icv)
matrix_data <- cbind(matrix_data, SA_icv)
matrix_data <- cbind(matrix_data, thicknessAvg)
matrix_df <- as_tibble(matrix_data)

# Create a tibble for age and sex
ids <-read.csv(file.path(directory, paste("IDs.csv", sep="")))
names(ids)[names(ids)=="x"] <- "RID"

# Demo info for those included in our analysis
demo <- merge(ptdemog, ids, by="RID")
# Remove duplicates by RID, keeping last value
demo <- demo[order(demo$VISCODE ,decreasing=FALSE),]
demo <- demo %>% distinct(RID, .keep_all = TRUE)
# Select age and sex (PTGENDER)
demo2 <- subset(demo, select=c(RID, AGE, PTGENDER))
names(demo2)[names(demo2) == 'RID'] <- 'participant'

allN <- merge(matrix_data, demo2, by="participant", all.x=TRUE)
covariates <- tibble(age = allN$AGE, sex = allN$PTGENDER)

#save RDA
save(covariates, file='SNF_covariates.rda')
save(matrix_df, file='SNF_matrix_df.rda')

# If there is a 'participant' column, use a left join or merge by the participant identifier. 
combined_df <- cbind(matrix_df, covariates) 
# Remove rows with NA values in age or sex (gender) 
combined_df <- combined_df %>% filter(!is.na(age) & !is.na(sex)) 
# Create a new matrix to store the residuals 
residual_matrix <- matrix(NA, nrow = nrow(combined_df), ncol = ncol(matrix_df) - 1) 
colnames(residual_matrix) <- colnames(matrix_df)[-1]  
# Exclude "participant" column (if it exists) 
rownames(residual_matrix) <- combined_df$participant 
# Loop through each column (except "participant") and regress it on age and sex (gender) 
for (i in 2:ncol(matrix_df)) {  
  # Start from 2 to skip the "participant" column if it exists   
  column_name <- colnames(matrix_df)[i]      
  # Create the formula for regression   
  formula <- as.formula(paste(column_name, "~ age + sex"))      
  # Fit the linear model   
  fit <- lm(formula, data = combined_df)      
  # Extract residuals and store in the new matrix   
  residual_matrix[, column_name] <- residuals(fit) } 
# Convert the residual matrix back to a data frame and add participant IDs as the first column 
residual_df <- data.frame(participant = combined_df$participant, residual_matrix) 
# Save the residuals to a new matrix format 
residual_matrix <- as.matrix(residual_df[-1])  


# End of option 2
#-------------------------------#

#--- Running SNF after params ---#

# Normalizing measures using function from SNF package for continuous variables but optional
thicknessAvg = standardNormalization(thicknessAvg)
corticalVolume = standardNormalization(CV_icv)
surfaceArea = standardNormalization(SA_icv)
behavioural = standardNormalization(behavioural)

# Setting the parameters (finalized after comparisons using SNF_ADNI_params.r )
K =18;		# number of neighbors, usually (10~30), usually sample size/10
alpha = 0.8;  	# hyperparameter, usually (0.3~0.8)
t = 10; 	# Number of Iterations, usually (10~20) 

# Creating participant distance matrices using euclidean distances
Dist_behavioural = dist2(as.matrix(behavioural),as.matrix(behavioural));
Dist_volume = dist2(as.matrix(corticalVolume),as.matrix(corticalVolume));
Dist_area = dist2(as.matrix(surfaceArea),as.matrix(surfaceArea));
Dist_thickness = dist2(as.matrix(thicknessAvg),as.matrix(thicknessAvg));

# Creating participant affinity matrices within each data type
AM_behavioural = affinityMatrix(Dist_behavioural,K, alpha)
AM_volume = affinityMatrix(Dist_volume, K, alpha)
AM_area = affinityMatrix(Dist_area, K, alpha)
AM_thickness = affinityMatrix(Dist_thickness, K, alpha)

# Setting output directory
directory <- ("C:\\Users\\bourq\\Documents\\0-Carleton\\6 - PhD\\Y2\\0-MR(John)\\SNF\\SNFcode")

# Setting cluster number
C = 4 
# calling function to integrate data types using SNF and cluster participants using spectral clustering across resampling 80% of participants 1000 times
robust.W = RobustCoreClusteringMatrix(feature.affinity.mat.list = list(AM_behavioural, AM_volume, AM_area, AM_thickness),
                                      exp.num.samples = 1000, num.clusts = C)
# Two matrices - Dense Core Cluster Matrix and Sparse Core Cluster Matrix
dense <- robust.W[1]
dense <- matrix(unlist(dense), ncol = 515, byrow = TRUE)
sparse <- robust.W[2]
sparse <- matrix(unlist(sparse), ncol = 515, byrow = TRUE)

# Displaying clusters
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))

# Saving an image of the cluster heatmap
png('SN_dense_0.8_18_1000perms.png')
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
dev.off()

# Calculating normalized mutual information (NMI) based off of the original data types and the clustering similarity matrix
SNF_NMIScores <-rankFeaturesByNMI(list(thicknessAvg, corticalVolume, surfaceArea, behavioural), dense) 
SNF_NMIScores 

# Separating and organizing scores
TA_scores <- as.data.frame(SNF_NMIScores[[1]][1])
CV_scores <- as.data.frame(SNF_NMIScores[[1]][2])
SA_scores <- as.data.frame(SNF_NMIScores[[1]][3])
behav_scores <- as.data.frame(SNF_NMIScores[[1]][4])
names(TA_scores) <- c("NMI")
names(CV_scores) <- c("NMI")
names(SA_scores) <- c("NMI")
names(behav_scores) <- c("NMI")

# Saving csv of scores
directory <- ("output/1_SNF_analysis")
write.csv(TA_scores, file=file.path(directory, paste("TA_scores_k18_0.8_1000perms.csv", sep="")))
write.csv(CV_scores, file=file.path(directory, paste("CV_scores_k18_0.8_1000perms.csv", sep="")))
write.csv(SA_scores, file=file.path(directory, paste("SA_scores_k18_0.8_1000perms.csv", sep="")))
write.csv(behav_scores, file=file.path(directory, paste("behav_scores_k18_0.8_1000perms.csv", sep="")))

# Saving csv of NMI scores and clustering similarity matrix
all_scores <- rbind(TA_scores, CV_scores)
all_scores <- rbind(all_scores, SA_scores)
all_scores <- rbind(all_scores, behav_scores)
write.csv(all_scores, file=file.path(directory, paste("all_scores_k18_0.8_1000perms.csv", sep="")))
write.matrix(dense, file=file.path(directory, paste("dense_similarity_matrix.csv", sep="")))

# Find cluster labels of individuals using the robust clustering similarity matrix
robust.groups.df = RobustCoreClusteringClusters(core.clustering.list = robust.W,num.clusts = C,verbose = T)
clusters <- cbind(participant, robust.groups.df)
table(clusters$groups)

#--- Silhouette plot ---#

# Calculating silhouette width for each participant and the silhouette plot
dissim <- 1 - dense
dissim <- as.matrix(dissim)

# Format
clusters$groups <- as.integer(clusters$groups)

# Switching group names for better ordering
newClusters <- clusters
newClusters$groups[newClusters$groups == '3'] <- '5'
newClusters$groups[newClusters$groups == '4'] <- '3'
newClusters$groups[newClusters$groups == '5'] <- '4'
newClusters$groups <- as.integer(newClusters$groups)

sil <- silhouette(newClusters$groups, dmatrix = dissim)

# Adding silhouette widths for each individual
newClusters$silhouette_width <- 0
for (i in 1:515){
  clusters[i, 4] <- sil[i ,3]
}

write.csv(newClusters, file=file.path(directory, paste("4clust_groups_k18_0.8_1000perms.csv", sep="")))

# Saving the silhouette plot
name=("Silhouette_plot_0.8_18.png")
png(name)
plot(sil, col = c("#bc3c29", "#e18727", "#20854e", "#0072b5"), border = NA)
dev.off()

# Check if anyone doesn't reliably cluster 
unreliable.patients = UnrelClustPatientsFullCohort(core.clustering.list = robust.W,verbose = TRUE)

# Check if anyone doesn't reliably cluster within their given cluster group
unreliable.patients.by.grp = UnrelClustPatientsByGrp(core.clustering.list = robust.W,id.group.df = robust.groups.df,verbose = TRUE)




