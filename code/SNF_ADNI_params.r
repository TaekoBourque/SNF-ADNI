#--- Code to include libraries ---#

if (!require("SNFtool")) install.packages("SNFtool")
library(SNFtool)

if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# Install ADNIMERGE
library(ADNIMERGE)

### Adapted from Grace Jacobs
#--- Checking iterations of parameters and their effect on cluster number and NMI scores ---#

# Creating output folders
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.3")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.4")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.5")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.6")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.7")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.8")

# Setting up empty dataframes to compare nearest neighbors (K) and alpha parameters
ktests=seq(10,30)

resultdf_1 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_2 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_3 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_4 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))

# Setting up empty dataframes to save NMI scores across parameters
TA_NMI = data.frame(matrix(0, ncol = 20, nrow = 68))
CV_NMI = data.frame(matrix(0, ncol = 20, nrow = 69))
SA_NMI = data.frame(matrix(0, ncol = 20, nrow = 70))
behav_NMI = data.frame(matrix(0, ncol = 20, nrow = 4))

t = 10; 	# Number of Iterations, usually (10~20)

# Range of parameters tested
alphas_list <- c("0.3", "0.4", "0.5", "0.6", "0.7", "0.8")

# Importing the different data types as individual participant matrices 
participants <- read.csv("IDs.csv", header=FALSE)
TA <- read.csv("TA.csv", header=FALSE)
CV <- read.csv("CV.csv", header=FALSE)
SA <- read.csv("SA.csv", header=FALSE)
behav <- read.csv("Behavioural.csv", header=FALSE)

# Cleaning dfs
TA <- TA[,-1]
names(TA) <- TA[1,]
TA <- TA[-1,]
TA2 <- mutate_all(TA, function(x) as.numeric(as.character(x)))

CV <- CV[,-1]
names(CV) <- CV[1,]
CV <- CV[-1,]
CV2 <- mutate_all(CV, function(x) as.numeric(as.character(x)))

SA <- SA[,-1]
names(SA) <- SA[1,]
SA <- SA[-1,]
SA2 <- mutate_all(SA, function(x) as.numeric(as.character(x)))

behav <- behav[,-1]
names(behav) <- behav[1,]
behav <- behav[-1,]
behav2 <- mutate_all(behav, function(x) as.numeric(as.character(x)))

# Normalizing measures within each data type using a function from the SNF package
TA2 = standardNormalization(TA2)
CV2 = standardNormalization(CV2)
SA2 = standardNormalization(SA2)
behav2 = standardNormalization(behav2)

# Creating participant distance matrices using euclidean distances
Dist_TA = dist2(as.matrix(TA2),as.matrix(TA2));
Dist_CV = dist2(as.matrix(CV2),as.matrix(CV2));
Dist_SA = dist2(as.matrix(SA2),as.matrix(SA2));
Dist_behav = dist2(as.matrix(behav2),as.matrix(behav2));

### for each parameter combination: SNF is used to create a similarity matrix and (1) the optimal number of clusters as well as the normalized mutual information (NMI) between each measure and the matrix are calculated and recorded
for(param in 1:length(alphas_list)){
  alpha <- as.numeric(as.character(alphas_list[param]))
  print(alpha)
  for(test in 1:nrow(resultdf_1)){
    K <- resultdf_1$K[test]
    print(test)
    
    # Creating participant affinity matrices for each data type
    AM_TA = affinityMatrix(Dist_TA,K, alpha)
    AM_CV = affinityMatrix(Dist_CV,K,alpha)
    AM_SA = affinityMatrix(Dist_SA,K, alpha)
    AM_behav = affinityMatrix(Dist_behav,K, alpha)
    
    # Using SNF to create fused similarity matrix
    SNF1 = SNF(list(AM_TA, AM_CV, AM_SA, AM_behav), K, t)
    
    # Calculating optimal number of clusters given the created similarity matrix
    EstClust_1 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[1]][1] # or [[1]][2]
    EstClust_2 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[2]][1] # or [[1]][2]
    EstClust_3 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[3]][1] # or [[1]][2]
    EstClust_4 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[4]][1] # or [[1]][2]
    
    # Adding cluster number to tracking file for the given parameter
    resultdf_1[test,(param +1)]<-EstClust_1
    resultdf_2[test,(param +1)]<-EstClust_2
    resultdf_3[test,(param +1)]<-EstClust_3
    resultdf_4[test,(param +1)]<-EstClust_4
    
    # Calculating NMI scores based on similarity matrix
    SNF1_NMIScores <-rankFeaturesByNMI(list(TA2, CV2, SA2, behav2), SNF1)
    
    # Organizing NMI scores
    TA_NMI[, test] <- SNF1_NMIScores[[1]][1]
    CV_NMI[, test] <- SNF1_NMIScores[[1]][2]
    SA_NMI[, test] <- SNF1_NMIScores[[1]][3]
    behav_NMI[, test] <- SNF1_NMIScores[[1]][4]
    
    # Calculating groups using spectral clustering
    Group <-spectralClustering(SNF1, 4)
    
  }
  
  # Setting directory and saving NMI score files
  directory <- (paste("output/1_SNF_analysis/parameter_iterations/alpha_", alpha, sep=""))
  write.csv(TA_NMI, file=file.path(directory, paste("/TA_snf.csv", sep="")))
  write.csv(CV_NMI, file=file.path(directory, paste("/CV_snf.csv", sep="")))
  write.csv(SA_NMI, file=file.path(directory, paste("/SA_snf.csv", sep="")))
  write.csv(behav_NMI, file=file.path(directory, paste("/behav_snf.csv", sep="")))
  all_scores <- rbind(TA_NMI, CV_NMI)
  all_scores <- rbind(all_scores, SA_NMI)
  all_scores <- rbind(all_scores, behav_NMI)
  write.csv(all_scores, file=file.path(directory, paste("/all_scores.csv", sep="")))
}

# Saving estimated number of clusters files
directory <- ("output/1_SNF_analysis/parameter_iterations/")

write.csv(resultdf_1, file=file.path(directory, paste("Estnumclus_1.csv", sep="")))
write.csv(resultdf_2, file=file.path(directory, paste("Estnumclus_2.csv", sep="")))
write.csv(resultdf_3, file=file.path(directory, paste("Estnumclus_3.csv", sep="")))
write.csv(resultdf_4, file=file.path(directory, paste("Estnumclus_4.csv", sep="")))
