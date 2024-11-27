# Including libraries

library(effsize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(data.table)
library(ggsignif)
library(dunn.test)
library(knitr)
library(png)
library(reshape2)
library(haven)
library(ggpubr)
library(lsr)
library(car)
library(sjstats)
library(agricolae)

install.packages("kableExtra")
library(kableExtra)

install.packages("effectsize")
library(effectsize)

# Import file
NMI <- read.csv("output/1_SNF_analysis/all_scores_k18_0.8_1000perms.csv")
colnames(NMI) <- c("Measure", "NMI")

#--- Cluster analysis, n=515 ---#

# Measures included: TA, CV, SA from freesurfer , 4 clinical scores from 4 scales   
# Parameters below chosen based on SNF guidelines and increasing the NMI scores each data type has with the fused similarity network   
# K=18; iterations=10; hyperparameter=0.8; No individuals cluster unreliably 

#--- Similarity Network Fusion Clusters ---#

# Demographics of groups and their symptom severity for clinical symptoms

Demographics1 <- filter(all) %>%
  group_by(groups) %>% 
  summarise(N=length(groups), Age=round(mean(AGE, na.rm = TRUE), 2), Age_sd=round(sd(AGE, na.rm=TRUE), 2), Education= round(mean(PTEDUCAT, na.rm=TRUE), 2), EducationSD= round(sd(PTEDUCAT, na.rm=TRUE), 2),Depres_score = round(mean(GDTOTAL, na.rm = TRUE), 2),Depres_scoreSD = round(sd(GDTOTAL, na.rm = TRUE), 2), Hachinski =round(mean(HMSCORE, na.rm = TRUE), 2), HachinskiSD =round(sd(HMSCORE, na.rm = TRUE), 2))

Demographics1
directory <- ("C:\\Users\\bourq\\Documents\\0-Carleton\\6 - PhD\\Y2\\0-MR(John)\\SNF\\SNFcode")
write.csv(Demographics1, file=file.path(directory, paste("/Demos_DD.csv", sep="")))

# Categorical demos
all$groups <- as.factor(all$groups)
gr1 <- filter(all, groups == "1")
gr2 <- filter(all, groups == "2")
gr3 <- filter(all, groups == "3")
gr4 <- filter(all, groups == "4")

# Gender
table(gr1$PTGENDER, exclude = NULL)
table(gr2$PTGENDER, exclude = NULL)
table(gr3$PTGENDER, exclude = NULL)
table(gr4$PTGENDER, exclude = NULL)

# Smoking
table(gr1$MH16SMOK, exclude = NULL)
table(gr2$MH16SMOK, exclude = NULL)
table(gr3$MH16SMOK, exclude = NULL)
table(gr4$MH16SMOK, exclude = NULL)

# Alcohol abuse
table(gr1$MH14ALCH, exclude = NULL)
table(gr2$MH14ALCH, exclude = NULL)
table(gr3$MH14ALCH, exclude = NULL)
table(gr4$MH14ALCH, exclude = NULL)

# Cardiovascular
table(gr1$MH4CARD, exclude = NULL)
table(gr2$MH4CARD, exclude = NULL)
table(gr3$MH4CARD, exclude = NULL)
table(gr4$MH4CARD, exclude = NULL)

# Stat differences between groups...

# ...and gender
chisq.test(all$groups, all$PTGENDER, correct=FALSE)

# ...and cardiovascular
chisq.test(all$groups, all$MH4CARD, correct=FALSE)

# ...and smoking
chisq.test(all$groups, all$MH16SMOK, correct=FALSE)

# ...and alcohol abuse
chisq.test(all$groups, all$MH14ALCH, correct=FALSE)

# ...and age
kruskal.test(AGE ~ groups, data = all)

# ...and education
kruskal.test(PTEDUCAT ~ groups, data = all)

# ...and depression
kruskal.test(GDTOTAL ~ groups, data = all)

# ...and Hachinski score
kruskal.test(HMSCORE ~ groups, data = all)

# Demographics of clinical diagnoses

Demographics2 <- filter(all) %>%
  group_by(DIAGNOSIS) %>%
  summarise(N=length(groups), Age=round(mean(AGE, na.rm = TRUE), 2), Age_sd=round(sd(AGE, na.rm=TRUE), 2), Education= round(mean(PTEDUCAT, na.rm=TRUE), 2), EducationSD= round(sd(PTEDUCAT, na.rm=TRUE), 2),Depres_score = round(mean(GDTOTAL, na.rm = TRUE), 2),Depres_scoreSD = round(sd(GDTOTAL, na.rm = TRUE), 2), Hachinski =round(mean(HMSCORE, na.rm = TRUE), 2), HachinskiSD =round(sd(HMSCORE, na.rm = TRUE), 2))

Demographics2
write.csv(Demographics2, file=file.path(directory, paste("/Demos_Clin.csv", sep="")))

# Clinical demos
all$DIAGNOSIS <- as.factor(all$DIAGNOSIS)
grMCI <- filter(all, DIAGNOSIS == "MCI")
grCN <- filter(all, DIAGNOSIS == "CN")
grAD <- filter(all, DIAGNOSIS == "Dementia")

# Gender
table(grMCI$PTGENDER, exclude = NULL)
table(grCN$PTGENDER, exclude = NULL)
table(grAD$PTGENDER, exclude = NULL)

# stat differences between groups...

# ...and gender
chisq.test(all$DIAGNOSIS, all$PTGENDER, correct=FALSE)

# ...and age
kruskal.test(AGE ~ DIAGNOSIS, data = all)

# ...and education
kruskal.test(PTEDUCAT ~ DIAGNOSIS, data = all)

# ...and depression
kruskal.test(GDTOTAL ~ DIAGNOSIS, data = all)

# ...and Hachinski score
kruskal.test(HMSCORE ~ DIAGNOSIS, data = all)

#--- Creating statistics table for all SNF features for data-driven groups ---#

DataDriven <- filter(all_features) %>%
  group_by(Measure) %>%
  do(tidy(Anova(lm(brain_features ~AGE + GDTOTAL + PTGENDER + groups, data=.,), type=2)))

# Calculating the Eta sq for each feature
names <- unique(all_features$Measure)
effects <- c("term", "etasq")
for (idx in 1:length(unique(all_features$Measure))){
  print(names[idx])
  stats <- effectsize::eta_squared(Anova(lm(all[ ,c(names[idx])] ~AGE + GDTOTAL + PTGENDER + groups, data=all)))
  effects <- rbind(effects, stats)
}

effects <- effects[which(effects$Parameter == "groups"), ]
NMI$Measure <- names

NE <- cbind(names, effects)
colnames(NE)[1] <- "Measure"
colnames(NE)[3] <- "etasq"

NE$etasq <- signif(as.numeric(NE$etasq), digits = 2)

# Organizing table with other statistics for each measure
All_results <- DataDriven
All_results <- All_results %>% filter(term == "groups") 
All_results$p.fdr <- p.adjust(All_results$p.value, method="fdr")

All_results <- All_results %>% arrange(p.value)
All_results$statistic <- signif(All_results$statistic, digits = 3)
All_results$p.value <- signif(All_results$p.value, digits = 3)
All_results$df <- NULL

top_contributing_DD <- All_results

top_contributing_DD <- merge(NMI, top_contributing_DD, by="Measure")
top_contributing_DD <- merge(top_contributing_DD, NE, by="Measure")
top_contributing_DD <- subset(top_contributing_DD, select = c(-term, -Parameter))

# Order by NMI
top_contributing_DD <- top_contributing_DD[order(-top_contributing_DD$NMI), ]

# Formatting
top_contributing_DD$NMI <- format(round(top_contributing_DD$NMI,3), nsmall = 3)
top_contributing_DD$p.fdr <- format(round(top_contributing_DD$p.fdr,3), nsmall = 3)
top_contributing_DD$etasq <- format(round(top_contributing_DD$etasq,3), nsmall = 3)

# Add labels to all measures
all_labels <- read.csv("output/2_Comparing_subgroups/All_Labels2.csv")
DDL <- merge(all_labels, top_contributing_DD, by="Measure")

# to make DF for Table2
DD_df <- subset(DDL, select = c(-Label, -Mlab, -sumsq, -p.value, - CI, -CI_low, -CI_high))
write.csv(DD_df, file="output/2_Comparing_subgroups/DataDriven.csv")

#--- Creating statistics table for all SNF features for clinical diagnoses - ALL FEATURES ---#

Clinical <- filter(all_features) %>%
  group_by(Measure) %>%
  do(tidy(Anova(lm(brain_features ~AGE + GDTOTAL + PTGENDER + DIAGNOSIS, data=.,), type=2)))

# Calculating the Eta sq for each feature
names <- unique(all_features$Measure)
effects <- c("term", "etasq")
for (idx in 1:length(unique(all_features$Measure))){
  print(names[idx])
  stats <- effectsize::eta_squared(Anova(lm(all[ ,c(names[idx])] ~AGE + GDTOTAL + PTGENDER + DIAGNOSIS, data=all)))
  effects <- rbind(effects, stats)
}

effects <- effects[which(effects$Parameter == "DIAGNOSIS"), ]
NMI$Measure <- names

NE <- cbind(names, effects)
colnames(NE)[1] <- "Measure"
colnames(NE)[3] <- "etasq"

NE$etasq <- signif(as.numeric(NE$etasq), digits = 2)

# Organizing table with other statistics for each measure
All_results <- Clinical
All_results <- All_results %>% filter(term == "DIAGNOSIS") 
All_results$p.fdr <- p.adjust(All_results$p.value, method="fdr")

All_results <- All_results %>% arrange(p.value)
All_results$statistic <- signif(All_results$statistic, digits = 3)
All_results$p.value <- signif(All_results$p.value, digits = 3)
All_results$df <- NULL

top_contributing_Clin <- All_results

top_contributing_Clin <- merge(NMI, top_contributing_Clin, by="Measure")
top_contributing_Clin <- merge(top_contributing_Clin, NE, by="Measure")
top_contributing_Clin <- subset(top_contributing_Clin, select = c(-term, -Parameter))

# Order by NMI
top_contributing_Clin <- top_contributing_Clin[order(-top_contributing_Clin$NMI), ]

# Formatting
top_contributing_Clin$NMI <- format(round(top_contributing_Clin$NMI,3), nsmall = 3)
top_contributing_Clin$p.fdr <- format(round(top_contributing_Clin$p.fdr,3), nsmall = 3)
top_contributing_Clin$etasq <- format(round(top_contributing_Clin$etasq,3), nsmall = 3)

# Add labels to all measures
all_labels <- read.csv("output/2_Comparing_subgroups/All_Labels2.csv")
Clin <- merge(all_labels, top_contributing_Clin, by="Measure")

Clin_df <- subset(Clin, select = c(-Label, -Mlab, -sumsq, -p.value, - CI, -CI_low, -CI_high))
write.csv(Clin_df, file="output/2_Comparing_subgroups/Clinical.csv")

#--- Grouping by DD groups ---#
DDgroup1 <- filter(all, groups == "1")
DDgroup2 <- filter(all, groups == "2")
DDgroup3 <- filter(all, groups == "3")
DDgroup4 <- filter(all, groups == "4")

# Group 1
DD1_CV <- colMeans(DDgroup1 %>% dplyr::select(ends_with('CV')))
DD1_SA <- colMeans(DDgroup1 %>% dplyr::select(ends_with('SA')))            
DD1_TA <- colMeans(DDgroup1 %>% dplyr::select(ends_with('TA')))

# Group 2
DD2_CV <- colMeans(DDgroup2 %>% dplyr::select(ends_with('CV')))
DD2_SA <- colMeans(DDgroup2 %>% dplyr::select(ends_with('SA')))            
DD2_TA <- colMeans(DDgroup2 %>% dplyr::select(ends_with('TA')))

# Group 3
DD3_CV <- colMeans(DDgroup3 %>% dplyr::select(ends_with('CV')))
DD3_SA <- colMeans(DDgroup3 %>% dplyr::select(ends_with('SA')))           
DD3_TA <- colMeans(DDgroup3 %>% dplyr::select(ends_with('TA')))

# Group 4
DD4_CV <- colMeans(DDgroup4 %>% dplyr::select(ends_with('CV')))
DD4_SA <- colMeans(DDgroup4 %>% dplyr::select(ends_with('SA')))             
DD4_TA <- colMeans(DDgroup4 %>% dplyr::select(ends_with('TA')))

write.csv(DD1_CV, file="output/2_Comparing_subgroups/DD1_CV.csv")
write.csv(DD2_CV, file="output/2_Comparing_subgroups/DD2_CV.csv")
write.csv(DD3_CV, file="output/2_Comparing_subgroups/DD3_CV.csv")
write.csv(DD4_CV, file="output/2_Comparing_subgroups/DD4_CV.csv")

write.csv(DD1_TA, file="output/2_Comparing_subgroups/DD1_TA.csv")
write.csv(DD2_TA, file="output/2_Comparing_subgroups/DD2_TA.csv")
write.csv(DD3_TA, file="output/2_Comparing_subgroups/DD3_TA.csv")
write.csv(DD4_TA, file="output/2_Comparing_subgroups/DD4_TA.csv")

#--- Grouping by Clin diagnoses ---#
Clin_MCI <- filter(all, DIAGNOSIS == "MCI")
Clin_CN <- filter(all, DIAGNOSIS == "CN")
Clin_Dementia <- filter(all, DIAGNOSIS == "Dementia")

# MCI
MCI_CV <- colMeans(Clin_MCI %>% dplyr::select(ends_with('CV')))
MCI_SA <- colMeans(Clin_MCI %>% dplyr::select(ends_with('SA')))           
MCI_TA <- colMeans(Clin_MCI %>% dplyr::select(ends_with('TA')))

# CN
CN_CV <- colMeans(Clin_CN %>% dplyr::select(ends_with('CV')))
CN_SA <- colMeans(Clin_CN %>% dplyr::select(ends_with('SA')))             
CN_TA <- colMeans(Clin_CN %>% dplyr::select(ends_with('TA')))

# Dementia
Dementia_CV <- colMeans(Clin_Dementia %>% dplyr::select(ends_with('CV')))
Dementia_SA <- colMeans(Clin_Dementia %>% dplyr::select(ends_with('SA')))           
Dementia_TA <- colMeans(Clin_Dementia %>% dplyr::select(ends_with('TA')))

all_CV <- all %>% dplyr::select(ends_with('CV'))
all_SA <- all %>% dplyr::select(ends_with('SA'))
all_TA <- all %>% dplyr::select(ends_with('TA'))

write.csv(all_CV, file="output/2_Comparing_subgroups/all_CV.csv")
write.csv(CN_CV, file="output/2_Comparing_subgroups/CN_CV.csv")
write.csv(Dementia_CV, file="output/2_Comparing_subgroups/Dementia_CV.csv")
write.csv(MCI_CV, file="output/2_Comparing_subgroups/MCI_CV.csv")

#--- Calculating the Eta sq for each feature ---#
names <- unique(all_features$Measure)
effects <- c("term", "etasq")
for (idx in 1:length(unique(all_features$Measure))){
  print(names[idx])
  stats <- effectsize::eta_squared(Anova(lm(all[ ,c(names[idx])] ~AGE + GDTOTAL + PTGENDER + groups, data=all)))
  effects <- rbind(effects, stats)
}

effects <- effects[which(effects$Parameter == "groups"), ]
NMI$Measure <- names

NE <- cbind(names, effects)
colnames(NE)[1] <- "Measure"
colnames(NE)[3] <- "etasq"

NE$etasq <- signif(as.numeric(NE$etasq), digits = 2)

# Organizing table with other statistics for each measure
All_results <- DataDriven
All_results <- All_results %>% filter(term == "groups") 
All_results$p.fdr <- p.adjust(All_results$p.value, method="fdr")

All_results <- All_results %>% arrange(p.value)
All_results$statistic <- signif(All_results$statistic, digits = 3)
All_results$p.value <- signif(All_results$p.value, digits = 3)
All_results$df <- NULL

top_contributing_DD <- All_results

top_contributing_DD <- merge(NMI, top_contributing_DD, by="Measure")
top_contributing_DD <- merge(top_contributing_DD, NE, by="Measure")
top_contributing_DD <- subset(top_contributing_DD, select = c(-term, -Parameter))

# Order by NMI
top_contributing_DD <- top_contributing_DD[order(-top_contributing_DD$NMI), ]

# Formatting
top_contributing_DD$NMI <- signif(top_contributing_DD$NMI, digits = 3)
top_contributing_DD$p.fdr <- signif(top_contributing_DD$p.fdr, digits = 3)
top_contributing_DD$p.fdr <- ifelse(top_contributing_DD$p.fdr > 0.05, "n.s", top_contributing_DD$p.fdr)
top_contributing_DD$p.fdr <- format(top_contributing_DD$p.fdr, scientific = TRUE)

# Add labels to all measures
all_labels <- read.csv("output/2_Comparing_subgroups/All_Labels2.csv")
DDL <- merge(all_labels, top_contributing_DD, by="Measure")

write.csv(DDL, file="output/2_Comparing_subgroups/all_snf_measures_groupdiffs.csv")

write.csv(Clin_df, file="output/2_Comparing_subgroups/Clin_df.csv")
write.csv(DD_df, file="output/2_Comparing_subgroups/DD_df.csv")
save(DD_df, file='DD_df.rda')
save(Clin_df, file='Clin_df.rda')
