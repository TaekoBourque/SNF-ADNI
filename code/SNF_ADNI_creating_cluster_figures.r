# Including libraries

library(ggseg)
library(dplyr)
library(tidyr)
library(haven)
library(ggpubr)
library(plotrix)
library(tidyverse)

install.packages("igraph")
library(igraph)

# Importing file and formatting
ALL <- read.csv("output/2_Comparing_subgroups/all_snf_measures_groupdiffs.csv")
ALL <- subset(ALL, select=c(-Label))
colnames(ALL)[5] <- "label"

# Cortical thickness average
aCT <- dplyr::filter(ALL, grepl("TA",Measure))
aCT$region <- as.character(aCT$region)
aCT$hemi <- as.character(aCT$hemi)
aCT$NMI <- round(aCT$NMI, digits=3)

# Cortical volume
aCV <- dplyr::filter(ALL, grepl("CV",Measure))
aCV$region <- as.character(aCV$region)
aCV$hemi <- as.character(aCV$hemi)
aCV$NMI <- round(aCV$NMI, digits=3)

# Subcortical volume
aSV <- dplyr::filter(ALL, grepl("SV",Measure))
aSV$region <- as.character(aSV$region)
aSV$hemi <- as.character(aSV$hemi)
aSV$NMI <- round(aSV$NMI, digits=3)

#--- Pie chart break down of the 4 clusters by diagnosis ---#

freq = table(all$DIAGNOSIS)
lbls = paste(names(freq), "-", round((freq/515)*100, 1), "%", sep="")

setwd("output/2_Comparing_subgroups/Figures/Figure1")
png('pie_all.png')
pie(freq, main="All participants", labels = lbls, cex=1.5, col=c("gray28", "darkgray", "white"))
dev.off()

# DD1
freq = table(clust1$DIAGNOSIS)
lbls = paste(names(freq), "-", round((freq/39)*100, 1), "%", sep="")
png('pie_1.png')
pie(freq, main="Group 1", labels = lbls,cex=1.5, col=c("gray28", "darkgray", "white"), lwd= 3)
dev.off()

# DD2
freq = table(clust2$DIAGNOSIS)
lbls = paste(names(freq), "-", round((freq/179)*100, 1), "%", sep="")
png('pie_2.png')
pie(freq, main="Group 2", labels = lbls, cex=1.5, col=c("gray28", "darkgray", "white"))
dev.off()

# DD3
freq = table(clust3$DIAGNOSIS)
lbls = paste(names(freq), "-", round((freq/245)*100, 1), "%", sep="")
png('pie_3.png')
pie(freq, main="Group 3", labels = lbls,cex=1.5, col=c("gray28", "darkgray", "white"))
dev.off()

# DD4
freq = table(clust4$DIAGNOSIS)
lbls = paste(names(freq), "-", round((freq/52)*100, 1), "%", sep="")
png('pie_4.png')
pie(freq, main="Group 4", labels = lbls,cex=1.5, col=c("gray28", "darkgray", "white"))
dev.off()

#--- Visualization of the differences in clinical scores between groups ---#

# Plot of symptom scores by group
symptoms_stacked$Symptom <- factor(symptoms_stacked$Symptom, labels = c("ADAS", "MMSE", "MOCA", "CDR"))

# Getting the standard error of each behaviour measure
test <- filter(symptoms_stacked) %>%
  group_by(Symptom, groups) %>%
  summarise(std_err=std.error(severity, na.rm=TRUE), N=length(groups))

symptoms_stacked$std_error <- NA

for (i in 1:21){
  symp <- (test[i, c("Symptom")])
  gr <- as.character(test[i, 2])
  symptoms_stacked$std_error <- ifelse(symptoms_stacked$Symptom == symp$Symptom & symptoms_stacked$groups == gr, as.numeric(test[i, c("std_err")]), symptoms_stacked$std_error)
}

ggplot(symptoms_stacked, aes(x=Symptom, y=severity, fill=as.factor(groups))) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, face="bold", angle = 30, hjust = 1), axis.title.x=element_blank(), legend.title=element_text(size=18), axis.text.y=element_text(size=18), axis.title.y=element_text(size=20), legend.text=element_text(size=18)) +
  guides(fill=guide_legend(title="Groups")) +
  scale_fill_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
  labs(y= "Z-score") +
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group=groups), linetype = "dotted", width = 0.75, position = position_dodge())

ggsave("Group_behaviour_scores_residuals_stderrorbars.png", plot=last_plot(), width = 12, height = 6, units = c("in"), dpi=400, device = "png", path = "output/2_Comparing_subgroups/Figures/Figure1/")

#--- Plot of symptom scores by diagnosis ---#

symptoms_stacked$Symptom <- factor(symptoms_stacked$Symptom, labels = c("ADAS", "MMSE", "MOCA", "CDR"))

test <- filter(symptoms_stacked) %>%
  group_by(Symptom, DIAGNOSIS) %>%
  summarise(std_err=std.error(severity, na.rm=TRUE), N=length(DIAGNOSIS))

symptoms_stacked$std_error <- NA

for (i in 1:21){
  symp <-(test[i, c("Symptom")])
  gr <- as.character(test[i, 2])
  symptoms_stacked$std_error <- ifelse(symptoms_stacked$Symptom == symp$Symptom & symptoms_stacked$groups == gr, as.numeric(test[i, c("std_err")]), symptoms_stacked$std_error)
}

ggplot(symptoms_stacked, aes(x=Symptom, y=severity, fill=DIAGNOSIS)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(size=18, face="bold", angle = 30, hjust = 1), axis.title.x=element_blank(), legend.title=element_text(size=18), axis.text.y=element_text(size=18), axis.title.y=element_text(size=20), legend.text=element_text(size=18)) +
  guides(fill=guide_legend(title="Diagnosis")) +
scale_fill_manual(values=c("gray28", "darkgray", "white")) +
  labs(y= "Z-score") +
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group=DIAGNOSIS), linetype = "dotted", width = 0.75, position = position_dodge())

ggsave("Dx_behaviour_scores_residuals_stderrorbars.png", plot=last_plot(), width = 12, height = 6, units = c("in"), dpi=400, device = "png", path ="output/2_Comparing_subgroups/Figures/Figure1/")

#--- Supplementary Figures ---#

# Age of groups
my_comparisons <- list(  c("1", "3"), c("2", "3"), c("4", "3"))

ggplot(all, aes(x=groups, y=AGE, color=groups)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = groups), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "Age (years)", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("gray28", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#0072b5", "#e18727", "#20854e")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = TRUE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/age_groups.png", plot = last_plot())

#--- ggseg ---#

# CV
all_labels <- read.csv("output/2_Comparing_subgroups/All_Labels2.csv")

allCV_icv <-read.csv("output/1_SNF_analysis/CV.csv")
allTA <-read.csv("output/1_SNF_analysis/TA.csv")
allSA_icv <-read.csv("output/1_SNF_analysis/SA.csv")
  
gLab <- all$groups
dLab <- all$DIAGNOSIS
#1= CU; 2=MCI; 3=Dementia

zCV <- subset(allCV_icv, select = c(-X))
zCV <- apply(zCV,2,scale)
CVlabs <- as.data.frame(cbind(dLab, gLab, zCV))
names(CVlabs)[names(CVlabs) == 'gLab'] <- 'groups'
names(CVlabs)[names(CVlabs) == 'dLab'] <- 'DIAGNOSIS'

zTA <- subset(allTA, select = c(-X))
zTA <- apply(zTA,2,scale)
TAlabs <- as.data.frame(cbind(dLab,gLab,zTA))

zSA <- subset(allSA_icv, select = c(-X))
zSA <- apply(zSA,2,scale)
SAlabs <- as.data.frame(cbind(dLab,gLab,zSA))

# Grouping by DD groups
zDD1 <- filter(CVlabs, as.numeric(CVlabs$groups) == "1")
zDD1 <- subset(zDD1, select=c(-groups, -DIAGNOSIS))
zDD1 <- mutate_all(zDD1, function(x) as.numeric(as.character(x)))
zDD1 <- as.data.frame(colMeans(zDD1))
zDD1 <- tibble::rownames_to_column(zDD1, "Measure")
zDD1 <- merge(all_labels, zDD1, by="Measure")
names(zDD1)[names(zDD1) == 'colMeans(zDD1)'] <- 'zscore'

zDD2 <- filter(CVlabs, as.numeric(CVlabs$groups) == "2")
zDD2 <- subset(zDD2, select=c(-groups, -DIAGNOSIS))
zDD2 <- mutate_all(zDD2, function(x) as.numeric(as.character(x)))
zDD2 <- as.data.frame(colMeans(zDD2))
zDD2 <- tibble::rownames_to_column(zDD2, "Measure")
zDD2 <- merge(all_labels, zDD2, by="Measure")
names(zDD2)[names(zDD2) == 'colMeans(zDD2)'] <- 'zscore'

zDD3 <- filter(CVlabs, as.numeric(CVlabs$groups) == "3")
zDD3 <- subset(zDD3, select=c(-groups, -DIAGNOSIS))
zDD3 <- mutate_all(zDD3, function(x) as.numeric(as.character(x)))
zDD3 <- as.data.frame(colMeans(zDD3))
zDD3 <- tibble::rownames_to_column(zDD3, "Measure")
zDD3 <- merge(all_labels, zDD3, by="Measure")
names(zDD3)[names(zDD3) == 'colMeans(zDD3)'] <- 'zscore'

zDD4 <- filter(CVlabs, as.numeric(CVlabs$groups) == "4")
zDD4 <- subset(zDD4, select=c(-groups, -DIAGNOSIS))
zDD4 <- mutate_all(zDD4, function(x) as.numeric(as.character(x)))
zDD4 <- as.data.frame(colMeans(zDD4))
zDD4 <- tibble::rownames_to_column(zDD4, "Measure")
zDD4 <- merge(all_labels, zDD4, by="Measure")
names(zDD4)[names(zDD4) == 'colMeans(zDD4)'] <- 'zscore'

# Grouping by Clinical diagnosis

zCN <- filter(CVlabs, as.numeric(CVlabs$DIAGNOSIS) == "1")
zCN <- subset(zCN, select=c(-groups, -DIAGNOSIS))
zCN <- mutate_all(zCN, function(x) as.numeric(as.character(x)))
zCN <- as.data.frame(colMeans(zCN))
zCN <- tibble::rownames_to_column(zCN, "Measure")
zCN <- merge(all_labels, zCN, by="Measure")
names(zCN)[names(zCN) == 'colMeans(zCN)'] <- 'zscore'

zMCI <- filter(CVlabs, as.numeric(CVlabs$DIAGNOSIS) == "2")
zMCI <- subset(zMCI, select=c(-groups, -DIAGNOSIS))
zMCI <- mutate_all(zMCI, function(x) as.numeric(as.character(x)))
zMCI <- as.data.frame(colMeans(zMCI))
zMCI <- tibble::rownames_to_column(zMCI, "Measure")
zMCI <- merge(all_labels, zMCI, by="Measure")
names(zMCI)[names(zMCI) == 'colMeans(zMCI)'] <- 'zscore'

zDementia <- filter(CVlabs, as.numeric(CVlabs$DIAGNOSIS) == "3")
zDementia <- subset(zDementia, select=c(-groups, -DIAGNOSIS))
zDementia <- mutate_all(zDementia, function(x) as.numeric(as.character(x)))
zDementia <- as.data.frame(colMeans(zDementia))
zDementia <- tibble::rownames_to_column(zDementia, "Measure")
zDementia <- merge(all_labels, zDementia, by="Measure")
names(zDementia)[names(zDementia) == 'colMeans(zDementia)'] <- 'zscore'

# TA
                        
# Grouping by DD groups
TzDD1 <- filter(TAlabs, as.numeric(TAlabs$gLab) == "1")
TzDD1 <- subset(TzDD1, select=c(-gLab, -dLab))
TzDD1 <- mutate_all(TzDD1, function(x) as.numeric(as.character(x)))
TzDD1 <- as.data.frame(colMeans(TzDD1))
TzDD1 <- tibble::rownames_to_column(TzDD1, "Measure")
TzDD1 <- merge(all_labels, TzDD1, by="Measure")
names(TzDD1)[names(TzDD1) == 'colMeans(TzDD1)'] <- 'zscore'

TzDD2 <- filter(TAlabs, as.numeric(TAlabs$gLab) == "2")
TzDD2 <- subset(TzDD2, select=c(-gLab, -dLab))
TzDD2 <- mutate_all(TzDD2, function(x) as.numeric(as.character(x)))
TzDD2 <- as.data.frame(colMeans(TzDD2))
TzDD2 <- tibble::rownames_to_column(TzDD2, "Measure")
TzDD2 <- merge(all_labels, TzDD2, by="Measure")
names(TzDD2)[names(TzDD2) == 'colMeans(TzDD2)'] <- 'zscore'

TzDD3 <- filter(TAlabs, as.numeric(TAlabs$gLab) == "3")
TzDD3 <- subset(TzDD3, select=c(-gLab, -dLab))
TzDD3 <- mutate_all(TzDD3, function(x) as.numeric(as.character(x)))
TzDD3 <- as.data.frame(colMeans(TzDD3))
TzDD3 <- tibble::rownames_to_column(TzDD3, "Measure")
TzDD3 <- merge(all_labels, TzDD3, by="Measure")
names(TzDD3)[names(TzDD3) == 'colMeans(TzDD3)'] <- 'zscore'

TzDD4 <- filter(TAlabs, as.numeric(TAlabs$gLab) == "4")
TzDD4 <- subset(TzDD4, select=c(-gLab, -dLab))
TzDD4 <- mutate_all(TzDD4, function(x) as.numeric(as.character(x)))
TzDD4 <- as.data.frame(colMeans(TzDD4))
TzDD4 <- tibble::rownames_to_column(TzDD4, "Measure")
TzDD4 <- merge(all_labels, TzDD4, by="Measure")
names(TzDD4)[names(TzDD4) == 'colMeans(TzDD4)'] <- 'zscore'

# Grouping by Clinical diagnosis

TzCN <- filter(TAlabs, as.numeric(TAlabs$dLab) == "1")
TzCN <- subset(TzCN, select=c(-gLab, -dLab))
TzCN <- mutate_all(TzCN, function(x) as.numeric(as.character(x)))
TzCN <- as.data.frame(colMeans(TzCN))
TzCN <- tibble::rownames_to_column(TzCN, "Measure")
TzCN <- merge(all_labels, TzCN, by="Measure")
names(TzCN)[names(TzCN) == 'colMeans(TzCN)'] <- 'zscore'

TzMCI <- filter(TAlabs, as.numeric(TAlabs$dLab) == "2")
TzMCI <- subset(TzMCI, select=c(-gLab, -dLab))
TzMCI <- mutate_all(TzMCI, function(x) as.numeric(as.character(x)))
TzMCI <- as.data.frame(colMeans(TzMCI))
TzMCI <- tibble::rownames_to_column(TzMCI, "Measure")
TzMCI <- merge(all_labels, TzMCI, by="Measure")
names(TzMCI)[names(TzMCI) == 'colMeans(TzMCI)'] <- 'zscore'

TzDementia <- filter(TAlabs, as.numeric(TAlabs$dLab) == "3")
TzDementia <- subset(TzDementia, select=c(-gLab, -dLab))
TzDementia <- mutate_all(TzDementia, function(x) as.numeric(as.character(x)))
TzDementia <- as.data.frame(colMeans(TzDementia))
TzDementia <- tibble::rownames_to_column(TzDementia, "Measure")
TzDementia <- merge(all_labels, TzDementia, by="Measure")
names(TzDementia)[names(TzDementia) == 'colMeans(TzDementia)'] <- 'zscore'

# SA

# Grouping by DD groups
SzDD1 <- filter(SAlabs, as.numeric(SAlabs$gLab) == "1")
SzDD1 <- subset(SzDD1, select=c(-gLab, -dLab))
SzDD1 <- mutate_all(SzDD1, function(x) as.numeric(as.character(x)))
SzDD1 <- as.data.frame(colMeans(SzDD1))
SzDD1 <- tibble::rownames_to_column(SzDD1, "Measure")
SzDD1 <- merge(all_labels, SzDD1, by="Measure")
names(SzDD1)[names(SzDD1) == 'colMeans(SzDD1)'] <- 'zscore'

SzDD2 <- filter(SAlabs, as.numeric(SAlabs$gLab) == "2")
SzDD2 <- subset(SzDD2, select=c(-gLab, -dLab))
SzDD2 <- mutate_all(SzDD2, function(x) as.numeric(as.character(x)))
SzDD2 <- as.data.frame(colMeans(SzDD2))
SzDD2 <- tibble::rownames_to_column(SzDD2, "Measure")
SzDD2 <- merge(all_labels, SzDD2, by="Measure")
names(SzDD2)[names(SzDD2) == 'colMeans(SzDD2)'] <- 'zscore'

SzDD3 <- filter(SAlabs, as.numeric(SAlabs$gLab) == "3")
SzDD3 <- subset(SzDD3, select=c(-gLab, -dLab))
SzDD3 <- mutate_all(SzDD3, function(x) as.numeric(as.character(x)))
SzDD3 <- as.data.frame(colMeans(SzDD3))
SzDD3 <- tibble::rownames_to_column(SzDD3, "Measure")
SzDD3 <- merge(all_labels, SzDD3, by="Measure")
names(SzDD3)[names(SzDD3) == 'colMeans(SzDD3)'] <- 'zscore'

SzDD4 <- filter(SAlabs, as.numeric(SAlabs$gLab) == "4")
SzDD4 <- subset(SzDD4, select=c(-gLab, -dLab))
SzDD4 <- mutate_all(SzDD4, function(x) as.numeric(as.character(x)))
SzDD4 <- as.data.frame(colMeans(SzDD4))
SzDD4 <- tibble::rownames_to_column(SzDD4, "Measure")
SzDD4 <- merge(all_labels, SzDD4, by="Measure")
names(SzDD4)[names(SzDD4) == 'colMeans(SzDD4)'] <- 'zscore'

# Grouping by Clinical diagnosis

SzCN <- filter(SAlabs, as.numeric(SAlabs$dLab) == "1")
SzCN <- subset(SzCN, select=c(-gLab, -dLab))
SzCN <- mutate_all(SzCN, function(x) as.numeric(as.character(x)))
SzCN <- as.data.frame(colMeans(SzCN))
SzCN <- tibble::rownames_to_column(SzCN, "Measure")
SzCN <- merge(all_labels, SzCN, by="Measure")
names(SzCN)[names(SzCN) == 'colMeans(SzCN)'] <- 'zscore'

SzMCI <- filter(SAlabs, as.numeric(SAlabs$dLab) == "2")
SzMCI <- subset(SzMCI, select=c(-gLab, -dLab))
SzMCI <- mutate_all(SzMCI, function(x) as.numeric(as.character(x)))
SzMCI <- as.data.frame(colMeans(SzMCI))
SzMCI <- tibble::rownames_to_column(SzMCI, "Measure")
SzMCI <- merge(all_labels, SzMCI, by="Measure")
names(SzMCI)[names(SzMCI) == 'colMeans(SzMCI)'] <- 'zscore'

SzDementia <- filter(SAlabs, as.numeric(SAlabs$dLab) == "3")
SzDementia <- subset(SzDementia, select=c(-gLab, -dLab))
SzDementia <- mutate_all(SzDementia, function(x) as.numeric(as.character(x)))
SzDementia <- as.data.frame(colMeans(SzDementia))
SzDementia <- tibble::rownames_to_column(SzDementia, "Measure")
SzDementia <- merge(all_labels, SzDementia, by="Measure")
names(SzDementia)[names(SzDementia) == 'colMeans(SzDementia)'] <- 'zscore'

##### FINAL GGSSEG ##### --------------------------------------------------------------------------------

### ----- CV -------

#for DD1------------------------------
results <- data.frame(cbind(region=as.character(zDD1$region), z=zDD1$zscore, hemi=as.character(zDD1$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_DD1.png", plot = last_plot())

#for DD2------------------------------
results <- data.frame(cbind(region=as.character(zDD2$region), z=zDD2$zscore, hemi=as.character(zDD2$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_DD2.png", plot = last_plot()) 

#for DD3------------------------------
results <- data.frame(cbind(region=as.character(zDD3$region), z=zDD3$zscore, hemi=as.character(zDD3$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_DD3.png", plot = last_plot())

#for DD4------------------------------
results <- data.frame(cbind(region=as.character(zDD4$region), z=zDD4$zscore, hemi=as.character(zDD4$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_DD4.png", plot = last_plot())

#for CN------------------------------
results <- data.frame(cbind(region=as.character(zCN$region), z=zCN$zscore, hemi=as.character(zCN$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_CN.png", plot = last_plot())

#for MCI------------------------------
results <- data.frame(cbind(region=as.character(zMCI$region), z=zMCI$zscore, hemi=as.character(zMCI$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_MCI.png", plot = last_plot())

#for Dementia------------------------------
results <- data.frame(cbind(region=as.character(zDementia$region), z=zDementia$zscore, hemi=as.character(zDementia$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/CV_Dementia.png", plot = last_plot())


### ----- TA -------

#for DD1------------------------------
results <- data.frame(cbind(region=as.character(TzDD1$region), z=TzDD1$zscore, hemi=as.character(TzDD1$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_DD1.png", plot = last_plot())

#for DD2------------------------------
results <- data.frame(cbind(region=as.character(TzDD2$region), z=TzDD2$zscore, hemi=as.character(TzDD2$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_DD2.png", plot = last_plot())

#for DD3------------------------------
results <- data.frame(cbind(region=as.character(TzDD3$region), z=TzDD3$zscore, hemi=as.character(TzDD3$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_DD3.png", plot = last_plot())

#for DD4------------------------------
results <- data.frame(cbind(region=as.character(TzDD4$region), z=TzDD4$zscore, hemi=as.character(TzDD4$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_DD4.png", plot = last_plot())

#for CN------------------------------
results <- data.frame(cbind(region=as.character(TzCN$region), z=TzCN$zscore, hemi=as.character(TzCN$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_CN.png", plot = last_plot())

#for MCI------------------------------
results <- data.frame(cbind(region=as.character(TzMCI$region), z=TzMCI$zscore, hemi=as.character(TzMCI$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_MCI.png", plot = last_plot())

#for Dementia------------------------------
results <- data.frame(cbind(region=as.character(TzDementia$region), z=TzDementia$zscore, hemi=as.character(TzDementia$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/TA_Dementia.png", plot = last_plot())

### ----- SA -------

#for DD1------------------------------
results <- data.frame(cbind(region=as.character(SzDD1$region), z=SzDD1$zscore, hemi=as.character(SzDD1$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_DD1.png", plot = last_plot())

#for DD2------------------------------
results <- data.frame(cbind(region=as.character(SzDD2$region), z=SzDD2$zscore, hemi=as.character(SzDD2$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_DD2.png", plot = last_plot())

#for DD3------------------------------
results <- data.frame(cbind(region=as.character(SzDD3$region), z=SzDD3$zscore, hemi=as.character(SzDD3$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_DD3.png", plot = last_plot())

#for DD4------------------------------
results <- data.frame(cbind(region=as.character(SzDD4$region), z=SzDD4$zscore, hemi=as.character(SzDD4$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
   scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_DD4.png", plot = last_plot())

#for CN------------------------------
results <- data.frame(cbind(region=as.character(SzCN$region), z=SzCN$zscore, hemi=as.character(SzCN$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_CN.png", plot = last_plot())

#for MCI------------------------------
results <- data.frame(cbind(region=as.character(SzMCI$region), z=SzMCI$zscore, hemi=as.character(SzMCI$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_MCI.png", plot = last_plot())

#for Dementia------------------------------
results <- data.frame(cbind(region=as.character(SzDementia$region), z=SzDementia$zscore, hemi=as.character(SzDementia$hemi)))
results$z <- as.numeric(as.character(results$z))

results %>% 
  ggseg(atlas=dk, mapping=aes(fill=z), colour = "white") +
  scale_fill_gradientn(colours = c("#01273D", "#003657", "#004670", "#025080","#005C92","#0072B5", "#1C79B4", "#308BC4","#4E97C4", "#65A8D1","#8BC8ED","#A3D4F1","#C3E9FF", "white","lightgoldenrod1", "yellow2","gold","orange", "darkorange2","orangered3","red3"), na.value="grey", breaks=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7),labels=c(-1.5, -1.4, -1.3, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(-1.5, 0.7))
ggsave("output/2_Comparing_subgroups/Figures/Supplement/SA_Dementia.png", plot = last_plot())


#--- Age of groups ---#

my_comparisons <- list(  c("1", "2"), c("1", "3"),c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))

ggplot(all, aes(x=group, y=AGE, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "Age (years)", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/age_groups.png", plot = last_plot())

#--- ADAS ---#
my_comparisons <- list(  c("1", "2"), c("1", "3"),c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))

ggplot(all, aes(x=group, y=TOTAL13, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "ADAS", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/ADAS_groups.png", plot = last_plot())

#--- MMSE ----#
ggplot(all, aes(x=group, y=MMSCORE, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "MMSE", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/MMSE_groups.png", plot = last_plot())

#--- MoCA ----#
ggplot(all, aes(x=group, y=MOCA, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "MoCA", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/MoCA_groups.png", plot = last_plot())

#--- CDR ---#
ggplot(all, aes(x=group, y=CDRSB, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "CDR", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/CDR_groups.png", plot = last_plot())

#--- GENDER ---#

ggplot(all, aes(x=DIAGNOSIS, fill=PTGENDER)) +
  theme_classic() +
  geom_bar(stat = "Count") +
  labs(y= "Count", x="Group") +
  facet_wrap(~group, nrow = 1)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/sex_distribution.png", plot = last_plot(), width = 7, height =3)

#--- p-tau ----#
ptau <- merge(ugotptau181, ids, by="RID")
#remove duplicates by RID, keeping last value
ptau <- ptau[order(ptau$VISCODE ,decreasing=FALSE),]
ptau <- ptau %>% distinct(RID, .keep_all = TRUE)
ptau_val <- subset(ptau, select=c(RID,PLASMAPTAU181))
names(ptau_val)[names(ptau_val) == 'RID'] <- 'participant'

allTest <- all
# merging together all the data
allTest <- merge(allTest, ptau_val, by="participant", all.x=TRUE)

# ptau
ggplot(allTest, aes(x=group, y=PLASMAPTAU181, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "UoG P-Tau 181", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/UoGPtau181_groups.png", plot = last_plot())

#--- amyloid ---#

amrat <- merge(fujirebioabeta, ids, by="RID")
#remove duplicates by RID, keeping last value
amrat <- amrat[order(amrat$VISCODE ,decreasing=FALSE),]
amrat <- amrat %>% distinct(RID, .keep_all = TRUE)
am_ratio <- subset(amrat, select=c(RID,ABETA42_40))
names(am_ratio)[names(am_ratio) == 'RID'] <- 'participant'

allTest <- all
# merging together all the data
allTest <- merge(allTest, am_ratio, by="participant", all.x=TRUE)

ggplot(allTest, aes(x=group, y=ABETA42_40, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "Fujirebio Beta-Amyloid Ratio", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/fuji_amrat_groups.png", plot = last_plot())

#--- amrat euroimmun ---#

amrat2 <- merge(euroimmun, ids, by="RID")
#remove duplicates by RID, keeping last value
amrat2 <- amrat2[order(amrat2$VISCODE ,decreasing=FALSE),]
amrat2 <- amrat2 %>% distinct(RID, .keep_all = TRUE)
am_ratio2 <- subset(amrat2, select=c(RID,BETA.AMYLOID.42.40))
names(am_ratio2)[names(am_ratio2) == 'RID'] <- 'participant'

allTest <- all
# merging together all the data
allTest <- merge(allTest, am_ratio2, by="participant", all.x=TRUE)

ggplot(allTest, aes(x=group, y=BETA.AMYLOID.42.40, color=group)) +
  geom_boxplot() +
    theme_classic() +
    geom_point(aes(fill = DIAGNOSIS, color = group), size = 2, colour= "black", shape = 21, position = position_jitterdodge()) +
   labs(y= "EUROIMMUN CSF Beta-Amyloid Ratio", x="Group", legend.title = " Clinical Diagnosis") +
  scale_fill_manual(values=c("black", "darkgray", "white")) +
  scale_color_manual(values=c("#bc3c29", "#e18727", "#20854e", "#0072b5")) +
   theme(axis.text.y=element_text(size=14), axis.text.x=element_text(size=20), 
  axis.title.y=element_text(size=22),axis.title.x=element_blank(), legend.text=element_text(size=16),plot.title = element_text(size=20),legend.title=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", p.adjust.method = "fdr", hide.ns = FALSE)

ggsave("output/2_Comparing_subgroups/Figures/Supplement/euro_amrat_groups.png", plot = last_plot())

#--- network-visualization ---#

install.packages("igraph")
library(igraph)

install.packages("ggraph")
library(ggraph)

install.packages("tidygraph")
library(tidygraph)

install.packages("patchwork")
library(patchwork)

load("denseCoreClusterMatrix.rda")

ADNISNF_Groups <- as.data.frame(all$groups)
names(ADNISNF_Groups)[1] <- "groups"
ADNISNF_Diagnoses <- as.data.frame(all$DIAGNOSIS)
names(ADNISNF_Diagnoses)[1] <- "DIAGNOSIS"

set.seed(123) # For reproducible results

# Convert the correlation matrix into a graph object
# Note: We are using absolute values of correlations to consider the strength, regardless of the sign.
graph <- graph_from_adjacency_matrix(abs(dense), weighted=TRUE, mode="undirected", diag=FALSE)

# Plot the graph using a spring-loaded layout (force-directed layout)
# Adjusting the layout parameter can change the appearance. Here, we use layout_with_fr (Fruchterman-Reingold)
plot(graph, layout=layout_with_fr(graph), edge.width=E(graph)$weight*100, main="Spring-Loaded Graph from Correlation Matrix")

set.seed(101)
final_layout <- ggraph::create_layout(graph, 'igraph', algorithm = 'fr')

V(graph)$degree <- degree(graph)
first <- ggraph::ggraph(final_layout) + 
  ggraph::geom_edge_link(aes(width = weight), edge_alpha = 0.01) + 
  geom_point(aes(x = x, y = y, size = 5, fill = as.factor(ADNISNF_Groups$groups)), shape = 21, stroke = 1, color = "black") + 
  ggforce::theme_no_axes() + 
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#bc3c29", "#e18727", "#20854e", "#0072b5"))


V(graph)$degree <- degree(graph)
second <- ggraph::ggraph(final_layout) + 
    ggraph::geom_edge_link(aes(width = weight), edge_alpha = 0.01) + 
    geom_point(aes(x = x, y = y, size = 5, fill = as.factor(ADNISNF_Diagnoses$DIAGNOSIS)), shape = 21, stroke = 1, color = "black") +
    ggforce::theme_no_axes() + 
    theme_void() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("black", "darkgray", "white"))

library(patchwork)

first + second 
ggsave("output/2_Comparing_subgroups/Figures/Figure1/networkVis.png", plot = last_plot())

#save RDA
save(dense, file='denseCoreClusterMatrix.rda')
save(ADNISNF_Groups, file='ADNISNF_Groups.rda')
save(ADNISNF_Diagnoses, file='ADNISNF_Diagnoses.rda')
