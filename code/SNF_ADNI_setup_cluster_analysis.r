### Setting up data for the analysis and comparison of SNF clusters

clusters <- read.csv("output/1_SNF_analysis/4clust_groups_k18_0.8_1000perms.csv")

clusters$X <- NULL
clusters$id <- NULL
names(clusters)[names(clusters)=="V1"] <- "ID"

# Importing individual data type matrices with measure headers
directory <- ("output/1_SNF_analysis")
TA <-read.csv(file.path(directory, paste("TA.csv", sep="")))
CV <-read.csv(file.path(directory, paste("CV.csv", sep="")))
SA <-read.csv(file.path(directory, paste("SA.csv", sep="")))
behav <-read.csv(file.path(directory, paste("Behavioural.csv", sep="")))
ids <-read.csv(file.path(directory, paste("IDs.csv", sep="")))
names(ids)[names(ids)=="x"] <- "RID"

#--- Retrieving patient demo info for those included in our analysis ---#

demo <- merge(ptdemog, ids, by="RID")

# Remove duplicates by RID, keeping last value
demo <- demo[order(demo$VISCODE ,decreasing=FALSE),]
demo <- demo %>% distinct(RID, .keep_all = TRUE)

# Select date of birth, education, gender, marital status, recent job
demo2 <- subset(demo, select=c(RID, AGE, PTEDUCAT, PTGENDER, PTMARRY, PTWRECNT))
names(demo2)[names(demo2) == 'RID'] <- 'participant'

#--- Retrieving patient diagnostic summary for those included in our analysis ---#

dx <- dxsum %>% filter(COLPROT == "ADNI3")
demoDX <- merge(dx, ids, by="RID") # includes multiple entries for same participants

# Remove duplicates by RID, keeping last value
demoDX <- demoDX[order(demoDX$VISCODE ,decreasing=FALSE),]
demoDX <- demoDX %>% distinct(RID, .keep_all = TRUE)

# Merge and reduce dfs
ridvis <-read.csv(file.path(directory, paste("ridvis.csv", sep="")))
demoDX2<- Reduce(function(x, y) merge(x, y, by = c("RID", "VISCODE"), all.x=TRUE), list(ridvis, dx))
dx2 <- subset(demoDX2, select=c(RID, DIAGNOSIS, DXMDES))
names(dx2)[names(dx2) == 'RID'] <- 'participant'

#--- Retrieving patient demo info for those included in our analysis ---#

# Depression
gdep <- merge(gdscale, ids, by="RID")
#remove duplicates by RID, keeping last value
gdep <- gdep[order(gdep$VISCODE ,decreasing=FALSE),]
gdep <- gdep %>% distinct(RID, .keep_all = TRUE)
depression <- subset(gdep, select=c(RID,GDTOTAL))
names(depression)[names(depression) == 'RID'] <- 'participant'

# Hachinski score
ghach <- merge(modhach, ids, by="RID")
#remove duplicates by RID, keeping last value
ghach <- ghach[order(ghach$VISCODE ,decreasing=FALSE),]
ghach <- ghach %>% distinct(RID, .keep_all = TRUE)
hachinski <- subset(ghach, select=c(RID,HMSCORE))
names(hachinski)[names(hachinski) == 'RID'] <- 'participant'

# Smoking
gsmoke <- merge(medhist, ids, by="RID")
#remove duplicates by RID, keeping last value
gsmoke <- gsmoke[order(gsmoke$VISCODE ,decreasing=FALSE),]
gsmoke <- gsmoke %>% distinct(RID, .keep_all = TRUE)
smoking <- subset(gsmoke, select=c(RID,MH16SMOK))
names(smoking)[names(smoking) == 'RID'] <- 'participant'

# Alcohol
galc <- merge(medhist, ids, by="RID")
#remove duplicates by RID, keeping last value
galc <- galc[order(galc$VISCODE ,decreasing=FALSE),]
galc <- galc %>% distinct(RID, .keep_all = TRUE)
alcohol <- subset(galc, select=c(RID,MH14ALCH))
names(alcohol)[names(alcohol) == 'RID'] <- 'participant'

# Cardiovascular
gcard <- merge(medhist, ids, by="RID")
#remove duplicates by RID, keeping last value
gcard <- gcard[order(gcard$VISCODE ,decreasing=FALSE),]
gcard <- gcard %>% distinct(RID, .keep_all = TRUE)
cardioVasc <- subset(gcard, select=c(RID,MH4CARD))
names(cardioVasc)[names(cardioVasc) == 'RID'] <- 'participant'

#--- Merging together all the data ---#
all <- merge(clusters, dx2, by="participant", all.x=TRUE)
all <- merge(all, demo2, by="participant", all.x=TRUE)
all <- merge(all, depression, by="participant", all.x=TRUE)
all <- merge(all, hachinski, by="participant", all.x=TRUE)
all <- merge(all, smoking, by="participant", all.x=TRUE)
all <- merge(all, alcohol, by="participant", all.x=TRUE)
all <- merge(all, cardioVasc, by="participant", all.x=TRUE)

all <- cbind(all, TA)
all <- subset(all, select = -X)
all <- cbind(all, CV)
all <- subset(all, select = -X)
all <- cbind(all, SA)
all <- subset(all, select = -X)
all <- cbind(all, behav)
all <- subset(all, select = -X)

# Ensuring classes of variables are correct
all$group <- as.factor(all$groups)
all$DXMDES <- as.factor(all$DXMDES)
all$PTGENDER <- as.factor(all$PTGENDER)
all$PTMARRY <- as.factor(all$PTMARRY)
all$PTWRECNT <- as.factor(all$PTWRECNT)

all$AGE <- as.numeric(all$AGE)
all$PTEDUCAT <- as.numeric(all$PTEDUCAT)

clust1 <- all[which(all$group=="1"),]
clust2 <- all[which(all$group=="2"),]
clust3 <- all[which(all$group=="3"),]
clust4 <- all[which(all$group=="4"),]

save(all, file='ADNISNF_all.rda')

orgMeasure <-read.csv(file.path(directory, paste("NMI_Names.csv", sep="")))
orgMeasure <- orgMeasure[order(orgMeasure$NMI,decreasing=TRUE),]

basic_demos <- all[,c("participant", "DIAGNOSIS",  "groups", "PTGENDER", "MOCA", "TOTAL13", "CDRSB", "MMSCORE", "GDTOTAL")]
TA_stacked <- cbind(basic_demos, TA)
CV_stacked <- cbind(basic_demos, CV)
SA_stacked <- cbind(basic_demos, SA)

TA_stacked <- all %>%
  dplyr::select(participant, DIAGNOSIS, groups, AGE, PTGENDER, ends_with("TA")) %>%
  gather(Region,thickness, -participant, -AGE, -DIAGNOSIS,-groups, -PTGENDER)

CV_stacked <- all %>%
  dplyr::select(participant, DIAGNOSIS, groups, AGE, PTGENDER, ends_with("CV")) %>%
  gather(Region,volume, -participant, -AGE, -DIAGNOSIS,-groups, -PTGENDER)

SA_stacked <- all %>%
  dplyr::select(participant, DIAGNOSIS, groups, AGE, PTGENDER, ends_with("SA")) %>%
  gather(Region,area, -participant, -AGE, -DIAGNOSIS,-groups, -PTGENDER)

all_features <- all %>%
  dplyr::select(participant, DIAGNOSIS, groups, AGE, PTGENDER, TOTAL13, MMSCORE, MOCA, CDRSB, GDTOTAL, ends_with("TA"), ends_with("CV"), ends_with("SA")) %>% 
   gather(Measure, brain_features, -participant, -AGE, -DIAGNOSIS, -groups, -PTGENDER, -GDTOTAL)

symps <- all
symps$MOCA <- (symps$MOCA - mean(symps$MOCA, na.rm = TRUE))/sd(symps$MOCA, na.rm = TRUE)
symps$TOTAL13 <- (symps$TOTAL13 - mean(symps$TOTAL13, na.rm = TRUE))/sd(symps$TOTAL13, na.rm = TRUE)
symps$CDRSB <- (symps$CDRSB - mean(symps$CDRSB, na.rm = TRUE))/sd(symps$CDRSB, na.rm = TRUE)
symps$MMSCORE <- (symps$MMSCORE - mean(symps$MMSCORE, na.rm = TRUE))/sd(symps$MMSCORE, na.rm = TRUE)

symptoms_stacked <- symps %>%
   dplyr::select(participant, DIAGNOSIS, groups, AGE, PTGENDER, TOTAL13, MMSCORE, MOCA, CDRSB, GDTOTAL) %>%
  gather(Symptom, severity, -participant, -AGE, -DIAGNOSIS, -groups, -PTGENDER, -GDTOTAL)

# Participant group and dx
groupDX <- all %>%
  dplyr::select(participant, DIAGNOSIS, groups)
directory <- ("output/1_SNF_analysis")
write.csv(groupDX, file=file.path(directory, paste("/groupDX.csv", sep="")))
