setwd("./PREDOC_COURSE/QBM/")
library(tidyverse)
library(rcompanion)
library(ape)
library(vegan)
library(labdsv)
library(limma)
library(PRROC)
library(ggrepel)

# pkgconfig::set_config("tibble::rownames" = NA)

# colorectal <- read.delim("Colorectal_cancer/feat_all.tsv", sep="\t")
# colorectal_tbl <- as_tibble(colorectal, rownames=NA)
# colnames(colorectal_tbl) <- str_replace_all(colnames(colorectal_tbl), fixed("."), "-")
# colorectal_tbl <- rownames_to_column(colorectal_tbl, var="Bacteria")

colorectal <- read_tsv("Colorectal_cancer/feat_genus.tsv")
colorectal_species <- read_tsv("Colorectal_cancer/feat_all_0removed.tsv")
colorectal_meta <- read.delim("Colorectal_cancer/meta_all.tsv", sep="\t")

# Select only samples for which meta-data info is available
colorectal <- select(colorectal, c(Genus,colorectal_meta$Sample_ID))
# Normalise for library size
colorectal_norm <- bind_cols(colorectal["Genus"], colorectal[,2:ncol(colorectal)]/colSums(colorectal[,2:ncol(colorectal)]))
# Transform to log
colorectal_log <- bind_cols(colorectal["Genus"], log(colorectal[,2:ncol(colorectal)]))
colorectal_lognorm <- bind_cols(colorectal_log["Genus"], colorectal_log[,2:ncol(colorectal_log)]/colSums(colorectal_log[,2:ncol(colorectal_log)]))

ggplot(colorectal_meta) +
  geom_violin(aes(x="Age", y=Age, col=Group, fill=Group), size=0.9, alpha=0.1) +
  scale_color_manual(labels = c("Colorectal cancer", "Control"), values = c("#CC0000", "#3399FF")) +
  scale_fill_manual(labels = c("Colorectal cancer", "Control"), values = c("#CC0000", "#3399FF")) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=12), legend.title = element_text(face="bold"), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())


CRC <- filter(colorectal_meta, Group=="CRC")
mean(CRC$Age, na.rm=TRUE)
CTR <- filter(colorectal_meta, Group=="CTR")
mean(CTR$Age, na.rm=TRUE)

# set.seed(32)
# wi <- wilcoxonR(x=colorectal_meta[,"Age"], g=colorectal_meta[,"Group"], ci = TRUE, R=1000)


### Principal Coordinate Analysis (PCoA) ###
long_colorectal <- gather(colorectal_norm, key=Sample_ID, value=Relative_abundance, colnames(colorectal_norm)[2:ncol(colorectal_norm)]) %>% 
  mutate(Group = ifelse(Sample_ID %in% CRC$Sample_ID, 1, 0))

s <- as.data.frame(bind_cols(c(long_colorectal["Sample_ID"], long_colorectal["Genus"], long_colorectal["Relative_abundance"]))) 
color_matrix <- matrify(s)

# Construct a (dis)similarity matrix based on Bray-Curtis distance
dist <- vegdist(color_matrix,  method = "bray")
PCOA <- pcoa(dist)
## Some distance measures may result in negative eigenvalues. In that case, add a correction:
# PCOA <- pcoa(dist, correction = "cailliez")

# Barplot of eigenvalues
# barplot(PCOA$values$Relative_eig[1:10])
ggplot() + 
  geom_bar(aes(x=seq(1:10), y=PCOA$values$Eigenvalues[1:10]), stat = "identity") +
  labs(y="Varience explained", x="PCo", title="Principle Coordinate Analysis") +
  scale_x_continuous(breaks = seq(1:10), labels=as.character(seq(1:10))) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=12))


# Make sure colname order is the same: 
sum(rownames(color_matrix) == rownames(PCOA$vectors[,1:2])) == 575

color_matrix_ext <- tibble(`Sample_ID` = rownames(PCOA$vectors[,1:5]), 
                           `PCo1` = PCOA$vectors[,1], 
                           `PCo2` = PCOA$vectors[,2],
                           `PCo3` = PCOA$vectors[,3],
                           `PCo4` = PCOA$vectors[,4],
                           `PCo5` = PCOA$vectors[,5],
                           `PCo6` = PCOA$vectors[,6])
color_matrix_ext <- merge(color_matrix_ext, select(colorectal_meta, c("Sample_ID", "Age","Country", "Gender", "BMI", "Study", "Group", "Localization")), by="Sample_ID")

ggplot(color_matrix_ext) + 
  geom_point(aes(x=PCo1, y=PCo2, color=Study, shape=Group)) +
  labs(x="PCo_1", y="PCo_2", title="PCo1 vs. PCo2") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=12))

# biplot(PCOA) # package plot
# biplot.pcoa(PCOA, color_matrix) # package plot

# REMOVE BATCH EFFECT OF STUDY

# Add pseudocount of 10^-10 and take the ln
colorectal_lognorm <- bind_cols(colorectal_norm["Genus"], log(colorectal_norm[,2:ncol(colorectal_norm)] + 10^(-8)))

long_colorectal_lognorm <- gather(colorectal_lognorm, key=Sample_ID, value=Relative_abundance, colnames(colorectal_lognorm)[2:ncol(colorectal_lognorm)]) %>% 
  mutate(Group = ifelse(Sample_ID %in% CRC$Sample_ID, 1, 0))  

color_matrix_log <- matrify(as.data.frame(bind_cols(c(long_colorectal_lognorm["Sample_ID"], long_colorectal_lognorm["Genus"], long_colorectal_lognorm["Relative_abundance"]))))

sum(rownames(color_matrix_log) == colorectal_meta$Sample_ID) == 575
study_rm <- removeBatchEffect(t(color_matrix_log),
                              batch = colorectal_meta$Study, 
                              covariates = matrix(c(colorectal_meta$Age), ncol = 1))

# save(study_rm, file = "Count_matrix_wo_Study_effect.Rdata")

dist <- vegdist(exp(t(study_rm)),  method = "bray")
PCOA <- pcoa(dist)

color_matrix_ext <- tibble(`Sample_ID` = rownames(PCOA$vectors[,1:5]), 
                           `PCo1` = PCOA$vectors[,1], 
                           `PCo2` = PCOA$vectors[,2],
                           `PCo3` = PCOA$vectors[,3],
                           `PCo4` = PCOA$vectors[,4],
                           `PCo5` = PCOA$vectors[,5],
                           `PCo6` = PCOA$vectors[,6])
color_matrix_ext <- merge(color_matrix_ext, select(colorectal_meta, c("Sample_ID", "Age", "Gender", "BMI", "Study", "Group", "Country", "Localization", "Sampling_rel_to_colonoscopy")), by="Sample_ID")

ggplot(color_matrix_ext) + 
  geom_point(aes(x=PCo1, y=PCo2, color=Study, shape=Group)) +
  labs(x="PCo_1", y="PCo_2", title="PCo1 vs. PCo2") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=12))


### After removing batch-effect
colorectal_wo_bs <- as_tibble(exp(study_rm))
colorectal_wo_bs <- bind_cols(tibble(`Genus`=rownames(study_rm)), colorectal_wo_bs)
long_colorectal_wo_bs <- gather(colorectal_wo_bs, key=Sample_ID, value=Relative_abundance, colnames(colorectal_wo_bs)[2:ncol(colorectal_wo_bs)]) %>% 
  mutate(Group = ifelse(Sample_ID %in% CRC$Sample_ID, 1, 0)) 
color_matrix_wo_bs <- matrify(as.data.frame(bind_cols(c(long_colorectal_wo_bs["Group"], long_colorectal_wo_bs["Genus"], long_colorectal_wo_bs["Relative_abundance"]))))


#### GENERELAZATION PERFORMANCE PLOTS ####

## ROC curve ##
rf_roc <- read_csv("./rf_roc.txt")
rf_roc <- select(rf_roc, -threshold) %>% mutate(Model = rep("RF", nrow(rf_roc)), Taxon=rep("Genus", nrow(rf_roc)))
logreg_roc <- read_csv("./logr_roc.txt")
logreg_roc <- select(logreg_roc, -threshold) %>% mutate(Model = rep("LogReg", nrow(logreg_roc)), Taxon=rep("Genus", nrow(logreg_roc)))

rf_roc_species <- read_csv("roc_species/rf_roc.txt")
rf_roc_species <- select(rf_roc_species, -threshold) %>% mutate(Model = rep("RF", nrow(rf_roc_species)), Taxon=rep("Species", nrow(rf_roc_species)))
logreg_roc_species <- read_csv("roc_species/logr_roc.txt")
logreg_roc_species <- select(logreg_roc_species, -threshold) %>% mutate(Model = rep("LogReg", nrow(logreg_roc_species)), Taxon=rep("Species", nrow(logreg_roc_species)))

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

simple_auc(logreg_roc$tpr, logreg_roc$fpr)
simple_auc(rf_roc$tpr, rf_roc$fpr)
simple_auc(logreg_roc_species$tpr, logreg_roc_species$fpr)
simple_auc(rf_roc_species$tpr, rf_roc_species$fpr)

ROC <- bind_rows(bind_rows(bind_rows(rf_roc, logreg_roc), rf_roc_species), logreg_roc_species)

ggplot(ROC) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed", colour="grey74") + 
  geom_line(aes(x=fpr, y=tpr, color=Model, linetype=Taxon)) +
  ylim(0, 1) +
  labs(x="FPR", y="TPR") +
  scale_color_manual(labels = c("LogReg", "RF"), values = c("hotpink3", "#3399FF")) +
  geom_label_repel(data=data.frame(x=0.258, y=0.825),aes(x=x, y=y), label="AUC = 0.83", color="hotpink3", segment.color="hotpink3", nudge_y = 0.1, nudge_x = -0.12) +
  geom_label_repel(data=data.frame(x=0.34, y=0.842),aes(x=x, y=y), label="AUC = 0.84", color="#3399FF", segment.color="#3399FF", nudge_y = 0.15)+
  geom_label_repel(data=data.frame(x=0.75, y=0.924),aes(x=x, y=y), label="AUC = 0.8", color="hotpink3", segment.color="hotpink3", nudge_y = -0.1, segment.size = 0.3, nudge_x = 0.1) +
  geom_label_repel(data=data.frame(x=0.5, y=0.86),aes(x=x, y=y), label="AUC = 0.79", color="#3399FF", segment.color="#3399FF", segment.size	= 0.3, nudge_y = -0.1) +
  theme_classic() +
  theme(title = element_text(face="bold", size=12), axis.text = element_text(size = 10))


## PR curve ##
rf_pr <- read_csv("./precision_recall_genus/rf_precision_recall.txt")
rf_pr <- select(rf_pr, -threshold) %>% mutate(Model = rep("RF", nrow(rf_pr)), Taxon=rep("Genus", nrow(rf_pr)))
logreg_pr <- read_csv("./precision_recall_genus/logr_precision_recall.txt")
logreg_pr <- select(logreg_pr, -threshold) %>% mutate(Model = rep("LogReg", nrow(logreg_pr)), Taxon=rep("Genus", nrow(logreg_pr)))

rf_pr_species <- read_csv("precision_recall_species/rf_precision_recall.txt")
rf_pr_species <- select(rf_pr_species, -threshold) %>% mutate(Model = rep("RF", nrow(rf_pr_species)), Taxon=rep("Species", nrow(rf_pr_species)))
logreg_pr_species <- read_csv("precision_recall_species/logr_precision_recall.txt")
logreg_pr_species <- select(logreg_pr_species, -threshold) %>% mutate(Model = rep("LogReg", nrow(logreg_pr_species)), Taxon=rep("Species", nrow(logreg_pr_species)))

PR <- bind_rows(bind_rows(bind_rows(rf_pr, logreg_pr), rf_pr_species), logreg_pr_species)

simple_auc(logreg_pr$recall, logreg_pr$precision)
simple_auc(rf_pr$recall, rf_pr$precision)
simple_auc(logreg_pr_species$recall, logreg_pr_species$precision)
simple_auc(rf_pr_species$recall, rf_pr_species$precision)


ggplot(PR) +
  geom_abline(slope=-1, intercept = 1, linetype = "dashed", colour="grey74") +
  #geom_line(data = data.frame(x=c(0,1), y=c(1,0)), aes(x=x, y=y), linetype = "dashed", colour="grey74") +
  geom_line(aes(x=precision, y=recall, color=Model, linetype=Taxon)) +
  ylim(0, 1) +
  xlim(0,1) +
  labs(x="Precision", y="Recall") +
  scale_color_manual(labels = c("LogReg", "RF"), values = c("hotpink3", "#3399FF")) +
  geom_label_repel(data=data.frame(x=0.7, y=0.88),aes(x=x, y=y), label="AUC = 0.34", color="hotpink3", segment.color="hotpink3", nudge_y = 0.15, nudge_x = 0.12) +
  geom_label_repel(data=data.frame(x=0.7, y=0.83),aes(x=x, y=y), label="AUC = 0.32", color="#3399FF", segment.color="#3399FF", nudge_y = 0.12, nudge_x = 0.15)+
  geom_label_repel(data=data.frame(x=0.95, y=0.4),aes(x=x, y=y), label="AUC = 0.35", color="hotpink3", segment.color="hotpink3", nudge_y = 0.08, nudge_x = 0.8, segment.size = 0.3) +
  geom_label_repel(data=data.frame(x=0.8, y=0.3),aes(x=x, y=y), label="AUC = 0.24", color="#3399FF", segment.color="#3399FF", segment.size	= 0.3, nudge_x = -0.3, nudge_y=0.1) +
  theme_classic() +
  theme(title = element_text(face="bold", size=12), axis.text = element_text(size = 10))

