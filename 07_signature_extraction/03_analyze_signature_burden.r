
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(ape)
library(ggtree)
library(reshape2)
library(ggplot2)
library(cowplot)


#tree_dir="/Users/so11/phd/psoriasis/scratch/phylogenics/MPBoot/snv/consensus_trees/"
#toutmeans <- read.table("/Users/so11/phd/psoriasis/scratch/snv/branch_exposures_w_prior.txt", h=T)
toutmeans <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/hdp/snv/branch_exposures_w_prior.txt")
tree_dir="/lustre/scratch119/humgen/projects/psoriasis/phylogenics/MPBoot/snv/consensus_trees/"

sample_meta <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/sample_info/sample_meta.txt", h=T, stringsAsFactors = F)
patients <- unique(sample_meta$patient_ID[sample_meta$Seq_status=="WES" & is.na(sample_meta$excl_criteria)])


## Hard-coded. May need to change if you re-run the HDP
sigNames <- c("Unassigned","SBS7b - UV","PUVA treatment", "SBS1/5", "UV-component N1",
              "UV-component N2","SBS2 - APOBEC","SBS7d - UV", "New component N3",
              "SBS13 - APOBEC")




all_sig_muts <- data.frame()
for(patient in patients) {
  
  tree <- read.tree(paste(tree_dir, patient, "/", patient, "_MutAssign_unadj.tree", sep=""))
  tree_original <- tree
  tree_df <- fortify(tree_original)
  
  ## Extract the branches for this patient. 
  b <- toutmeans[grep(paste(patient, "_", sep=""), rownames(toutmeans)),]
  b$node <- unlist(strsplit(rownames(b), split="_"))[c(FALSE, TRUE)]
  
  test <- merge(tree_df, b, by="node", all.x=T)
  
  ## For each signature, multiply the branch lengths with the fraction of mutations attributed to that signature
  ## in that branch. Then sum over all 'anchestor branches' to obtain total counts for each crypt. 
  for(i in 1:length(sigNames)) {
    tree_df <- fortify(tree_original)
    adj.branch <- test$branch.length
    
    adj.branch[!is.na(test[,(9+i)])] <- test$branch.length[!is.na(test[,(9+i)])]*test[!is.na(test[,(9+i)]), (9+i)]
    tree$edge.length=adj.branch[tree$edge[,2]]
    tree_df <- fortify(tree)
    if(i==1) {
      test_sigs <- tree_df[tree_df$isTip, c("label", "x")]
      
      ## Some branches have fewer than 50 mutations and so were not included in the signature extraction. 
      ## I want to find the total number of mutations belonging to such branches that precede each crypt
      ## and subtract that number from the number of mutations assigned to each signature. 
      no_extraction_branches <- adj.branch
      no_extraction_branches[!is.na(test[,(9+i)])] <- 0
      tree_tmp <- tree_original
      tree_tmp$edge.length=no_extraction_branches[tree_tmp$edge[,2]]
      tree_tmp_df <- fortify(tree_tmp)
      no_ext <- tree_tmp_df[tree_tmp_df$isTip, "x"]
    } else {
      test_sigs <- cbind(test_sigs, (tree_df[tree_df$isTip, c("x")]-no_ext))
    }
    #write.tree(tree, paste(tree_dir, patient, "/", patient, "_", sigNames[i], "_truncBinom_adj.tree", sep=""))
  }
  colnames(test_sigs) <- c("sample_ID", sigNames)
  all_sig_muts <- rbind(all_sig_muts, test_sigs)
  
}

sample_meta <- merge(sample_meta, all_sig_muts,by.x="True_name" ,by.y="sample_ID")
sample_meta$type <- "Lesional"
sample_meta$type[grep("H", sample_meta$True_name)] <- "Healthy"


## Plot histograms of the signature exposures
sbs7b <- ggplot(sample_meta, aes(x=SBS7b, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

unassigned <- ggplot(sample_meta, aes(x=Unassigned, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

puva <- ggplot(sample_meta, aes(x=PUVA, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

sbs7a <- ggplot(sample_meta, aes(x=SBS7a, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

new1 <- ggplot(sample_meta, aes(x=New1, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

new2 <- ggplot(sample_meta, aes(x=New2, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

sbs2 <- ggplot(sample_meta, aes(x=SBS2, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

sbs13 <- ggplot(sample_meta, aes(x=SBS13, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))


plot_grid(sbs7a, sbs7b, new1,new2, sbs2, sbs13, puva, unassigned, ncol=2)


## Not really a correlation between SBS7a and b.
ggplot(sample_meta, aes(x=SBS7a,y=SBS7b, fill=type)) + geom_point( colour="black", shape=21) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))


ggplot(sample_meta, aes(x=Unassigned,y=New1, fill=type)) + geom_point( colour="black", shape=21) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

ggplot(sample_meta, aes(x=New2,y=New1, fill=type)) + geom_point( colour="black", shape=21) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))



## Check the association of the mutation burden with age and disease
## duration after removing UV and PUVA associated mutations
sample_meta$burden_no_UV_PUVA <- sample_meta$Unassigned + sample_meta$New1 + sample_meta$New2 + sample_meta$SBS13 + sample_meta$SBS2

ggplot(sample_meta, aes(x=burden_no_UV_PUVA, fill=type)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

patient_meta <- read.table("/Users/so11/phd/psoriasis/sample_info/patient_meta.txt", h=T, stringsAsFactors = F)

df <- merge(sample_meta, patient_meta, by.x="patient_ID", by.y="Patient.ID")
df$Age_at_sampling <- df$Year_of_sampling - df$Year_of_birth
df$DiseaseDuration <- df$Year_of_sampling - df$Year_of_onset

ggplot(df, aes(x=Age_at_sampling,y=burden_no_UV_PUVA, fill=type)) + geom_point( colour="black", shape=21) + 
  scale_fill_manual(values=c("#FFDBAC","#e18377"))

df$Location <- "Torso"
df$Location[grep("_LL_", df$True_name)] <- "Lower_leg"
df$Location[grep("_UL_", df$True_name)] <- "Upper_leg"
df$Location[grep("_LA_", df$True_name)] <- "Lower_arm"
df$Location[grep("_UA_", df$True_name)] <- "Upper_arm"

df$Smoking <- factor(df$Smoking, levels=c("never_smoker", "ex_smoker", "smoker"))

library(nlme)

model_null.lme <- lme(fixed = burden_no_UV_PUVA ~ Age_at_sampling + Location + BMI + Smoking + PASI, 
                      random = list(patient_ID = pdSymm(form = ~ Age_at_sampling - 1), biopsy_ID = pdSymm(form = ~ Age_at_sampling - 1)), 
                      weights = varIdent(form= ~ 1 | type),
                      data = df, method="ML")


model_dur.lme <- lme(fixed = burden_no_UV_PUVA ~ Age_at_sampling + Location + BMI + Smoking + PASI + DiseaseDuration, 
                     random = list(patient_ID = pdSymm(form = ~ Age_at_sampling - 1), biopsy_ID = pdSymm(form = ~ Age_at_sampling - 1)), 
                     weights = varIdent(form= ~ 1 | type),
                     data = df, method="ML")

summary(model_null.lme)
summary(model_dur.lme)

anova(model_null.lme,model_dur.lme, test=T)$"p-value"[2]


SDs <- data.frame(patient_ID=unique(df$patient_ID), sd=tapply(df$burden_no_UV_PUVA, df$patient_ID, sd))
SDs <- merge(SDs, df, by="patient_ID")

ggplot(SDs, aes(x=Age_at_sampling,y=sd^2)) + geom_point() + 
  ggtitle("Variance after removing UV and PUVA")


df$burden_no_PUVA <- df$SBS13 + df$SBS2 + df$SBS7b + df$SBS7a + df$New1 + df$New2
SDs <- data.frame(patient_ID=unique(df$patient_ID), sd=tapply(df$burden_no_PUVA, df$patient_ID, sd))
SDs <- merge(SDs, df, by="patient_ID")

ggplot(SDs, aes(x=Age_at_sampling,y=sd^2)) + geom_point() + 
  ggtitle("Variance after removing PUVA only")
