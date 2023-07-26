
## The purpose of this script is to create the figures for the manuscript.
working_dir="/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/"

microd_meta <- read.table(paste(working_dir, "Microdissection_metaData.txt", sep=""), h=T)
patient_meta <- read.table(paste(working_dir, "Patient_metaData.txt", sep=""), h=T)
treatment_data <- read.table(paste(working_dir, "treatment_data_combined.txt", sep=""), h=T)

patient_meta$Shows_Psoralen_signature <- ifelse(patient_meta$anyCloneGT50PUVA=="Yes" | patient_meta$Patient.ID %in% c("patient89","patient56","patient105"), T, F)
patient_meta <- merge(patient_meta, treatment_data, by.x="Patient.ID", by.y="PatientID", all.x=T)
patient_meta$Sex[patient_meta$Sex=="female"] <- "Female"
patient_meta$Sex[patient_meta$Sex=="male"] <- "Male"


library(ggplot2)
library(reshape2)
library(cowplot)
library(ggsignif)
library(ggExtra)
library(scales)

## DEFINE PLOTTING VARIABLES
#############################
font_family <- "Helvetica"
text_size <- 7
pdf_width_mm <- 180
pdf_width_in <- pdf_width_mm / 25.4
pdf_height_mm <- 185
pdf_height_in <- pdf_height_mm / 25.4
BASESIZE=18

# Location colour vector
#Abdomen     Arm    Back   Flank     Leg
loc_colors <- c("#264653", "#2A9D8F","#E9C46A","#F4A261","#E76F51")

# Disease type (lesional vs non-lesional) colour vector
type_colours <- c("#FF7075", "#5DB4EA")


## Plot patient demographics

ggplot(patient_meta, aes(x=Age_at_sampling, y=Disease_duration, colour=Shows_Psoralen_signature)) + geom_point() + 
  theme_classic(base_size = BASESIZE) + scale_color_brewer(palette="Set2")

ggplot(patient_meta, aes(x=Age_at_sampling, y=Disease_duration, colour=EverPuva)) + geom_point() + 
  theme_classic(base_size = BASESIZE) + scale_color_brewer(palette="Set2")

tapply(patient_meta$Age_at_sampling, patient_meta$Shows_Psoralen_signature, median, na.rm=T)



puva <- patient_meta[,c("Patient.ID", "AmountPuva")]
puva$Treatment <- "PUVA"
colnames(puva)[2] <- "Amount"
uvb <- patient_meta[,c("Patient.ID", "AmountUVB")]
uvb$Treatment <- "UVB"
colnames(uvb)[2] <- "Amount"
tmp <- rbind(puva, uvb)
tmp$Amount[tmp$Amount=="unknown"] <- "Has history but \n dose unknown"
tmp$Amount <- factor(tmp$Amount, levels=c("None", "<=50", "51-200", ">200","Has history but \n dose unknown"))

ggplot(tmp, aes(x=Amount)) + geom_bar(fill="gray") + theme_bw(base_size = BASESIZE) +
  facet_grid(~Treatment) + labs(x="Cycles of phototherapy", y="# Patients") + 
  theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(expand = c(0,0), limits = c(0,111))


##############
##
## FIGURE 1
###############

## This figure is mostly made in Inkscape. Just want the VAF histogram
## but I'm also planning to print out some of the info that will be manually added
## to Inkscape

#Print info
biopsy_meta <- unique(microd_meta[microd_meta$ExclusionCriteria=="PASS",c("BiopsyID", "MetaLocation", "DiseaseStatus")])
table(biopsy_meta$MetaLocation, biopsy_meta$DiseaseStatus)
table(sample_meta$Location.y[sample_meta$Seq_status=="WES"], sample_meta$type[sample_meta$Seq_status=="WES"])
table(patient_meta$Sex)


## Variant allele frequencies:
ggplot(microd_meta, aes(x=MedianVAF, fill=DiseaseStatus)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + theme_classic(base_size = BASESIZE) +
  scale_fill_manual(values=type_colours) + labs(x="Median VAF of microbiopsy",y="# Microbiopsies") +
  theme(legend.title = element_blank(), legend.position = "top") +
  geom_vline(xintercept=median(microd_meta$MedianVAF[microd_meta$DiseaseStatus=="Lesional"],na.rm=T), size=1,colour="#E00007", linetype="dashed") + 
  geom_vline(xintercept=median(microd_meta$MedianVAF[microd_meta$DiseaseStatus=="Non-lesional"],na.rm=T), size=1,colour="#156CA2", linetype="dashed")


write.table(microd_meta[,c("SampleID","MedianVAF","DiseaseStatus")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig1c.txt", sep=""), quote=F, row.names = F, sep = "\t")

pdf(paste(working_dir, "Olafsson_Fig1c.pdf", sep=""), height=2.7, width = (pdf_width_in/2)-0.5, fonts = font_family)
Main_fig1c <- ggplot(microd_meta, aes(x=MedianVAF, fill=DiseaseStatus)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + theme_classic()+  
  theme(text=element_text(size=text_size, family = font_family)) +
  scale_fill_manual(values=type_colours) + labs(x="Median VAF of microbiopsy",y="# Microbiopsies") +
  theme(legend.title = element_blank(), legend.position = "top") +
  facet_wrap(~DiseaseStatus, nrow=2, scales = "free_y") + theme(strip.background = element_blank()) +
  theme(legend.position = "none")
print(Main_fig1c)
dev.off()

tapply(microd_meta$MedianVAF, microd_meta$DiseaseStatus, median, na.rm=T)


shared_by_distance <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/shared_by_distance.txt", h=T)
shared_by_distance$type <- "Lesional"
shared_by_distance$type[grep("H", shared_by_distance$pairs)] <- "Non-lesional"

write.table(shared_by_distance[,c("pairs", "frac_shared","distances")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig1d.txt", sep=""), quote=F, row.names = F, sep = "\t")

pdf(paste(working_dir, "Olafsson_Fig1d.pdf", sep=""), height=2.7, width = (pdf_width_in/2)-0.5, fonts = font_family)
Main_fig1d <-ggplot(shared_by_distance, aes(x=distances, y=frac_shared, colour=type)) + 
  geom_point(size=0.5) + theme_classic()+  
  theme(text=element_text(size=text_size, family = font_family)) +
  scale_color_manual(values=type_colours) + labs(x="Micrometers separating microbiopsies",y="Fraction of shared mutations") +
  theme(legend.title = element_blank(), legend.position = "top") +
  facet_wrap(~type, nrow=1) + theme(strip.background = element_blank()) +
  theme(legend.position = "none")
print(Main_fig1d)
dev.off()


col_vec <- c("white", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")


p81 <- read.tree("/lustre/scratch126/humgen/projects/psoriasis/phylogenics/hdp_clustering/trees/patient81.cluster_tree.phylo")

p81$node.label[1] <- ""


pdf(paste(working_dir, "Olafsson_Fig1e.pdf", sep=""), height=2.36, width = (pdf_width_in/2)-0.5, fonts = font_family)
p <- ggtree(p81) + geom_nodepoint(aes(fill=label), shape=21, size=3) + 
  geom_tippoint(aes(fill=label), shape=21, size=3) + theme_classic()+  
  theme(text=element_text(size=text_size, family = font_family)) +
  theme(axis.line.x=element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "none") + labs(x="Number of mutations") + 
  coord_flip() + scale_x_reverse() + 
  scale_fill_manual(values=col_vec)

## Change the order of the branches to match the biopsy
p2 <- flip(p, 14,11) %>% flip(8, 12) %>% flip(7, 5) %>% flip(13,11) %>% flip(13,12) %>% flip(3,4)
print(p2)
dev.off()

###################################################
##
## SUPPLEMENTARY FIGURE 1 - related to Figure 1
###################################################

write.table(patient_meta[,c("Patient.ID", "Age_at_sampling","Disease_duration","Sex")], file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1a.txt", sep=""), quote=F, row.names = F, sep = "\t")

## Patient demographics:
p <- ggplot(patient_meta, aes(x = Age_at_sampling, y=Disease_duration, colour=Sex)) +
  geom_point() + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + scale_color_brewer(palette="Set2") + 
  theme(legend.position = c(0.2,0.8)) + labs(x="Age",y="Disease Duration", colour="")

# Changing the relative size
p1 <- ggMarginal(p, type = "densigram", 
           size = 3, xparams = list(fill = "#8DA0CB"),
           yparams = list(fill = "#A6D854"))


## Plot the treatment information
######################
puva <- patient_meta[,c("Patient.ID", "AmountPuva", "EverPuva")]
puva$Treatment <- "PUVA"
colnames(puva)[2] <- "NrSessions"
colnames(puva)[3] <- "Treated"
uvb <- patient_meta[,c("Patient.ID", "AmountUVB", "EverUVB")]
uvb$Treatment <- "UVB"
colnames(uvb)[2] <- "NrSessions"
colnames(uvb)[3] <- "Treated"

methotrexate <- patient_meta[, c("Patient.ID","EverMethotrexate")]
colnames(methotrexate)[2] <- "Treated"
methotrexate$Treatment <- "Methotrexate"

steroids <- patient_meta[, c("Patient.ID","EverSteroids")]
colnames(steroids)[2] <- "Treated"
steroids$Treatment <- "Topical Steroids"

photo <- rbind(puva, uvb)
photo$NrSessions[photo$NrSessions=="unknown"] <- "Unk"
photo$NrSessions <- factor(photo$NrSessions, levels=c("None","Unk", "<=50", "51-200", ">200"))

treatments <- rbind(puva[,c(1,3,4)], uvb[,c(1,3,4)],methotrexate, steroids)

treatments$Treated <- ifelse(treatments$Treated, "Treatment history", "No history")
treatments$Treated <- factor(treatments$Treated, levels=c("Treatment history", "No history"))
treatments$Treatment[treatments$Treatment=="Topical Steroids"] <- "Topical \n Steroids"

write.table(treatments[,c("Patient.ID", "Treated","Treatment")], file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1b.txt", sep=""), quote=F, row.names = F, sep = "\t")

p2 <- ggplot(treatments, aes(x=Treatment,fill=Treated)) + geom_bar() + theme_classic()+  
  theme(text=element_text(size=text_size, family = font_family)) +
  labs(x="", y="# Patients", fill="") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position = "top") + scale_fill_brewer(palette = "Paired") + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) + guides(fill = guide_legend(nrow = 1))

write.table(photo, file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1c.txt", sep=""), quote=F, row.names = F, sep = "\t")

p3 <- ggplot(photo, aes(x=NrSessions, fill=NrSessions)) + geom_bar() +theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + 
  facet_wrap(~Treatment) + labs(x="# Treatment sessions", y="# Patients") + 
  theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(expand = c(0,0), limits = c(0,111)) + 
  scale_fill_brewer(palette = "Set2") + theme(legend.position = "none")

plot_grid(pA, pB, pC, nrow = 1, rel_widths = c(3,2,2))




sensitivity <- read.table(paste(working_dir, "Sensitivity_technical_duplicates.txt", sep=""), h=T)

## The Volume of the microbiopsies
write.table(microd_meta[,c("SampleID", "CutVolume","DiseaseStatus")], file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1d.txt", sep=""), quote=F, row.names = F, sep = "\t")

pA <- ggplot(microd_meta, aes(x=CutVolume/1000000, fill=DiseaseStatus)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) +
  scale_fill_manual(values=type_colours) + labs(x=expression(paste("Microbiopsy Volume [", mm^3,"]", sep="")),y="# Microbiopsies") +
  theme(legend.title = element_blank(), legend.position = "top") +
  facet_wrap(~DiseaseStatus, nrow=2, scales = "free_y") + theme(strip.background = element_blank()) +
  theme(legend.position = "none")


## The Area of the microbiopsies
write.table(microd_meta[,c("SampleID", "CutArea","DiseaseStatus")], file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1e.txt", sep=""), quote=F, row.names = F, sep = "\t")

pB <- ggplot(microd_meta, aes(x=CutArea/1000000, fill=DiseaseStatus)) + geom_histogram(position="dodge") + 
  scale_y_continuous(expand = c(0, 0)) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) +
  scale_fill_manual(values=type_colours) + labs(x=expression(paste("Microbiopsy Surface Area [", mm^2,"]", sep="")),y="# Microbiopsies") +
  theme(legend.title = element_blank(), legend.position = "top") + xlim(c(0,0.02)) + 
  facet_wrap(~DiseaseStatus, nrow=2, scales = "free_y") + theme(strip.background = element_blank()) +
  theme(legend.position = "none")

## VAF as a function of Volume
pC <- ggplot(microd_meta, aes(y=MedianVAF, x=CutVolume/1000000, fill=DiseaseStatus)) + 
  geom_point(colour="black", shape=21) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) +
  scale_fill_manual(values=type_colours) + labs(y="Median VAF of sample",x=expression(paste("Microbiopsy Volume [", mm^3,"]", sep=""))) +
  theme(legend.title = element_blank(), legend.position = "top") +
  facet_wrap(~DiseaseStatus, nrow=2) + theme(strip.background = element_blank()) +
  theme(legend.position = "none")

## VAF as a function of Area
pD <- ggplot(microd_meta, aes(y=MedianVAF, x=CutArea/1000000, fill=DiseaseStatus)) + 
  geom_point(colour="black", shape=21) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) +
  scale_fill_manual(values=type_colours) + labs(y="Median VAF of sample",x=expression(paste("Microbiopsy Surface Area [", mm^2,"]", sep=""))) +
  theme(legend.title = element_blank(), legend.position = "top") + xlim(c(0,0.015))
## The Coverage distribution

write.table(microd_meta[,c("SampleID", "MedianCoverage","DiseaseStatus")], file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1f.txt", sep=""), quote=F, row.names = F, sep = "\t")

pE <- ggplot(microd_meta, aes(x=MedianCoverage, fill=DiseaseStatus)) + geom_histogram(position="dodge") +
  labs(x="Median On-target coverage", y="# Microbiopsies") + scale_y_continuous(expand = c(0, 0)) + 
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) +
  scale_fill_manual(values=type_colours) + 
  theme(legend.title = element_blank(), legend.position = "top")  + 
  facet_wrap(~DiseaseStatus, nrow=2, scales="free_y") + theme(strip.background = element_blank()) +
  theme(legend.position = "none")

## Sensitivity of mutation calls
write.table(sensitivity, file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig1g.txt", sep=""), quote=F, row.names = F, sep = "\t")

pF <- ggplot(sensitivity, aes(x=lower,y=S)) + geom_point() + ylim(0,1) + 
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + labs(x="Lower coverage of pair", y="Sensitivity") + 
  theme(legend.position = "top") + geom_vline(xintercept=median(microd_meta$MedianCoverage[microd_meta$ExclusionCriteria=="PASS"]),size=1, linetype="dashed") 

pdf(paste(working_dir, "Olafsson_ED_Fig1.pdf", sep=""),paper="a4", height = pdf_height_in, width = pdf_width_in-1.5, fonts = font_family)
ED_fig1 <- plot_grid(p1, plot_grid(p2,p3, pA, pB, pE, pF, ncol=2, labels=c("b","c","d", "e", "f","g")), nrow=2, rel_heights = c(1:4),labels = c("a",""))
print(ED_fig1)
dev.off()




############################################
##
##  FIGURE 2 - MUTATION BURDEN AND SIGNATURES
############################################

## Mutation burden as a function of age.
clone_burden <- read.table(paste(working_dir, "clone_mutation_burden.txt", sep=""), h=T)
clone_burden <- merge(clone_burden, microd_meta[,c("SampleID", "BiopsyID", "MetaLocation",
                                                   "PatientID", "DiseaseStatus")], by.x="HighCellFrac_sample", by.y="SampleID")
clone_burden <- merge(clone_burden, patient_meta[,c("Patient.ID", "Age_at_sampling", "Disease_duration",
                                                 "Sex", "BMI", "Smoking", "PASI")], by.x="PatientID", by.y="Patient.ID")

clone_burden$Disease_duration[clone_burden$DiseaseStatus=="Non-lesional"] <- 0
clone_burden_original <- clone_burden
## Remove some patients where the meta-data is uncertain
clone_burden <- clone_burden[!(clone_burden$PatientID %in% c("patient09", "patient44", "patient13", "patient50", "patient58")),]


library(nlme)
model_null.lme <- lme(fixed = TotalSBS_adj ~ Age_at_sampling  + MetaLocation, 
                      random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                      weights = varIdent(form= ~ 1 | DiseaseStatus),
                      data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

model_dur.lme <- lme(fixed = TotalSBS_adj ~ Age_at_sampling  + MetaLocation + Disease_duration, 
                      random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                      weights = varIdent(form= ~ 1 | DiseaseStatus),
                      data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

summary(model_null.lme)
summary(model_dur.lme)
anova(model_null.lme,model_dur.lme, test=T)$"p-value"[2]
lme.ints <- intervals(model_null.lme, which="fixed")$fixed

maxAge=max(clone_burden$Age_at_sampling)
ageEff=lme.ints["Age_at_sampling", "est."]
low <-lme.ints["Age_at_sampling", "lower"] 
upp <- lme.ints["Age_at_sampling", "upper"] 


pB <- ggplot(clone_burden, aes(y=TotalSBS_adj, x=Age_at_sampling, fill=DiseaseStatus)) + geom_point( colour="black", shape=21, size=0.5) + 
  scale_fill_manual(values=type_colours)  + 
  labs(y="Total SBS", x="Age at sample donation") +
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + theme(legend.title = element_blank(), legend.position = "top") +
  geom_ribbon(aes(ymin=Age_at_sampling*low, ymax=Age_at_sampling*upp, x=Age_at_sampling), alpha = 0.3, show.legend=F) +
  geom_line(aes(y=Age_at_sampling*ageEff, x=Age_at_sampling)) + 
  guides(fill = guide_legend(override.aes = list(size=2.5)))

## Plot the burden of non-PUVA related mutations
clone_burden$noPUVA <- clone_burden$TotalSBS_adj - clone_burden$Psoralens

model_noPUVA.null <- lme(fixed = noPUVA ~ Age_at_sampling  + MetaLocation + Disease_duration, 
                         random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                         weights = varIdent(form= ~ 1 | DiseaseStatus),
                         data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

summary(model_noPUVA.null)
noPUVA.ints <- intervals(model_noPUVA.null, which="fixed")$fixed
noPUVA.ints
noPUVA.ageEff=noPUVA.ints["Age_at_sampling", "est."]
noPUVA.low <-noPUVA.ints["Age_at_sampling", "lower"] 
noPUVA.upp <- noPUVA.ints["Age_at_sampling", "upper"] 

write.table(clone_burden[,c("CloneID", "noPUVA", "Age_at_sampling","DiseaseStatus")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig2b.txt", sep=""), quote=F, row.names = F, sep = "\t")

pB <- ggplot(clone_burden, aes(y=noPUVA, x=Age_at_sampling, fill=DiseaseStatus)) + geom_point( colour="black", shape=21, size=0.5) + 
  scale_fill_manual(values=type_colours)  + 
  labs(y="Number of substitutions", x="Age at sample donation") +
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + theme(legend.title = element_blank(), legend.position = "none") +
  geom_ribbon(aes(ymin=Age_at_sampling*noPUVA.low, ymax=Age_at_sampling*noPUVA.upp, x=Age_at_sampling), alpha = 0.3, show.legend=F) +
  geom_line(aes(y=Age_at_sampling*noPUVA.ageEff, x=Age_at_sampling)) + 
  guides(fill = guide_legend(override.aes = list(size=2.5))) + 
  facet_wrap(~DiseaseStatus, nrow = 2)

noPUVA.ints
coef(summary(model_noPUVA.null))

## Plot the burden considering only UV-light mutations
clone_burden$UV <- clone_burden$SBS7b + clone_burden$SBS7c
model_UV.null <- lme(fixed = UV ~ Age_at_sampling  + MetaLocation, 
                     random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                     weights = varIdent(form= ~ 1 | DiseaseStatus),
                     data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

model_UV.dur <- lme(fixed = UV ~ Age_at_sampling  + MetaLocation + Disease_duration, 
                    random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                    weights = varIdent(form= ~ 1 | DiseaseStatus),
                    data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

summary(model_UV.null)
summary(model_UV.dur)
anova(model_UV.null,model_UV.dur, test=T)$"p-value"[2]
UV.ints <- intervals(model_UV.dur, which="fixed")$fixed
UV.ints

UV.ageEff=UV.ints["Age_at_sampling", "est."]
UV.low <-UV.ints["Age_at_sampling", "lower"] 
UV.upp <- UV.ints["Age_at_sampling", "upper"] 

pC <- ggplot(clone_burden, aes(y=UV, x=Age_at_sampling, fill=DiseaseStatus)) + geom_point( colour="black", shape=21, size=0.5) + 
  scale_fill_manual(values=type_colours)  + 
  labs(y="#UV SBS", x="Age at sample donation") +
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + theme(legend.title = element_blank(), legend.position = "top") +
  geom_ribbon(aes(ymin=Age_at_sampling*UV.low, ymax=Age_at_sampling*UV.upp, x=Age_at_sampling), alpha = 0.3, show.legend=F) +
  geom_line(aes(y=Age_at_sampling*UV.ageEff, x=Age_at_sampling)) + 
  guides(fill = guide_legend(override.aes = list(size=2.5))) + 
  facet_wrap(~DiseaseStatus, nrow=2)

## Look at the mutation burden of the clock-like SBS1 and SBS5

model_clock.null <- lme(fixed = SBS1.5 ~ Age_at_sampling  + MetaLocation, 
                        random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                        weights = varIdent(form= ~ 1 | DiseaseStatus),
                        data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

model_clock.dur <- lme(fixed = SBS1.5 ~ Age_at_sampling  + MetaLocation + Disease_duration, 
                       random = list(PatientID = pdSymm(form = ~ Age_at_sampling - 1), BiopsyID = pdSymm(form = ~ Age_at_sampling - 1)), 
                       weights = varIdent(form= ~ 1 | DiseaseStatus),
                       data = clone_burden[!is.na(clone_burden$Disease_duration),], method="ML")

summary(model_clock.null)
summary(model_clock.dur)
anova(model_clock.null,model_clock.dur, test=T)$"p-value"[2]
clock.ints <- intervals(model_clock.dur, which="fixed")$fixed
clock.ints

clock.ageEff=clock.ints["Age_at_sampling", "est."]
clock.low <-clock.ints["Age_at_sampling", "lower"] 
clock.upp <- clock.ints["Age_at_sampling", "upper"] 

write.table(clone_burden[,c("CloneID", "SBS1.5", "Age_at_sampling","DiseaseStatus")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig2c.txt", sep=""), quote=F, row.names = F, sep = "\t")

pC <- ggplot(clone_burden, aes(y=SBS1.5, x=Age_at_sampling, fill=DiseaseStatus)) + geom_point( colour="black", shape=21, size=0.5) + 
  scale_fill_manual(values=type_colours)  + 
  labs(y="SBS1/5 mutations", x="Age at sample donation") +
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + theme(legend.title = element_blank(), legend.position = "none") +
  geom_ribbon(aes(ymin=Age_at_sampling*clock.low, ymax=Age_at_sampling*clock.upp, x=Age_at_sampling), alpha = 0.3, show.legend=F) +
  geom_line(aes(y=Age_at_sampling*clock.ageEff, x=Age_at_sampling)) + 
  guides(fill = guide_legend(override.aes = list(size=2.5))) + 
  facet_wrap(~DiseaseStatus, nrow=2)

## Collect the fixed effects for plotting: 

effects <- data.frame(Type=c("noPUVA", "UV", "clock","noPUVA", "UV", "clock"), 
                      Effect=c(noPUVA.ageEff, UV.ageEff, clock.ageEff,noPUVA.ints[7,2], UV.ints[7,2], clock.ints[7,2]), 
                      lower=c(noPUVA.low, UV.low, clock.low, noPUVA.ints[7,1], UV.ints[7,1], clock.ints[7,1]),
                      upper=c(noPUVA.upp, UV.upp, clock.upp, noPUVA.ints[7,3], UV.ints[7,3], clock.ints[7,3]),
                      group=c("Age","Age","Age","Disease \nDuration","Disease \nDuration","Disease \nDuration"))

effects$Type[effects$Type=="noPUVA"] <- "Total SBS \n without psoralen"
effects$Type[effects$Type=="clock"] <- "Clock-like"
effects$Type <- factor(effects$Type, levels=c("Clock-like","UV", "Total SBS \n without psoralen"))

write.table(effects, file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig2d.txt", sep=""), quote=F, row.names = F, sep = "\t")

pD <- ggplot(effects, aes(y=Effect, x=Type, colour=group)) + geom_point(position=position_dodge(.25), size=2) + 
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family))+ 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.25), size=1.2) + 
  coord_flip() + labs(x="", y="Substitutions per year", colour="") + 
  geom_hline(yintercept = 0, linetype="dashed") + theme(legend.position = "top") +
  scale_color_brewer(palette = "Paired")




noPUVA.ageEff=noPUVA.ints["Age_at_sampling", "est."]
noPUVA.low <-noPUVA.ints["Age_at_sampling", "lower"] 
noPUVA.upp <- noPUVA.ints["Age_at_sampling", "upper"] 
UV.ints[2,]    ## UV Age effect
UV.ints[7,]    ## UV disease duration effect
clock.ints[2,] ## Clock-like age effect
clock.ints[7,] ## Clock-like disease duration effect

## Mutational signature exposure
sig_colours <- c("#F5BB00","#53ADFC","#FCECC9","#3BB273", "#DD3B46","#C46BAE", "#0566C6","#B5ABAB", 
                 "#F17300", "#A7CECB", "#2D2A32",
                 "#247BA0")


signature_burden <- clone_burden_original
## Merge the unknown components with the Unassigned
## component because we don't really believe they are real
signature_burden$Unassigned <- signature_burden$Unassigned + signature_burden$Unknown.N1 + signature_burden$Unknown.N2 + signature_burden$Unknown.N3
signature_burden$Unknown.N1 <- NULL
signature_burden$Unknown.N2 <- NULL 
signature_burden$Unknown.N3 <- NULL
signature_burden <- signature_burden[order(signature_burden$TotalSBS_adj),]
# Clones with very few mutations were not included in the signature extraction
# Assign all mutations to the Unassigned component.
signature_burden$Unassigned[signature_burden$TotalSBS_adj<20]<- signature_burden$TotalSBS_adj[signature_burden$TotalSBS_adj<20]

#### Looking at the psoralen signature

tmp <- signature_burden[signature_burden$PatientID %in% patient_meta$Patient.ID[patient_meta$Shows_Psoralen_signature],c(colnames(signature_burden)[c(1,3:11)])]
tmp$nr <- c(1:nrow(tmp))

signature_burden_m <- melt(tmp, id.vars=c("CloneID", "nr", "TotalSBS_adj", "PatientID"))

## Order patients by highest mutation burden:
signature_burden_m$PatientID <- factor(signature_burden_m$PatientID, levels=names(tapply(signature_burden_m$TotalSBS_adj, signature_burden_m$PatientID, max))
                                       [order(tapply(signature_burden_m$TotalSBS_adj, signature_burden_m$PatientID, max), decreasing = T)])

sigNames <- as.character(unique(signature_burden_m$variable))
sigNames[4] <- "SBS1/5"
levels(signature_burden_m$variable)[4] <- "SBS1/5"

write.table(clone_burden, file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig2a.txt", sep=""), quote=F, row.names = F, sep = "\t")

pA1 <- ggplot(signature_burden_m, aes(x=reorder(CloneID, -nr), y =TotalSBS_adj/(length(sigNames))))  + 
  geom_bar(stat="identity", width=1) + theme_bw()+  theme(text=element_text(size=text_size, family = font_family)) +theme(axis.text.x = element_blank()) + 
  facet_grid(~PatientID, scale="free",space="free_x") + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank()) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),plot.margin = unit(c(0.5,0.5,-0.25,0.5), "cm")) + 
  labs(x="", y="# SBS") 



pA2 <- ggplot(signature_burden_m, aes(fill=variable, y=value, x=reorder(CloneID, -nr))) + 
  theme_bw()+  theme(text=element_text(size=text_size, family = font_family)) +
  geom_bar(position="fill", stat="identity", width = 1) + theme(axis.text.x = element_blank()) + 
  labs(x="Cell clones grouped by patient", y="Proportion", fill="Mutation Signature") +
  facet_grid(~PatientID, scale="free",space="free_x") +
  scale_fill_manual(values=sig_colours[-5]) + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position = "bottom", plot.margin = unit(c(0,0.5,0,0.5), "cm"))



plot_grid(plot_grid(pA1, pA2,rel_heights = c(1,2), nrow=2,align = "v", labels = c("a", "")),plot_grid(pB,pC,pD, nrow=1, labels=c("b","c", "d")), nrow=2)


### 

tmp <- signature_burden[signature_burden$MetaLocation=="Leg",c(colnames(signature_burden)[c(1,3:11)])]
tmp$nr <- c(1:nrow(tmp))

signature_burden_m <- melt(tmp, id.vars=c("CloneID", "nr", "TotalSBS_adj", "PatientID"))

## Order patients by highest mutation burden:
signature_burden_m$PatientID <- factor(signature_burden_m$PatientID, levels=names(tapply(signature_burden_m$TotalSBS_adj, signature_burden_m$PatientID, max))
                              [order(tapply(signature_burden_m$TotalSBS_adj, signature_burden_m$PatientID, max), decreasing = T)])

sigNames <- as.character(unique(signature_burden_m$variable))
sigNames[4] <- "SBS1/5"
levels(signature_burden_m$variable)[4] <- "SBS1/5"

pA1 <- ggplot(signature_burden_m, aes(x=reorder(CloneID, -nr), y =TotalSBS_adj/(length(sigNames)))) + theme_bw() + 
  geom_bar(stat="identity", width=1) + theme_bw()+  theme(text=element_text(size=text_size, family = font_family)) +theme(axis.text.x = element_blank()) + 
  facet_grid(~PatientID, scale="free",space="free_x") + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank()) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),plot.margin = unit(c(0.5,0.5,-0.25,0.5), "cm")) + 
  labs(x="", y="# SBS") 



pA2 <- ggplot(signature_burden_m, aes(fill=variable, y=value, x=reorder(CloneID, -nr))) + theme_bw()+  theme(text=element_text(size=text_size, family = font_family)) + 
  geom_bar(position="fill", stat="identity", width = 1) + theme(axis.text.x = element_blank()) + 
  labs(x="Cell clones grouped by patient", y="Proportion", fill="Mutation Signature") +
  facet_grid(~PatientID, scale="free",space="free_x") +
  scale_fill_manual(values=sig_colours[-5]) + 
  scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"), plot.margin = unit(c(0,0.5,0,0.5), "cm"))



plot_grid(plot_grid(pA1, pA2,rel_heights = c(1,2), nrow=2,align = "v", labels = c("a", "")),plot_grid(pB,pC,pD, nrow=1, labels=c("b","c", "d")), nrow=2)


pdf(paste(working_dir, "Olafsson_Fig2.pdf", sep=""), height=pdf_height_in-2.5, width = pdf_width_in, fonts = font_family)
Main_fig2 <- plot_grid(plot_grid(pA1, pA2,rel_heights = c(1,2), nrow=2,align = "v", labels = c("a", "")),plot_grid(pB,pC,pD, nrow=1, labels=c("b","c", "d")), nrow=2)
print(Main_fig2)
dev.off()


######################################################################
##
##  SUPPLEMENTARY FIGURE 3 - Mutational signatures
######################################################################
write.table(clone_burden, file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig3.txt", sep=""), quote=F, row.names = F, sep = "\t")

for(loc in unique(signature_burden$MetaLocation)) {
  tmp <- signature_burden[signature_burden$MetaLocation==loc,c(colnames(signature_burden)[c(1,3:11)])]
  tmp$nr <- c(1:nrow(tmp))
  
  signature_burden_m <- melt(tmp, id.vars=c("CloneID", "nr", "TotalSBS_adj", "PatientID"))
  levels(signature_burden_m$variable)[4] <- "SBS1/5"
  ## Order patients by highest mutation burden:
  signature_burden_m$PatientID <- factor(signature_burden_m$PatientID, levels=names(tapply(signature_burden_m$TotalSBS_adj, signature_burden_m$PatientID, max))
                                         [order(tapply(signature_burden_m$TotalSBS_adj, signature_burden_m$PatientID, max), decreasing = T)])
  
  sigNames <- as.character(unique(signature_burden_m$variable))
  
  pE1 <- ggplot(signature_burden_m, aes(x=reorder(CloneID, -nr), y =TotalSBS_adj/(length(sigNames)))) + theme_bw() + 
    geom_bar(stat="identity", width=1) + theme_bw()+  theme(text=element_text(size=text_size, family = font_family)) +theme(axis.text.x = element_blank()) + 
    facet_grid(~PatientID, scale="free",space="free_x") + 
    scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
    theme(strip.background = element_blank(),strip.text.x = element_blank()) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.ticks.x = element_blank(),
          axis.line = element_line(colour = "black"),plot.margin = unit(c(0.5,0.5,-0.25,0.5), "cm")) + 
    labs(x="", y="# SBS") + ggtitle(label=paste("Microbiopsies from the", loc)) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ## Set legend.position to "none" for all locations except Abdomen. Set to right for Abdomen.
  pE2 <- ggplot(signature_burden_m, aes(fill=variable, y=value, x=reorder(CloneID, -nr))) + theme_bw()+  theme(text=element_text(size=text_size, family = font_family)) + 
    geom_bar(position="fill", stat="identity", width = 1) + theme(axis.text.x = element_blank()) + 
    labs(x="Cell clones grouped by patient", y="Proportion", fill="Mutation Signature") +
    facet_grid(~PatientID, scale="free",space="free_x") +
    scale_fill_manual(values=sig_colours[-5]) + 
    scale_y_continuous(expand = c(0,0)) + theme(axis.ticks.x = element_blank()) + 
    theme(strip.background = element_blank(),strip.text.x = element_blank()) +
    theme(legend.position = "right", plot.margin = unit(c(0,0.5,0,0.5), "cm"))+
    guides(fill=guide_legend(nrow=4,byrow=TRUE)) + 
    theme(legend.key.size = unit(0.3, "cm"))
}


pdf(paste(working_dir, "Olafsson_ED_Fig3.pdf", sep=""), height=pdf_height_in, width = pdf_width_in, fonts = font_family)
ED_Fig3 <- plot_grid(plot_grid(pA1, pA2,rel_heights = c(4,4), nrow=2,align = "v", labels = c("a", "")),
plot_grid(pB1, pB2,rel_heights = c(4,4), nrow=2,align = "v", labels = c("b", "")),
plot_grid(pC1, pC2,rel_heights = c(4,4), nrow=2,align = "v", labels = c("c", "")),
plot_grid(pD1, pD2,rel_heights = c(4,4), nrow=2,align = "v", labels = c("d", "")),
plot_grid(pE1, pE2,rel_heights = c(4,4), nrow=2,align = "v",axis = "rl", labels = c("e", "")), 
nrow=5)
print(ED_Fig3)
dev.off()

########################################################
##
##  SUPPLEMENTARY FIGURE 2 - HDP COMPONENTS AND SIGNATURE INFO - RELATED TO FIG 2
#######################################################
load(paste(working_dir, "wes_hdp_components.RData", sep=""))

mut_colours <- c("#16bdebff", "#000000ff", "#e22926ff", "#999999ff", "#9fce67ff", "#ecc6c5ff")
# labels along botom x-axis
trinuc_context <- sapply(strsplit(colnames(mut_count), '\\.'), `[`, 4)
# group labels along the top (and controls colour grouping)
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
                              each=16))

comp_distn <- comp_categ_distn(luad_multi)
ncomp <- nrow(comp_distn$mean) - 1
comp_to_plot <- rownames(comp_distn$mean)
all_clusters <- data.frame(trinuc_context, group_factor)
for(i in seq_along(comp_to_plot)) {
  cname <- comp_to_plot[i]
  sig <- comp_distn$mean[cname, ]
  all_clusters <- cbind(all_clusters, sig)
}
colnames(all_clusters) <- c("Context","Mutation","Unassigned","SBS7b - UV","Psoralens", "SBS1/5", "Unknown N1",
                            "Unknown N2","SBS2 - APOBEC","Unknown N3", "SBS7c",
                            "SBS13 - APOBEC")

write.table(all_clusters, file=paste(working_dir, "source_data/Olafsson_SourceDate_ED_Fig2.txt", sep=""), quote=F, row.names = F, sep = "\t")

pdf(paste(working_dir,"Olafsson_ED_Fig2.pdf", sep=""), height = pdf_height_in+1.5, width=pdf_width_in, fonts = font_family)
par(mfrow=c(5,2), mar=c(3, 3, 1, 1), cex=0.6)
plot_comp_distn(luad_multi, cat_names=trinuc_context,
                grouping=group_factor, col=mut_colours,
                show_group_labels=TRUE, 
                plot_title = c("Unassigned","SBS7b - UV","Psoralens", "SBS1/5", "Unknown N1",
                               "Unknown N2","SBS2 - APOBEC","Unknown N3", "SBS7c",
                               "SBS13 - APOBEC"), cred_int=F)

dev.off()

##########################################################
##
##  FIGURE 3 - CHARACTERISATION OF THE PUVA SIGNATURE
##########################################################

load(paste(working_dir,"puva_characterization.data.RData", sep=""))
library(MutationalPatterns)

# The 192 class profile of the PUVA signature.
#puva_hdp_component <- read.table("/lustre/scratch119/humgen/projects/psoriasis/signature_extraction/puva_characterization/PUVA_signature_hdp_component.txt", h=T)

hartwig_sample <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/HMF003357A_192_trinucs.txt", h=T)
psoralen_196_counts <- merge(mut_mat_s[,16], hartwig_sample, by=0)
rownames(psoralen_196_counts) <- psoralen_196_counts$Row.names
psoralen_196_counts$Row.names <- NULL
colnames(psoralen_196_counts)[1] <- "P34H_LL_4"

psoralen_196_counts$strand <- "Transcribed"
psoralen_196_counts$strand[grep("untranscribed", rownames(psoralen_196_counts))] <- "Untranscribed"
psoralen_196_counts$context <- rownames(psoralen_196_counts)
psoralen_196_counts$context <- gsub("-untranscribed", "", psoralen_196_counts$context)
psoralen_196_counts$context <- gsub("-transcribed", "", psoralen_196_counts$context)
psoralen_196_counts$substitution <- substr(psoralen_196_counts$context, 3,5)
psoralen_196_counts$context <- gsub("\\[", "", psoralen_196_counts$context)
psoralen_196_counts$context <- gsub("\\]", "", psoralen_196_counts$context)
psoralen_196_counts$context <- gsub(">", "", psoralen_196_counts$context)
psoralen_196_counts$context <-paste(substr(psoralen_196_counts$context, 1,2), substr(psoralen_196_counts$context, 4,4), sep="")

## Convert counts to probability
psoralen_196_counts$P34H_LL_4 <- psoralen_196_counts$P34H_LL_4/sum(psoralen_196_counts$P34H_LL_4)
psoralen_196_counts$HMF003357A <- psoralen_196_counts$HMF003357A/sum(psoralen_196_counts$HMF003357A)

ymax <- max(c(psoralen_196_counts$P34H_LL_4, psoralen_196_counts$HMF003357A))
width <- 1
spacing <- 0
colors <- c("#16bdebff", "#000000ff", "#e22926ff", "#999999ff", "#9fce67ff", "#ecc6c5ff")
psoralen_196_counts

psoralens_m <- melt(psoralen_196_counts)

write.table(psoralen_196_counts, file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig3a.txt", sep=""), quote=F, sep = "\t")

pA <- ggplot(psoralens_m, aes(x=context, y=value, fill=substitution, alpha=strand)) +
  geom_bar(stat = "identity", colour = "black", size = 0.1, position = "dodge",width = width) + ylab("Relative contribution") + 
  scale_y_continuous(breaks = seq(0, ymax, 0.03),expand = c(0, 0), limits = c(0,0.13)) +
  scale_alpha_discrete(range = c(0.25, 1)) + scale_fill_manual(values = colors) + 
  facet_grid(variable~substitution, scales="free_x")  + theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  +
  theme(axis.title.y = element_text(size = text_size, vjust = 1), axis.text.y = element_text(size = text_size), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text.x = element_text(size = text_size), strip.text.y = element_text(size = text_size), 
        panel.grid.major.x = element_blank(), panel.spacing.x = unit(spacing, "lines")) + 
  theme(axis.ticks.x = element_blank()) + theme(strip.background = element_rect(fill="white")) + 
  theme(legend.position = "none")



# Extended sequence context plot
## Now a supplementary figure
#pB <- plot_profile_heatmap(mut_mat_ext_context)

# Transcription coupled damage
strand_colours <- c("#003049", "#D62828")

tcd <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/puva_tcd.txt", h=T, sep = "\t")
tcd$DistanceFromTSS <- gsub("kb", "", tcd$DistanceFromTSS)
tcd$DistanceFromTSS <- as.numeric(tcd$DistanceFromTSS)

pB <- ggplot(tcd, aes(x=DistanceFromTSS, y=standardized, colour=expression)) + geom_point(size=0.7) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  +
  geom_line(aes(x=DistanceFromTSS,y=standardized, color = expression, group = class, linetype=MutationType)) + 
  labs(x="Distance from the TSS [kb]", y="Mutation rate relative \n to intergenic regions", colour="", linetype="") + 
  scale_y_continuous(trans="log2") + scale_colour_manual(values=c("#F79256","#FBD1A2", "#7DCFB6","#00B2CA","#1D4E89")) + 
  theme(axis.text.x = element_text(size=text_size)) + 
  scale_x_continuous(breaks=pretty(tcd$DistanceFromTSS, 10)) + 
  theme(legend.spacing.y = unit(-0.35, "cm"))

tcd$class <- gsub("\n","", tcd$class)
tcd$MutationType <- gsub("\n","", tcd$MutationType)
write.table(tcd, file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig3c.txt", sep=""), quote=F, row.names=F,sep = "\t")


## Mutation rate as a function of gene expression
expression_res <- read.table("/lustre/scratch126/humgen/projects/psoriasis/signature_extraction/puva_characterization/Mutation_rate_pr_expr_bin.txt", h=T)
expression_res$scaled <- expression_res$Mutation_Rate/expression_res$Mutation_Rate[1]

write.table(expression_res, file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig3d.txt", sep=""), quote=F, row.names=F,sep = "\t")

pC <- ggplot(expression_res, aes(x=Bin, y=scaled)) + geom_bar(position="dodge",stat="identity", fill="#E24E1B", colour="black") + 
  scale_y_continuous(expand=c(0,0), limits = c(0,1.2)) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  + 
  labs(y="Scaled mutation rate \n at TpA sites", x="Genes binned by \nincreasing expression") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  geom_hline(yintercept = 1, linetype=2)


## Mutation burden as a function of the number of PUVA cycles
clone_burden <- read.table(paste(working_dir, "Supplementary_material/Supplementary_Table3_clone_mutationBurden.txt", sep=""), h=T)
clone_burden <- merge(clone_burden, microd_meta[,c("SampleID", "BiopsyID", "MetaLocation",
                                                   "PatientID", "DiseaseStatus")], by.x="HighCellFrac_sample", by.y="SampleID")
clone_burden <- merge(clone_burden, patient_meta[,c("Patient.ID", "Age_at_sampling", "Disease_duration",
                                                    "Sex", "BMI", "Smoking", "PASI")], by.x="PatientID", by.y="Patient.ID")

clone_burden$Disease_duration[clone_burden$DiseaseStatus=="Non-lesional" & !is.na(clone_burden$Disease_duration)] <- 0

clone_burden <- merge(clone_burden, patient_meta[,c(1,22:27)], by.y="Patient.ID", by.x="PatientID")
clone_burden$AmountUVB <- factor(clone_burden$AmountUVB, levels=c("None","unknown", "<=50","51-200"))
clone_burden$AmountPuva <- factor(clone_burden$AmountPuva, levels=c("None","unknown", "<=50","51-200",">200"))

write.table(clone_burden[,c("PatientID","HighCellFrac_sample","CloneID","AmountPuva","PUVA")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig3b.txt", sep=""), quote=F, row.names=F,sep = "\t")

pD <- ggplot(clone_burden, aes(x=AmountPuva, y=PUVA+1, fill=AmountPuva))  +
  geom_boxplot(outlier.size = 0.2) +  theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  + scale_fill_brewer(palette="Set1") + theme(legend.position = "none") + 
  labs(x="Number of PUVA cycles", y="Psoralen \n mutation burden") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 

puva_colours <- c("#FCBF49", "#40476D")

microd_meta <- merge(microd_meta, patient_meta[, c("Patient.ID","Shows_Psoralen_signature","EverPuva","AmountPuva","EverUVB","AmountUVB","EverMethotrexate","EverSteroids")], by.x="PatientID", by.y="Patient.ID", all.x=T)

microd_meta$AmountUVB <- factor(microd_meta$AmountUVB, levels=c("None","unknown", "<=50","51-200"))
microd_meta$AmountPuva <- factor(microd_meta$AmountPuva, levels=c("None","unknown", "<=50","51-200",">200"))

mVAFs <- data.frame(tapply(microd_meta$MedianVAF, microd_meta$PatientID, median, na.rm=T))
mVAFs$PatientID <- rownames(mVAFs)
colnames(mVAFs)[1] <- "MedianVAF"
patient_meta <- merge(patient_meta, mVAFs, by.x="Patient.ID", by.y="PatientID")
patient_meta$labels <- ifelse(patient_meta$Shows_Psoralen_signature, "Psoralen signature seen", "Signature not seen")

write.table(patient_meta[,c("Patient.ID","Age_at_sampling","MedianVAF","labels")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig3f.txt", sep=""), quote=F, row.names=F,sep = "\t")

pE <- ggplot(patient_meta, aes(x=Age_at_sampling, y=MedianVAF, colour=labels)) + geom_point(size=0.5) +
  geom_smooth(method="lm") + theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  + 
  labs(y="Median VAF \n of microbiopsies", x="Patient Age", colour="") +
  theme(legend.position = "top") + scale_colour_manual(values=rev(puva_colours)) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) + theme(legend.key.size = unit(0.25, "cm"), legend.box.spacing = unit(0, "pt"))


## Signature potential damage: 
sig_damage <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/signature_potential_damage.txt", h=T, sep="\t")
sig_damage$Reference <- factor(sig_damage$Reference, levels=c("Uniform mutation rate", "SBS7b"))
## Nonsense,Splice, Missense, Indels, synonymous
#annotation_cols <- c("#5D378E","#7155A5","#6A9A9E","#C0843E","#A9ACAB")
# Missense, Nonsense, Splice, synonymous
annotation_cols <- c("#6A9A9E","#5D378E","#7155A5","#A9ACAB")

write.table(sig_damage[,c("type", "sig","ratio_by_background", "Reference")], file=paste(working_dir, "source_data/Olafsson_SourceDate_Fig3e.txt", sep=""), quote=F, row.names=F,sep = "\t")

pF <- ggplot(sig_damage, aes(x=sig, y=ratio_by_background, fill=type)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  + 
  geom_hline(yintercept = 1, linetype="dashed") + 
  scale_y_continuous(expand=c(0,0), limits = c(0,1.8)) + labs(x="", y="Ratio of mutation type \n caused by signature \n vs background ", fill="") + facet_wrap(~Reference) + 
  scale_fill_manual(values=annotation_cols) + theme(legend.position = "bottom") + 
  theme(legend.key.size = unit(0.25, "cm"), legend.box.spacing = unit(-10, "pt")) 

pdf(paste(working_dir, "Olafsson_Fig3.pdf", sep=""), height=pdf_height_in, width = pdf_width_in, fonts = font_family)

main_fig3 <- plot_grid(pA,plot_grid(pD,pE, rel_widths = c(2,1), labels=c("b","f")), pB, plot_grid(pC, pF, nrow=1, labels=c("d","e"), rel_widths = c(1,2.5)), labels=c("a", "", "c", ""), ncol = 1, rel_heights = c(1,1,1.35,1))
print(main_fig3)
dev.off()



# Effects of PUVA on clonal expansions
puva_colours <- c("#FCBF49", "#40476D")
clone_sigs <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/clone_mutation_burden.txt", h=T)
x <- clone_sigs[clone_sigs$Psoralens>100,]
puva_exp_patients <- unique(unlist(strsplit(x$CloneID, spli="_"))[c(T,F)])
microd_meta$PUVA_exposed <- F
microd_meta$PUVA_exposed[microd_meta$PatientID %in% puva_exp_patients] <- T

puva_vaf_test <- wilcox.test(microd_meta$MedianVAF~microd_meta$PUVA_exposed)

pF <- ggplot(microd_meta, aes(x=PUVA_exposed, y=MedianVAF)) + geom_boxplot(aes(fill=PUVA_exposed)) +
  theme_classic(base_size = BASESIZE) + labs(x="", y="Median VAF of microbiopsy") +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=c("Not exposed", "Psoralen Exposed")) + 
  scale_fill_manual(values=puva_colours) + scale_y_continuous(limits = c(0,0.6)) +
  geom_signif(annotation = paste("P = ", formatC(puva_vaf_test$p.value, digits=2), sep=""),
              y_position = c(0.55), xmin = c(1), xmax = c(2),
              tip_length = c(0.03, 0.03))


chisq.test(data_frame$treatment, data_frame$improvement, correct=FALSE)

shared_by_distance <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/shared_by_distance.txt", h=T)
shared_by_distance$inSquare <- "No"
shared_by_distance$inSquare[shared_by_distance$distances>250 & shared_by_distance$frac_shared>0.2] <- "Yes"

chisq.test(shared_by_distance$PUVA_exposed, shared_by_distance$inSquare, correct=T)


sq_test <- chisq.test(shared_by_distance$PUVA_exposed, shared_by_distance$inSquare, correct=T)
table(shared_by_distance$PUVA_exposed, shared_by_distance$inSquare)

pG <- ggplot(shared_by_distance, aes(x=distances, y=frac_shared)) + geom_point(aes(colour=PUVA_exposed), alpha=0.75) + 
  labs(y="Fraction of shared mutations", x="Micrometers separating microbiopsies", colour="Psoralen exposed") + 
  theme_bw(base_size = BASESIZE) + theme(legend.position = "top") + 
  scale_colour_manual(values=puva_colours) +
  annotate("rect", xmin = 500, xmax = 3000, ymin = 0.1, ymax = 1,
           alpha = .2) 
#annotate(geom="text", x=2200, y=0.8, label=paste("P = ", format(sq_test$p.value, digits = 2), sep=""))

shared_by_distance$type <- "Lesional"
shared_by_distance$type[grep("H", shared_by_distance$pairs)] <- "Non-lesional"

ggplot(shared_by_distance, aes(x=distances, y=frac_shared, colour=type)) + 
  geom_point()


###########################################
###
### SUPPLEMENTARY FIGURE 4 - MORE ON PUVA
############################################
# a) Indel profile of exposed vs not
# b) Extended sequence context
# c) Leading vs lagging strand
# d) Replication timing bins

indel_counts <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/indel_counts_pr_patient.txt", h=T)

exp <- indel_counts[,colnames(indel_counts) %in% patient_meta$Patient.ID[patient_meta$Shows_Psoralen_signature]]
nexp <- indel_counts[,!(colnames(indel_counts) %in% patient_meta$Patient.ID[patient_meta$Shows_Psoralen_signature])]

sums <- data.frame(Exposed=rowSums(exp), NonExposed=rowSums(nexp))

write.table(sums, file=paste(working_dir, "source_data/Olafsson_SourceData_ED_Fig4a.txt", sep=""), quote=F, row.names=T,sep = "\t")

pA <- plot_indel_contexts(sums, condensed = TRUE) + 
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + 
  theme(legend.position = "bottom") + labs(fill="") + theme(legend.key.size = unit(0.25, "cm")) +
  theme(legend.box.spacing = unit(-10, "pt")) + scale_y_continuous(expand=c(0,0))


### Extended sequence context
library(dplyr)
fullcontext <- l_context <- r_context <- muttype <- NULL
nrmuts <- rel_nrmuts <- NULL
mut_matrix <- mut_mat_ext_context

tb <- mut_matrix %>% as.data.frame() %>% tibble::rownames_to_column("fullcontext") %>% 
  tidyr::pivot_longer(-fullcontext, names_to = "sample", 
                      values_to = "nrmuts") %>% tidyr::separate("fullcontext", 
                                                                into = c("l_context", "muttype", "r_context"), sep = "\\[|\\]") %>% 
  dplyr::mutate(mut = factor(muttype, levels = unique(muttype)), 
                r_context = factor(r_context, levels = unique(r_context)), 
                l_context = factor(l_context, levels = Biostrings::reverse(unique(l_context))))
tb <- tb %>% dplyr::group_by(sample) %>% dplyr::mutate(rel_nrmuts = nrmuts/sum(nrmuts)) %>% 
  dplyr::ungroup()

spacing=0.5
axis_size=8

tb$muttype <- factor(tb$muttype, levels=c("C>A","C>G","C>T","T>A","T>C","T>G"))
library(RColorBrewer)

write.table(tb, file=paste(working_dir, "source_data/Olafsson_SourceData_ED_Fig4b.txt", sep=""), quote=F, row.names=F,sep = "\t")


pB <- ggplot(tb, aes(x = r_context, y = l_context, fill = rel_nrmuts)) + 
  geom_raster() + scale_fill_distiller(palette = "YlOrRd", 
  direction = 1, name = "Relative contribution", limits = c(0,0.0362)) + 
  facet_grid( ~ muttype) + 
  labs(x = "3' context", y = "5' context") + theme_classic() + 
  theme(text=element_text(size=text_size, family = font_family)) +
  theme(axis.text.y = element_text(size = 5.5), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
        size = 5.5), panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(), panel.spacing.x = unit(spacing,"lines"), 
        panel.spacing.y = unit(spacing, "lines"),
        legend.position = "bottom") + theme(legend.box.spacing = unit(-10, "pt"))

# Replication strand bias
strand_counts_rep$strand[strand_counts_rep$strand=="left"] <- "Leading"
strand_counts_rep$strand[strand_counts_rep$strand=="right"] <- "Lagging"
strand_counts_rep$strand <- factor(strand_counts_rep$strand,levels=c("Leading", "Lagging"))
strand_bias_rep <- strand_bias_test(strand_counts_rep)

write.table(strand_counts_rep, file=paste(working_dir, "source_data/Olafsson_SourceData_ED_Fig4c.txt", sep=""), quote=F, row.names=F,sep = "\t")

pC <- plot_strand(strand_counts_rep, mode = "absolute") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,44000)) + 
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + theme(strip.background = element_blank()) +
  theme(strip.text = element_blank()) + guides(fill=guide_legend(title=""), alpha=guide_legend(title="Strand")) + 
  geom_signif(annotation = formatC(strand_bias_rep$p_poisson[c(4:6)], digits=2),
              y_position = c(28000,40000,19000), xmin = c(3.7, 4.7,5.7), xmax = c(4.3,5.3,6.3),
              tip_length = c(0.03, 0.03), size=0.25, textsize=2, family = font_family, parse=T) + theme(legend.key.size = unit(0.25, "cm")) + 
  labs(y="Total number \n of mutations") + theme(legend.box.spacing = unit(0, "pt")) +
  theme(legend.spacing.y = unit(0.05, "cm"))

# Replication timing across all samples
colnames(rt_results)[3] <- "MicrobiopsyID"
rt_results$V1 <- gsub("RT_", "Bin-", rt_results$V1)

write.table(rt_results, file=paste(working_dir, "source_data/Olafsson_SourceData_ED_Fig4d.txt", sep=""), quote=F, row.names=F,sep = "\t")


pD <- ggplot(rt_results, aes(x=V1, y=log10(V2))) + geom_point() +
  geom_line(aes(group=MicrobiopsyID, colour=MicrobiopsyID)) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family))  +
  labs(y="log10(Mutations per \n TpA/ApT site)", x="Bins of ascending replication timing") +
  theme(legend.position = "none") + theme(plot.margin = margin(5.5,5.5,5.5,15, "points"))


pdf(paste(working_dir, "Olafsson_ED_Fig4.pdf", sep=""), height=pdf_height_in, width = pdf_width_in, fonts = font_family)
ED_fig4 <- plot_grid(pA,pB,plot_grid(pC,pD, labels = c("c","d")), nrow = 3, labels=c("a","b",""), rel_heights = c(2.5,1.5,1))

print(ED_fig4)
dev.off()


#####################################
### FIGURE 4 - POSITIVE SELECTION
#####################################

## Nonsense,Splice, Missense, Indels, synonymous
annotation_cols <- c("#5D378E","#7155A5","#6A9A9E","#C0843E","#A9ACAB")
sel_cv <- read.table(paste(working_dir, "sel_cv_dNdS_results.txt", sep=""), h=T)
signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c(1:6,22,19)]


dnds_of_signif_genes <- sel_cv[sel_cv$gene_name %in% signif_genes$gene_name, c("gene_name", "wmis_cv", "wnon_cv",
                                                                                   "wind_cv", "wdbs_cv")]

## A) Number of mutations in significant genes
## Note: Use the annotation from https://www.cbioportal.org/mutation_mapper because
## dNdS doesn't properly annotate DBS mutations. 

annotated_muts <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/MutationMapper_annotated_drivers.tsv", sep="\t", h=T)
## dNdS and MutationMapper don't always use the same transcripts and sometimes the annotation
## isn't the same. This explains the Intron -> Missense below. 
annotated_muts$Mutation.Type[annotated_muts$Mutation.Type=="Intron"] <- "Missense"
annotated_muts$Mutation.Type[annotated_muts$Mutation.Type=="Missense_Mutation"] <- "Missense"
annotated_muts$Mutation.Type[annotated_muts$Mutation.Type %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del")] <- "Indel"
annotated_muts$Mutation.Type[annotated_muts$Mutation.Type %in% c("Splice_Region", "Splice_Site")] <- "Splice"
annotated_muts$Mutation.Type[annotated_muts$Mutation.Type=="Nonsense_Mutation"] <- "Nonsense"
annotated_muts$Mutation.Type[annotated_muts$Mutation.Type=="Silent"] <- "Synonymous"
annotated_muts$Mutation.Type <- factor(annotated_muts$Mutation.Type, levels=c("Nonsense","Splice", "Missense", "Indel", "Synonymous"))

class_counts <- data.frame(table(annotated_muts$Gene, annotated_muts$Mutation.Type))
class_counts$Var1 <- factor(class_counts$Var1, levels=c("NOTCH1", "FAT1", "PPM1D", "TP53", "NOTCH2", "CHEK2", "GXYLT1", "ZFP36L2","EEF1A1"))

write.table(class_counts, file=paste(working_dir, "source_data/Olafsson_SourceData_Fig4a.txt", sep=""), quote=F, row.names=F,sep = "\t")

pA <- ggplot(class_counts, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity") +
  scale_y_continuous(expand=c(0,1), limits=c(0,250)) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + labs(y="# Mutations", x="", fill="") + 
  scale_fill_manual(values=annotation_cols) + theme(axis.text.x = element_text(angle=90)) + 
  theme(legend.position = c(0.8,0.75)) + theme(panel.grid = element_blank()) + 
  theme(legend.key.size = unit(0.25, "cm"))

## B) dNdS ratios of different mutation classes
dnds_m <- melt(dnds_of_signif_genes)
dnds_m$variable <- factor(dnds_m$variable, labels = c("Missense", "Nonsense+Splice", "Indels", "DBS"))
dnds_m$gene_name <- factor(dnds_m$gene_name, levels=levels(class_counts$Var1))

write.table(dnds_m, file=paste(working_dir, "source_data/Olafsson_SourceData_Fig4b.txt", sep=""), quote=F, row.names=F,sep = "\t")

pB <- ggplot(dnds_m, aes(x=gene_name, y=value, fill=variable)) + geom_bar(position = "dodge", stat="identity") + 
  scale_y_continuous(expand=c(0,0), limits = c(0,60)) + theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + labs(y="Observed/Expected", x="", fill="") + 
  scale_fill_manual(values=c(annotation_cols[c(3,1,4)],"#D0A3BF")) + theme(axis.text.x = element_text(angle=90)) + 
  theme(legend.position = c(0.8,0.8)) + theme(panel.grid = element_blank()) + 
  theme(legend.key.size = unit(0.25, "cm"))

## C) Pathway-level dNdS
mis_and_trunc <- read.table("/lustre/scratch126/humgen/projects/psoriasis/selection_analyses/pathway_dnds_mis_and_trunc_w_covs.txt", h=T)

mis_and_trunc$name <- factor(mis_and_trunc$name, levels=c("wmis", "wtru"), labels = c("Missense", "Truncating"))
mis_and_trunc$geneList <- factor(mis_and_trunc$geneList, levels=c("Normal_skin_pos", "BCC", "GWAS_psoriasis","IBD_mucosa", "TNF", "IFNg", "IL12_23", "IL36_MyD88",
                                                                  "TLR","IL17", "MHC_classI"),
                                 c("Normal skin/SCCs", "BCCs", "Psoriasis GWAS hits", "IBD mucosa", "TNF", "IFNg", "IL12/23", "IL36/MyD88", "TLR", "IL17", "MHC-class I"))

library(scales)

pC <- ggplot(mis_and_trunc, aes(x=geneList, y=mle, ymin = cilow, ymax = cihigh, colour=name)) +
  geom_point(size=3, position=position_dodge(width=0.5)) + 
  geom_errorbar(size=1, position = "dodge", width=0.5) +theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + 
  geom_hline(yintercept = 1, aes(size=2)) + theme(axis.text.x = element_text(angle=75, hjust=1)) +labs(x="", y="dN/dS") + 
  scale_y_continuous(trans=log2_trans()) + theme(legend.position = "top") + labs(colour="") + 
  scale_colour_manual(values=annotation_cols[c(3,1)])


### Pathway-level dN/dS when implementing the pentanucleotide model 

geneLists <- read.table("/lustre/scratch126/humgen/projects/psoriasis/selection_analyses/pathway_dNdS_geneLists.txt")

df <- data.frame()
for(geneL in geneLists$V1) {
  results <- read.table(paste("/lustre/scratch126/humgen/projects/psoriasis/selection_analyses/pathway_pentamodel/", geneL, "_Full3075_1x2w_model_dNdSvals.txt", sep=""),h=T)
  results$pathway <- geneL
  df <- rbind(df, results)
}

mis_and_trunc <- df[df$omega %in% c("wmis_driv", "wnon_driv"),]
mis_and_trunc$q <- p.adjust(mis_and_trunc$P, method="BH")

mis_and_trunc$omega <- factor(mis_and_trunc$omega, levels=c("wmis_driv", "wnon_driv"), labels = c("Missense", "Nonsense"))
mis_and_trunc$pathway<- factor(mis_and_trunc$pathway, levels=c("Normal_skin_pos", "BCC", "GWAS_psoriasis","IBD_mucosa", "TNF", "IFNg", "IL12_23", "IL36_MyD88",
                                                                  "TLR","IL17", "MHC_classI"),
                                 c("Normal skin/SCCs", "BCCs", "Psoriasis GWAS hits", "IBD mucosa", "TNF", "IFNg", "IL12/23", "IL36/MyD88", "TLR", "IL17", "MHC-class I"))

write.table(mis_and_trunc, file=paste(working_dir, "source_data/Olafsson_SourceData_Fig4c.txt", sep=""), quote=F, row.names=F,sep = "\t")

pC <- ggplot(mis_and_trunc, aes(x=pathway, y=MLEs, ymin = lowbd, ymax = highbd, colour=omega)) +
  geom_point(size=1.5, position=position_dodge(width=0.5)) + 
  geom_errorbar(size=0.5, position = "dodge", width=0.5) +theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + 
  geom_hline(yintercept = 1, aes(size=2)) + theme(axis.text.x = element_text(angle=75, hjust=1)) +labs(x="", y="dN/dS") + 
  scale_y_continuous(trans=log2_trans()) + theme(legend.position = c(0.8,0.8)) + labs(colour="") + 
  scale_colour_manual(values=annotation_cols[c(3,1)])


## D) Plot the fraction of mutated cells 

cell_frac_lesional <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/mutant_cell_fractions_lesional.txt", h=T)
cell_frac_non_lesional <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/mutant_cell_fractions_non_lesional.txt", h=T)

cell_frac_lesional$type <- "Lesional"
cell_frac_non_lesional$type <- "Non-Lesional"
cell_frac_lesional <- cell_frac_lesional[!is.na(cell_frac_non_lesional$NOTCH1),]
cell_frac_non_lesional <- cell_frac_non_lesional[!is.na(cell_frac_non_lesional$NOTCH1),]


test <- melt(rbind(cell_frac_lesional, cell_frac_non_lesional))
test$variable <- factor(test$variable, levels=levels(class_counts$Var1))

write.table(test, file=paste(working_dir, "source_data/Olafsson_SourceData_Fig4d.txt", sep=""), quote=F, row.names=F,sep = "\t")

pD <- ggplot(test, aes(x=variable, y=value, fill=type)) + geom_boxplot(outlier.size=0.25) + ylim(c(0,1.05)) + 
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + scale_fill_manual(values=type_colours) +
  theme(axis.text.x = element_text(angle=90)) + labs(y="Fraction of mutated cells", x="") + 
  theme(legend.title = element_blank(), legend.position = c(0.8,0.8)) 



wilcox.test(cell_frac_lesional$NOTCH1, cell_frac_non_lesional$NOTCH1, paired = TRUE)
wilcox.test(cell_frac_lesional$NOTCH2, cell_frac_non_lesional$NOTCH2, paired = TRUE)
wilcox.test(cell_frac_lesional$PPM1D, cell_frac_non_lesional$PPM1D, paired = TRUE)
wilcox.test(cell_frac_lesional$FAT1, cell_frac_non_lesional$FAT1, paired = TRUE)
wilcox.test(cell_frac_lesional$TP53, cell_frac_non_lesional$TP53, paired = TRUE)
wilcox.test(cell_frac_lesional$ZFP36L2, cell_frac_non_lesional$ZFP36L2, paired = TRUE)
wilcox.test(cell_frac_lesional$EEF1A1, cell_frac_non_lesional$EEF1A1, paired = TRUE)



pdf(paste(working_dir, "Olafsson_Fig4.pdf", sep=""), width = pdf_width_in,height = pdf_height_in-1.5)
Main_fig4 <- plot_grid(plot_grid(pA, pC, nrow=1, labels=c("a", "c"), rel_widths = c(1,2)), plot_grid(pB,pD,nrow = 1,labels=c("b","d"), rel_widths = c(1,2)), nrow=2)
print(Main_fig4)
dev.off()


##########################################################
##
##  SUPPLEMENTARY FIGURE XX - FRACTION OF MUTANT CELLS
##########################################################

## I make two plots and then merge them with the lollipop plots in Inkscape. 
## the lollipop plots are made in https://www.cbioportal.org/mutation_mapper
mutant_cell_fraction <- read.table("/nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/manuscript_data_and_figures/mutant_cell_fractions.txt", h=T)
mutant_cell_m <- melt(mutant_cell_fraction)
mutant_cell_m$variable <- factor(mutant_cell_m$variable, levels=levels(class_counts$Var1))

write.table(mutant_cell_fraction, file=paste(working_dir, "source_data/Olafsson_SourceData_ED_Fig6a.txt", sep=""), quote=F, row.names=T,sep = "\t")

pdf(paste(working_dir, "Olafsson_ED_Fig6a.pdf", sep=""), width = pdf_width_in/2,height = 1.8, fonts = font_family)
ED_fig6 <- ggplot(mutant_cell_m, aes(x=value)) + geom_histogram(aes(fill=variable)) +
  theme_classic()+  theme(text=element_text(size=text_size, family = font_family)) + scale_y_continuous(expand=c(0,0)) + 
  facet_wrap(~variable, nrow = 2) + scale_fill_brewer(palette = "Set1") + 
  labs(x="Fraction of mutant cells", y="# Patients", fill="") + 
  theme(legend.position = "bottom") + theme(legend.position = "none")

print(ED_fig6)
dev.off()


mutant_cell_fraction$PatientID <- rownames(mutant_cell_fraction)

p <- merge(mutant_cell_fraction, patient_meta, by.x="PatientID", by.y="Patient.ID")

p_m <- melt(p[,c("NOTCH1", "FAT1", "PPM1D", "TP53","NOTCH2","CHEK2","GXYLT1", "ZFP36L2",
                 "Age_at_sampling", "Sex")], id=c("Age_at_sampling", "Sex"))

ggplot(p_m, aes(x=Age_at_sampling, y=value)) + geom_point(aes(colour=variable, shape=Sex)) +
  facet_wrap(.~variable, nrow=2) + geom_smooth(method="lm") + theme_bw(base_size = BASESIZE)+
  labs(x="Age at sampling", y="Fraction of mutant cells", colour="", shape="") + 
  theme(legend.position = "bottom")




