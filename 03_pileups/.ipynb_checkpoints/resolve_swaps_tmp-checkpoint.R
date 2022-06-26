
# bsub -o comb_test.out -e comb_test.err -R"select[mem>4000] rusage[mem=4000]" -M4000 "/software/R-3.6.1/bin/Rscript /nfs/users/nfs_s/so11/phd/psoriasis/bsub_jupyter_lab/psoriasis/03_pileups/resolve_swaps_tmp.R"
.libPaths("/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm5_R3.6_install/")
library(deepSNV)
library(vcfR)
library(VariantAnnotation)

#parameters <- commandArgs(TRUE)

#patient  <- parameters[1]
#sampleVector  <- parameters[2]
#project_vector  <- parameters[3]
#output_dir  <- parameters[4]
#pso_dataDir <- parameters[5]

# For testing:
#sampleVector = c("P01L_LL_11","P03H_T_5","P05L_T_7","P05L_T_6","P07H_UL_1","P07L_UL_1")
sampleVector=c("P01H_LL_6","P01L_LL_9","P01L_LL_10","P01L_LL_11","P01L_LL_12","P01L_LL_13","P01L_LL_14","P03H_T_1","P03H_T_3","P03H_T_5","P03L_T_1","P03L_T_2","P03L_T_3","P03L_T_4","P03L_T_5","P04H_T_6","P04H_T_7","P04L_T_5","P04L_T_6","P04L_T_7","P04L_T_8","P04L_T_9","P04L_T_10","P05H_T_1","P05H_T_3","P05H_T_6","P05L_T_6","P05L_T_7","P05L_T_8","P05L_T_9","P05L_T_10","P05L_T_11","P05L_T_12","P06H_T_1","P06H_T_3","P06L_T_1","P06L_T_2","P06L_T_3","P06L_T_4","P06L_T_5","P06L_T_6","P07H_UL_1","P07H_UL_3","P07H_UL_4","P09H_UA_1","P09H_UA_3","P07L_UL_1","P07L_UL_2","P07L_UL_3","P07L_UL_4","P07L_UL_5","P09L_LA_1","P09L_LA_2","P09L_LA_3","P09L_LA_4","P09L_LA_5","P09L_LA_6","P09L_LA_7","P11H_LA_1","P11H_LA_3","P11H_LA_5","P11L_LA_1","P11L_LA_2","P11L_LA_3","P11L_LA_4","P11L_LA_5","P11L_LA_6","P11L_LA_7","P11L_LA_8","P11L_LA_9","P11L_LA_10","P11L_LA_11","P12H_T_2","P12H_T_3","P12H_T_4","P13H_T_1","P13H_T_3","P13L_T_2","P13L_T_3","P13L_T_4","P13L_T_5","P13L_T_6","P13L_T_7","P17L_LL_2","P17L_LL_3","P17L_LL_4","P17H_LL_1","P23L_T_1","P23L_T_2","P23L_T_3","P23L_T_4","P23L_T_5","P23L_T_6","P23L_T_7","P23L_T_9","P23L_T_10","P23L_T_11","P23H_T_2","P23H_T_3","P23H_T_4","P29L_LL_1","P29L_LL_2","P29L_LL_3","P29L_LL_4","P29L_LL_5","P29L_LL_6","P29L_LL_7","P29L_LL_8","P29L_LL_9","P31L_UA_1","P31L_UA_2","P31L_UA_3","P31L_UA_4","P31L_UA_5","P31L_UA_6","P31L_UA_7","P31L_UA_8","P32L_T_1","P32L_T_2","P32L_T_3","P32L_T_4","P32L_T_6","P32L_T_7","P32H_T_1","P32H_T_2","P32H_T_5","P31H_UA_3","P31H_UA_5","P34L_LL_1","P34L_LL_2","P34L_LL_5","P34L_LL_6","P34L_LL_7","P34L_LL_8","P34L_LL_9","P34L_LL_10","P34L_LL_11","P34H_LL_2","P34H_LL_5","P34H_LL_6","P29H_LL_1","P29H_LL_3","P29H_LL_6","P03L_T_6","P03L_T_7","P03L_T_8","P03L_T_10")


patient <- "all_samples"
project_vector <- rep("2545", length(sampleVector))
output_dir="/lustre/scratch119/humgen/projects/psoriasis/pileups/"
pso_dataDir="/nfs/cancer_ref01/nst_links/live/"



vcf2tab = function(vcf_file, sample_name,PassFilter=TRUE,ASMD_value=140,CLMP_value=0){
   
   data = readVcf(vcf_file)
   reads = cbind((geno(data)$FAZ)[,2]+(geno(data)$RAZ)[,2],(geno(data)$FCZ)[,2]+(geno(data)$RCZ)[,2],
                 (geno(data)$FGZ)[,2]+(geno(data)$RGZ)[,2],(geno(data)$FTZ)[,2]+(geno(data)$RTZ)[,2])
   Var = data.frame(Chr=as.character(seqnames(rowRanges(data))),
                    Pos=start(ranges(rowRanges(data))),
                    Ref=as.character(ref(data)))
   Alt_tmp = CharacterList(alt(data))
   Delete = which(sum(alt(data)%in%c("A","C","T","G"))!=1) # check all variants have an alt allele
   Alt_tmp[Delete]="N"
   Var$Alt = as.character(unlist(Alt_tmp))
   Var$NR=rowSums(reads)
   Var$NV=NA
   colnames(reads)=c("A","C","G","T")
   for (k in c("A","C","G","T")){
     Var$NV[Var$Alt==k] = reads[Var$Alt==k,k]
   }
   
   select <- logical(length=nrow(info(data)))
   select <- !select
   if(PassFilter) select[rowRanges(data)$FILTER!="PASS"] <- FALSE
 
   if(!is.na(ASMD_value)) select[info(data)$ASMD<ASMD_value ] <- FALSE
   
   if(!is.na(CLMP_value)) select[info(data)$CLPM>CLMP_value ] <- FALSE
 
   Var <- Var[select,]
   Var$ID = paste(Var$Chr,Var$Pos,Var$Ref,Var$Alt,sep = "_")
   Var$VAF=Var$NV/Var$NR
   Var$SampleID <- sample_name
   return(Var)
 }


pindel_vcf_to_tab <- function(vcf_file, sample_name) {
  indel_data <- readVcf(vcf_file)
  # See this for explanation of the info fields:
  # https://www.ncbi.nlm.nih.gov/pubmed/26678383
  Var = data.frame(Chr=as.character(seqnames(rowRanges(indel_data))),
                   Pos=start(ranges(rowRanges(indel_data))),
                   Ref=as.character(ref(indel_data)),
                   Alt=as.character(unlist(alt(indel_data))))
  Var$ID = paste(Var$Chr,Var$Pos,as.character(Var$Ref),as.character(Var$Alt),sep = "_")
  Var$VAF <- geno(indel_data)$VAF[, 2]
  Var$Sample <- sample_name
  Var$type <- as.character(unlist(info(indel_data)$VT))

  gene_info <- as.data.frame(do.call(rbind, strsplit(as.character(info(indel_data)$VD),split=c('|'),fixed=TRUE)))
  Var$gene <- as.character(gene_info[,1])
  Var$protein_syntax <- as.character(gene_info$V5)
  Var$CDS_syntax <- as.character(gene_info$V4)
  #Var$annotation <- as.character(info(vcf_file)$VC)
  Var <- Var[rowRanges(indel_data)$FILTER=="PASS",]
  return(Var)
}


## Start doing stuff:
samples_of_interest <- unlist(strsplit(sampleVector, split=","))
project_numbers <- unlist(strsplit(project_vector, split=","))


all_muts = data.frame()
i <- 1
for(sample in samples_of_interest) {
    
    path=paste(pso_dataDir, project_numbers[i], "/", sep="")
    
    caveman_file <- paste(path,sample,"/",sample,".caveman_c.annot.vcf.gz",sep="")
#    pindel_file <- paste(path,sample,"/",sample,".pindel.annot.vcf.gz",sep="")
    vcf <- vcf2tab(caveman_file, sample)
	snvs <- vcf[,c("Chr", "Pos", "Ref", "Alt", "SampleID")]
#    vcf <- pindel_vcf_to_tab(pindel_file, sample)
#    indels <- vcf[,c("Chr", "Pos", "Ref", "Alt", "Sample")]
#    colnames(indels) <- c("Chr", "Pos", "Ref", "Alt", "SampleID")
	if(nrow(all_muts)==0) {
		all_muts = snvs
#        all_muts = rbind(all_muts,indels)
	} else {
		colnames(snvs) = colnames(all_muts)
#        colnames(indels) = colnames(snvs)
		all_muts = rbind(all_muts,snvs)
#        all_muts = rbind(all_muts,indels)
	}
    i <- i+1
}

cat("Mutations:\n")
print(table(all_muts$sampleID))

unique_muts = all_muts[!duplicated(all_muts[,c("Chr","Pos","Ref","Alt")]),]
cat(nrow(unique_muts), "unique_muts\n")

## Here's the twist: Randomly choose just 20k mutations to analyze:
unique_muts <- unique_muts[sample(1:nrow(unique_muts), size=10000),]




unique_muts_ids = paste(unique_muts$Chr,unique_muts$Pos,unique_muts$Ref,unique_muts$Alt,sep=":")
alts = matrix(nrow=nrow(unique_muts),ncol=length(samples_of_interest))
rownames(alts) = unique_muts_ids
colnames(alts) = samples_of_interest
covs = alts

## There is a problem with indels here. See if you can work out how to 
## fix the total alt
for(mut in c(1:nrow(unique_muts))) {
    if (mut%%500==0){
        cat("Mut",mut,"out of", nrow(unique_muts), "...\n")
    }

	ref = unique_muts[mut,"Ref"]
	alt = unique_muts[mut,"Alt"]
	alt_lower = tolower(alt)
	for(sample in samples_of_interest) {	
		bam = paste(path,sample,"/",sample,".sample.dupmarked.bam",sep="")
        # Need to change the second Pos to End, depending on the size of
        # the indel. 
		kk = bam2R(bam,unique_muts[mut,"Chr"],unique_muts[mut,"Pos"],unique_muts[mut,"Pos"], q=30, mask=3844, mq=30) 
        if(length(ref)==1 & length(alt)==1) {
            ## Variant is a snv
        total_bases = sum(kk[,c("A","C","G","T","a","c","g","t")],na.rm=T)
		total_alt   = sum(kk[,c(alt,alt_lower)],na.rm=T)
        } else {
            # Variant is an indel
        total_bases = sum(kk[,c("A","C","G","T","INS", "DEL","a","c","g","t","ins","del")],na.rm=T)
            if(length(alt)==1) {
                total_alt   = sum(kk[,c(alt,alt_lower)],na.rm=T)
            }

        }

		alts[mut,sample] = total_alt
		covs[mut,sample] = total_bases
	}
}


## Merge DBS mutations



write.table(covs,file=paste(output_dir, patient,"_covs.tsv", sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(alts,file=paste(output_dir, patient,"_alts.tsv", sep=""),sep="\t",row.names=T,col.names=T,quote=F)

