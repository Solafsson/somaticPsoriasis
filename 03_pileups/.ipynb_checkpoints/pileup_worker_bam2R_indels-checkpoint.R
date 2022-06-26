
# bsub -o test.out -e test.err -R"select[mem>2000] rusage[mem=2000]" -M2000 "/software/R-3.6.1/bin/Rscript pileup_worker_bam2R.R"
#.libPaths("/lustre/scratch114/projects/crohns/somatic_ibd_p1/R_packages_farm5_R3.6_install/")
.libPaths("/lustre/scratch119/humgen/projects/psoriasis/R_packages_farm5_R4.1.0_install/")
library(deepSNV)
library(vcfR)
library(VariantAnnotation)

parameters <- commandArgs(TRUE)

patient  <- parameters[1]
sampleVector  <- parameters[2]
project_vector  <- parameters[3]
output_dir  <- parameters[4]
pso_dataDir <- parameters[5]

options(stringsAsFactors=F)

# For testing:
#sampleVector = c("P01H_LL_6","P01L_LL_14", "P03L_T_2", "P05H_T_3", "P05L_T_10", "P06L_T_3", "P09H_UA_1", "P07L_UL_5")
#patient <- "patient01"
#project_vector <- rep("2545", length(sampleVector))
#output_dir="/lustre/scratch119/humgen/projects/psoriasis/pileups/"
#pso_dataDir="/nfs/cancer_ref01/nst_links/live/"
#path=paste("/nfs/cancer_ref01/nst_links/live/", projectNr, "/", sep="")



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
    pindel_file <- paste(path,sample,"/",sample,".pindel.annot.vcf.gz",sep="")
    # pindel_file <- paste(pso_dataDir,sample,".pindel.annot.vcf.gz",sep="")
    vcf <- vcf2tab(caveman_file, sample)
	snvs <- vcf[,c("Chr", "Pos", "Ref", "Alt", "SampleID")]
    vcf <- pindel_vcf_to_tab(pindel_file, sample)
    indels <- vcf[,c("Chr", "Pos", "Ref", "Alt", "Sample")]
    colnames(indels) <- c("Chr", "Pos", "Ref", "Alt", "SampleID")
	if(nrow(all_muts)==0) {
		all_muts = snvs
#        all_muts=indels
        all_muts = rbind(all_muts,indels)
	} else {
		colnames(snvs) = colnames(all_muts)
        colnames(indels) = colnames(snvs)
		all_muts = rbind(all_muts,snvs)
        all_muts = rbind(all_muts,indels)
	}
    i <- i+1
}

cat("Mutations:\n")
print(table(all_muts$SampleID))

unique_muts = all_muts[!duplicated(all_muts[,c("Chr","Pos","Ref","Alt")]),]
cat(nrow(unique_muts), "unique_muts\n")

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

	ref = as.character(unique_muts[mut,"Ref"])
	alt = as.character(unique_muts[mut,"Alt"])
	alt_lower = tolower(alt)
    i <- 1
	for(sample in samples_of_interest) {	
        path=paste(pso_dataDir, project_numbers[i], "/", sep="")
		bam = paste(path,sample,"/",sample,".sample.dupmarked.bam",sep="")

		kk = bam2R(bam,unique_muts[mut,"Chr"],unique_muts[mut,"Pos"],unique_muts[mut,"Pos"], q=30, mask=3844, mq=30, keepflag=3) 
        ##kk = bam2R(bam,unique_muts[mut,"Chr"],unique_muts[mut,"Pos"],unique_muts[mut,"Pos"], q=30, mask=3844, mq=30)
        if(nchar(ref)==1 & nchar(alt)==1) {
            ## Variant is a snv
        total_bases = sum(kk[,c("A","C","G","T","a","c","g","t")],na.rm=T)
		total_alt   = sum(kk[,c(alt,alt_lower)],na.rm=T)
        } else {
            # Variant is an indel
        total_bases = sum(kk[,c("A","C","G","T","INS", "DEL","a","c","g","t","ins","del")],na.rm=T)
        total_alt = sum(kk[,c("INS", "DEL","ins","del")],na.rm=T)

        }

		alts[mut,sample] = total_alt
		covs[mut,sample] = total_bases
        i <- i+1
	}
}


## Merge DBS mutations


## The purpose of this function is to loop over the read-count matricies and merge adjacent SBS that look like 
## they are actually DBS mutations. The rownames of both input objects are mutations on the form chr1:123345:A:C
## and the column names are the sample names. 
## The function looks for instances where adjacent substitutions occur. Once it finds it, it looks at the sample 
## that has the highest coverage of the first substitution and extracts the number of alt/ref reads for that sample.
## It then tests by fisher test if the read counts of both subs are significantly different. If not, it merges them
## into one DBS. 
## The function returns a list the first element of which is the alternative allele read count matrix and the second
## is the coverage matrix. 
merge_sbs_to_dbs <- function(alts, covs) {
  
  test_alts <- data.frame(alts[order(rownames(alts)),])
  test_cov <- data.frame(covs[order(rownames(covs)), ])
  
  test_alts$chr <- unlist(strsplit(rownames(test_alts), split=":"))[c(T,F,F,F)]
  test_alts$pos <- as.numeric(unlist(strsplit(rownames(test_alts), split=":"))[c(F,T,F,F)])
  test_alts$ref <- unlist(strsplit(rownames(test_alts), split=":"))[c(F,F,T,F)]
  test_alts$alt <- unlist(strsplit(rownames(test_alts), split=":"))[c(F,F,F,T)]
  
  test_cov$chr <- unlist(strsplit(rownames(test_cov), split=":"))[c(T,F,F,F)]
  test_cov$pos <- as.numeric(unlist(strsplit(rownames(test_cov), split=":"))[c(F,T,F,F)])
  test_cov$ref <- unlist(strsplit(rownames(test_cov), split=":"))[c(F,F,T,F)]
  test_cov$alt <- unlist(strsplit(rownames(test_cov), split=":"))[c(F,F,F,T)]
  
  candidate_dbs <- data.frame(sbs1=as.character(), sbs2=as.character(), sbs1_alt=as.numeric(),sbs2_alt=as.numeric(), sbs1_ref=as.numeric(),
                              sbs2_ref=as.numeric(), fisher_p=as.numeric(), sample=as.character())
  
  for(i in 1:(nrow(test_cov)-1)) {
    if(test_alts$chr[i]==test_alts$chr[i+1] & test_alts$pos[i]==(test_alts$pos[i+1]-1)) {
      # Subs are adjacent. 
      # Select the sample with the highest number of alt alleles to carry out the test
      sampleNr <- which(test_alts[i,c(1:(ncol(test_alts)-4))]==max(test_alts[i,c(1:(ncol(test_alts)-4))]))[1]
      
      #Create a contingency table for fisher.test
      x <- matrix(c(test_alts[c(i,i+1),sampleNr], c(test_cov[c(i,i+1),sampleNr]-test_alts[c(i,i+1),sampleNr])), nrow=2)
      fishy <- fisher.test(x)
      
      entry <- c(rownames(test_alts)[i], rownames(test_alts)[i+1],test_alts[c(i,i+1),sampleNr], c(test_cov[c(i,i+1),sampleNr]-test_alts[c(i,i+1),sampleNr]), fishy$p.value, colnames(test_cov[sampleNr]))
      candidate_dbs <- rbind(candidate_dbs, entry)
    }
  }
  
  colnames(candidate_dbs) <- c("sbs1", "sbs2","sbs1_alt","sbs2_alt","sbs1_ref","sbs2_ref", "fisher_p", "sample_used_for_testing")
  candidate_dbs$q <- p.adjust(candidate_dbs$fisher_p, method = "hochberg")
  candidate_dbs <- candidate_dbs[candidate_dbs$q > 0.05,]
  
  # Replace the 'first' sbs entry with a dbs entry
  for(i in 1:(nrow(test_cov)-1)) {
    if(rownames(test_cov)[i] %in% candidate_dbs$sbs1) {
      test_cov$ref[i] <- paste(test_cov$ref[i],test_cov$ref[i+1], sep="")
      test_cov$alt[i] <- paste(test_cov$alt[i],test_cov$alt[i+1], sep="")
      rownames(test_cov)[i] <- paste(test_cov$chr[i],test_cov$pos[i], test_cov$ref[i], test_cov$alt[i], sep=":")
      
      test_alts$ref[i] <- paste(test_alts$ref[i],test_alts$ref[i+1], sep="")
      test_alts$alt[i] <- paste(test_alts$alt[i],test_alts$alt[i+1], sep="")
      rownames(test_alts)[i] <- paste(test_alts$chr[i],test_alts$pos[i], test_alts$ref[i], test_alts$alt[i], sep=":")
    } 
  }
  
  # Remove the 'second' sbs entry
  test_cov <- test_cov[!(rownames(test_cov) %in% candidate_dbs$sbs2),]
  test_alts <- test_alts[!(rownames(test_alts) %in% candidate_dbs$sbs2),]
    
  test_cov$chr <- NULL
  test_cov$pos <- NULL
  test_cov$ref <- NULL
  test_cov$alt <- NULL
  
  test_alts$chr <- NULL
  test_alts$pos <- NULL
  test_alts$ref <- NULL
  test_alts$alt <- NULL    
    
  return(list(test_alts,test_cov))
}

dbs_merged <- merge_sbs_to_dbs(alts=alts, covs=covs)
alts <- dbs_merged[[1]]
covs <- dbs_merged[[2]]

write.table(covs,file=paste(output_dir, patient,"/", patient,"_covs.tsv", sep=""),sep="\t",row.names=T,col.names=T,quote=F)
write.table(alts,file=paste(output_dir, patient,"/", patient,"_alts.tsv", sep=""),sep="\t",row.names=T,col.names=T,quote=F)

