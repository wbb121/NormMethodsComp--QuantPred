########################################################################################################################
### obtain available data from curatedMetagenomicData
### 2023/10/20
########################################################################################################################

setwd("/Users/beibeiwang/Desktop/Normalization_quantitative/")

### packages
#BiocManager::install("curatedMetagenomicData")
all(sapply(c("curatedMetagenomicData","dplyr"),require, character.only=TRUE))

#======================================================================================================================#
### meta, stool samples, healthy status, known BMI, no duplicated samples belonging to the same subject, reads num > 1250
meta <- sampleMetadata %>%
  filter(body_site=="stool") %>%
  filter(disease == "healthy") %>%
  filter(!is.na(BMI)) %>%
  filter(!duplicated(subject_id)) %>%
  filter(number_reads>1250) 


### count, stool samples, healthy status, known BMI, no duplicated samples belonging to the same subject
# no available count data for MetaCardis_2020_a
count <- sampleMetadata %>%
  filter(body_site=="stool") %>%
  filter(disease == "healthy") %>%
  filter(!is.na(BMI)) %>%
  filter(!duplicated(subject_id)) %>%
  filter(number_reads>1250) %>%
  filter(study_name!="MetaCardis_2020_a") %>%
  returnSamples(dataType="relative_abundance",counts=T)
count_lst <- assays(count)
count_tab <- as.data.frame(count_lst$relative_abundance)


#======================================================================================================================#
### modify the metadata
# no available count data for MetaCardis_2020_a
meta <- meta %>% filter(study_name!="MetaCardis_2020_a") 
# remove samples not having relative abundance
meta <- meta[meta$sample_id%in%colnames(count_tab),]
# YachidaS_2019 dataset contains the ThomasAM_2019_c data, remove ThomasAM_2019_c
meta <- meta[meta$study_name!="ThomasAM_2019_c",]
# combine ThomasAM_2018a and ThomasAM_2018b to ThomasAM_2018
meta$study_name <- gsub("ThomasAM_2018a","ThomasAM_2018",meta$study_name)
meta$study_name <- gsub("ThomasAM_2018b","ThomasAM_2018",meta$study_name)
# remove datasets with sample size < 30
sample_sizes <- as.data.frame(table(meta$study_name))
colnames(sample_sizes) <- c("study_name","sample_size")
studies <- as.character(sample_sizes[sample_sizes$sample_size>=30,"study_name"])
meta <- meta[meta$study_name%in%studies,]
saveRDS(meta,"data/meta.rds")


### modify the count table according to the metadata
count_tab <- count_tab[,colnames(count_tab)%in%meta$sample_id]
count_tab <- count_tab[rowSums(count_tab)>0,]
saveRDS(count_tab,"data/count.rds")


#======================================================================================================================#
### summarize the profiles of different datasets
datasets <- unique(meta$study_name)
profiles <- data.frame(Dataset=character(0),Country=character(0),Sample_size=numeric(0),DNA_Exk=character(0),Seq_Plat=character(0),Reference=character(0))
for(i in 1:length(datasets)){
  df <- meta[meta$study_name==datasets[i],]
  profiles[i,"Dataset"] <- datasets[i]
  profiles[i,"Country"] <- paste(unique(df$country),collapse="/")
  profiles[i,"Sample_size"] <- sum(df$disease=="healthy")
  profiles[i,"DNA_Exk"] <- paste(unique(df$DNA_extraction_kit),collapse="/")
  profiles[i,"Seq_Plat"] <- paste(unique(df$sequencing_platform),collapse="/")
  profiles[i,"Reference"] <- paste(unique(df$PMID),collapse="/")
}
# DhakanDB_2019 and GuptaA_2019 are from the same paper, check the original paper and found healthy controls from GuptaA_2019 are overlapped with DhakanDB_2019
# remove GuptaA_2019 from analysis
meta <- meta[meta$study_name!="GuptaA_2019",]
count <- count[,colnames(count)%in%meta$sample_id]
count <- count[rowSums(count)>0,]
profiles <- profiles[profiles$Dataset!="GuptaA_2019",]
# save the data
saveRDS(meta,"data/meta.rds")
saveRDS(count,"data/count.rds")
write.csv(profiles,"results/curatedMetagenomicData_profiles.csv")


