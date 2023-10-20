########################################################################################################################
### quantitative phenotype prediction using data from curatedMetagenomicData
### repeat 10 times each
### 2023/10/20
########################################################################################################################

setwd("/home/wangbb/normalization_quantative_phenotype/curatedMetagenomicData")

all(sapply(c("foreach","doParallel","caret"), require, character.only=TRUE))
source("helper.R")

#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c("data/meta.rds","data/count.rds","TSS","rfr",2)

### parameters
meta <- readRDS(command_args[1])    # metadata
count <- readRDS(command_args[2])   # count table
norm_method <- command_args[3]      # normalization method
pred_method <- command_args[4]      # prediction method
pred_cluster <- as.numeric(command_args[5])  # number of clusters for predictions
studies <- unique(meta$study_name)


#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("norm_data")) dir.create("norm_data")
if(!dir.exists(paste0("norm_data/",norm_method))) dir.create(paste0("norm_data/",norm_method))
for(study1 in studies){
  for(study2 in studies){
    if(study1!=study2){
      # data to be normalized
      meta1 <- meta[meta$study_name==study1,]
      meta2 <- meta[meta$study_name==study2,]
      data1 <- count[,colnames(count)%in%meta1$sample_id]
      data2 <- count[,colnames(count)%in%meta2$sample_id]
      # normalization
      norm_data <- norm.func(p1=data1,p2=data2,norm_method=norm_method)
      # save the results
      saved_file <- paste0("norm_data/",norm_method,"/trn_",study1,"_tst_",study2,"_",norm_method,".rds")
      saveRDS(norm_data, saved_file)
      print(paste0("trn=",study1,",tst=",study2,",norm_method=",norm_method))
    }
  }
}


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("pred_results")) dir.create("pred_results")
if(!dir.exists(paste0("pred_results/",norm_method))) dir.create(paste0("pred_results/",norm_method))
for(study1 in studies){
  for(study2 in studies){
    if(study1!=study2){
      # normalized data
      norm_data <- readRDS(paste0("norm_data/",norm_method,"/trn_",study1,"_tst_",study2,"_",norm_method,".rds"))
      # prediction
      cl<- makeCluster(pred_cluster,setup_strategy="sequential") 
      registerDoParallel(cl)
      pred_res <- foreach(x=1:10,.packages=c("caret","pROC"),.combine=rbind,.errorhandling="pass") %dopar% {
        real.pred.func(trn=norm_data[[1]],tst=norm_data[[2]],meta=meta,pred_method=pred_method)
      }
      stopCluster(cl)
      # summarize the results
      res_df <- data.frame(study1=study1,study2=study2,norm_method=norm_method,pred_method=pred_method,
                           rmse=pred_res[,"rmse"],mape=pred_res[,"mape"])
      saved_file <- paste0("pred_results/",norm_method,"/trn_",study1,"_tst_",study2,"_",norm_method,"_",pred_method,".rds")
      saveRDS(res_df, saved_file)
      print(paste0("study1=",study1,",study2=",study2,",norm_method=",norm_method,",pred_method=",pred_method))
    }
  }
}











