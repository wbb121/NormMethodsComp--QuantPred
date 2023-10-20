########################################################################################################################
### summarize the prediction results and plot the figures for simulations and real data (curatedMetagenomicData)
### 2023/10/20
########################################################################################################################

setwd("/Users/beibeiwang/Desktop/Normalization_quantitative/")

### packages
all(sapply(c("vegan","ape","reshape2","ggplot2","ggtree","viridis","RColorBrewer","ggsci","aplot","ggpubr",
             "grid","ggplotify","patchwork"),require, character.only=TRUE))

#======================================================================================================================#
### normalization methods ###
#======================================================================================================================#
norm_methods <- data.frame(
  norm_method=c("TSS","UQ","MED","CSS","TMM","RLE_poscounts","GMPR","CLR+",
                "logcpm","LOG","AST","STD","rank","blom","VST","NPN",
                "QN","FSQN","BMC","limma","combat","conqur"),
  annotated_method=c("TSS","UQ","MED","CSS","TMM","RLE","GMPR","CLR",
                     "logCPM","LOG","AST","STD","Rank","Blom","VST","NPN",
                     "QN","FSQN","BMC","Limma","ComBat","ConQuR"),
  class=c(rep("Scaling",7),rep("CoDA",1),rep("Transformation",8),rep("Batch Correction",6))
)
norm_methods$annotated_method <- factor(norm_methods$annotated_method,
                                        levels=c("TSS","UQ","MED","CSS","TMM","RLE","GMPR","CLR",
                                                 "LOG","AST","STD","Rank","Blom","NPN","logCPM","VST",
                                                 "QN","FSQN","BMC","Limma","ComBat","ConQuR"))
norm_methods$class <- factor(norm_methods$class,levels=c("Scaling","CoDA","Transformation","Batch Correction"))
saveRDS(norm_methods,"data/norm_methods.rds")


#======================================================================================================================#
### simulation scenario1 ###
#======================================================================================================================#
### parameters
alphas <- c(0,0.2,0.4,0.6,0.8,1)
pred_method <- "rfr"
phenotypes <- c("linear","quadratic","inverse","logistic")
norm_methods <- readRDS("data/norm_methods.rds")

### prediction results
scenario1_df <- data.frame(alpha=character(0),norm_method=character(0),pred_method=character(0),phenotype=character(0),
                           rmse=numeric(0),mape=numeric(0))
scenario1_summ_df <- data.frame(alpha=character(0),norm_method=character(0),pred_method=character(0),phenotype=character(0),
                                average_rmse=numeric(0),average_mape=numeric(0),median_rmse=numeric(0),median_mape=numeric(0))
for(alpha in alphas){
  for(phenotype in phenotypes){
    for(norm_method in norm_methods$norm_method){
      #alpha=0;phenotype="linear";norm_method="TSS"
      df <- readRDS(paste0("scenario1/pred_results/",norm_method,"/norm_alpha",alpha,"_",norm_method,"_",pred_method,"_",phenotype,".rds"))
      summ_df <- data.frame(alpha=alpha,norm_method=norm_method,pred_method=pred_method,phenotype=phenotype,
                            average_rmse=round(mean(as.numeric(df$rmse)),3),average_mape=round(mean(as.numeric(df$mape)),3),
                            median_rmse=round(median(as.numeric(df$rmse)),3),median_mape=round(median(as.numeric(df$mape)),3))
      scenario1_df <- rbind(scenario1_df,df)
      scenario1_summ_df <- rbind(scenario1_summ_df,summ_df)
    }
  }
}
saveRDS(scenario1_df,"results/scenario1_res.rds")
saveRDS(scenario1_summ_df,"results/scenario1_res_summ.rds")


#======================================================================================================================#
### simulation scenario2 ###
#======================================================================================================================#
### parameters
parameters <- list(c(0,1),c(500,1),c(1000,1),c(0,2),c(0,4))
pred_method <- "rfr"
phenotypes <- c("linear","quadratic","inverse","logistic")
norm_methods <- readRDS("data/norm_methods.rds")


### prediction results
scenario2_df <- data.frame(batch_mean=numeric(0),batch_var=numeric(0),norm_method=character(0),
                           pred_method=character(0),phenotype=character(0),rmse=numeric(0),mape=numeric(0))
scenario2_summ_df <- data.frame(batch_mean=numeric(0),batch_var=numeric(0),norm_method=character(0),
                                pred_method=character(0),phenotype=character(0),
                                average_rmse=numeric(0),average_mape=numeric(0),
                                median_rmse=numeric(0),median_mape=numeric(0))
for(i in 1:length(parameters)){
  for(phenotype in phenotypes){
    for(norm_method in norm_methods$norm_method){
      #alpha=0;phenotype="linear";norm_method="TSS"
      df <- readRDS(paste0("scenario2/pred_results/",norm_method,"/pred_",pred_method,"_mean",parameters[[i]][1],"_var",parameters[[i]][2],"_",phenotype,".rds"))
      summ_df <- data.frame(batch_mean=parameters[[i]][1],batch_var=parameters[[i]][2],norm_method=norm_method,pred_method=pred_method,phenotype=phenotype,
                            average_rmse=round(mean(as.numeric(df$rmse)),3),average_mape=round(mean(as.numeric(df$mape)),3),
                            median_rmse=round(median(as.numeric(df$rmse)),3),median_mape=round(median(as.numeric(df$mape)),3))
      scenario2_df <- rbind(scenario2_df,df)
      scenario2_summ_df <- rbind(scenario2_summ_df,summ_df)
    }
  }
}
saveRDS(scenario2_df,"results/scenario2_res.rds")
saveRDS(scenario2_summ_df,"results/scenario2_res_summ.rds")


#======================================================================================================================#
### simulation scenario3 ###
#======================================================================================================================#
### parameters
overlaps <- c(2,4,6,8,10)
pred_method <- "rfr"
phenotypes <- c("linear","quadratic","inverse","logistic")
norm_methods <- readRDS("data/norm_methods.rds")

### prediction results
scenario3_df <- data.frame(overlap=character(0),norm_method=character(0),pred_method=character(0),phenotype=character(0),
                           rmse=numeric(0),mape=numeric(0))
scenario3_summ_df <- data.frame(overlap=character(0),norm_method=character(0),pred_method=character(0),phenotype=character(0),
                                average_rmse=numeric(0),average_mape=numeric(0),median_rmse=numeric(0),median_mape=numeric(0))
for(overlap in overlaps){
  for(phenotype in phenotypes){
    for(norm_method in norm_methods$norm_method){
      df <- readRDS(paste0("scenario3/pred_results/",norm_method,"/norm_overlap",overlap,"_",norm_method,"_",pred_method,"_",phenotype,".rds"))
      summ_df <- data.frame(overlap=overlap,norm_method=norm_method,pred_method=pred_method,phenotype=phenotype,
                            average_rmse=round(mean(as.numeric(df$rmse)),3),average_mape=round(mean(as.numeric(df$mape)),3),
                            median_rmse=round(median(as.numeric(df$rmse)),3),median_mape=round(median(as.numeric(df$mape)),3))
      scenario3_df <- rbind(scenario3_df,df)
      scenario3_summ_df <- rbind(scenario3_summ_df,summ_df)
    }
  }
}
saveRDS(scenario3_df,"results/scenario3_res.rds")
saveRDS(scenario3_summ_df,"results/scenario3_res_summ.rds")


#======================================================================================================================#
### curatedMetagenomicData ###
#======================================================================================================================#
### parameters
meta <- readRDS("data/meta.rds")
datasets <- unique(meta$study_name)
pred_method <- "rfr"
norm_methods <- readRDS("data/norm_methods.rds")


### prediction results
res_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),rmse=numeric(0),mape=numeric(0))
res_summ_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),average_rmse=numeric(0),average_mape=numeric(0))
for(study1 in datasets){
  for(study2 in datasets){
    for(norm_method in norm_methods$norm_method){
      if(study1 != study2){
        df <- readRDS(paste0("curatedMetagenomicData/pred_results/",norm_method,"/trn_",study1,"_tst_",study2,"_",norm_method,"_",pred_method,".rds"))
        summ_df <- data.frame(study1=study1,study2=study2,norm_method=norm_method,pred_method=pred_method,
                              average_rmse=round(mean(as.numeric(df$rmse)),3),average_mape=round(mean(as.numeric(df$mape)),3))
        res_df <- rbind(res_df,df)
        res_summ_df <- rbind(res_summ_df,summ_df)
      }
    }
  }
}
saveRDS(res_df,"results/curatedmetagenomicData_res.rds")
saveRDS(res_summ_df,"results/curatedmetagenomicData_res_summ.rds")


### order the rmse and mape for different normalization method
res_summ_order_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                                average_rmse=numeric(0),average_mape=numeric(0),rank_average_rmse=integer(0),rank_average_mape=integer(0))
for(study1 in datasets){
  for(study2 in datasets){
    if(study1 != study2){
      df <- res_summ_df[res_summ_df$study1==study1 & res_summ_df$study2==study2,]
      df <- df %>% mutate(rank_average_rmse=rank(average_rmse,ties.method="min"))
      df <- df %>% mutate(rank_average_mape=rank(average_mape,ties.method="min"))
      res_summ_order_df <- rbind(res_summ_order_df,df)
    }
  }
}
saveRDS(res_summ_order_df,"results/curatedmetagenomicData_res_summ_order.rds")



