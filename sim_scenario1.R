########################################################################################################################
### simulation in scenario 1
### run the code in SYS-7089P-TR4T sever
### 2023/10/20
########################################################################################################################

setwd("/home/wangbb/normalization_quantative_phenotype/scenario1")

suppressPackageStartupMessages(library(DirichletReg))
suppressPackageStartupMessages(library(caret))
source("helper.R")

#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c("TSS","rfr")

### parameters
alphas=c(0,0.2,0.4,0.6,0.8,1)               # population effects
sample_size=100                             # sample size
library_size=1000000                        # library size
num_genes=10                                # number of phenotype related gene 
count1=readRDS("data/GuptaA_2019_ctrl_count.rds")            # template data1
count2=readRDS("data/FengQ_2015_ctrl_count.rds")             # template data2
phenotypes=c("linear","quadratic","inverse","logistic")      # relationships between quantitative phenotype and phenotype-associated taxa
norm_method=command_args[1]                 # normalization method
pred_method=command_args[2]                 # prediction method

# make the names for genes in two populations
rownames(count1) <- make.names(rownames(count1))
rownames(count2) <- make.names(rownames(count2))
# choose phenotype related gene
union_genes <- union(rownames(count1),rownames(count2))
inter_genes <- intersect(rownames(count1),rownames(count2))
set.seed(1234); selected_genes <- sample(inter_genes,num_genes,replace=F)
# coefficients for selected genes
set.seed(1234); coefficients <- c(runif(num_genes-round(num_genes/2),3,5), runif(round(num_genes/2),-5,-3))


#======================================================================================================================#
### functions ###
#======================================================================================================================#
### function for simulating the count table
sim.count.func <- function(count1,count2,sample_size,library_size,alpha,seed){
  
  # the probability vectors from real data
  tab1 <- merge.func(count1,count2)[,colnames(count1)]
  tab2 <- merge.func(count1,count2)[,colnames(count2)]
  prob1 <- rowSums(tab1)/sum(rowSums(tab1))
  prob2 <- rowSums(tab2)/sum(rowSums(tab2))
  
  # pseudo probability vectors with population effect alpha
  prob1 <- alpha*prob1+(1-alpha)*prob2
  prob1 <- prob1/sum(prob1)
  prob1 <- prob1[prob1!=0]
  prob2 <- prob2[prob2!=0]
  
  # generate the dirichlet distributed probability vectors for each sample based on the base probability vectors
  set.seed(seed); prob1s <- rdirichlet(sample_size,1e6*prob1)
  set.seed(seed); prob2s <- rdirichlet(sample_size,1e6*prob2)
  
  # simulate the counts table using multinomial distribution MN(library_size,prob)
  sim1 <- matrix(NA,nrow=length(prob1),ncol=sample_size,dimnames=list(names(prob1),paste0("sim1_",1:sample_size)))
  sim2 <- matrix(NA,nrow=length(prob2),ncol=sample_size,dimnames=list(names(prob2),paste0("sim2_",1:sample_size)))
  for(i in 1:sample_size){
    set.seed(i); sim1[,i] <- rmultinom(1, size=library_size, prob=prob1s[i,])
    set.seed(i); sim2[,i] <- rmultinom(1, size=library_size, prob=prob2s[i,])
  }
  
  # return
  return(list(sim1,sim2))
  
}


### function for simulating the quantitative traits
sim.phenotype.func <- function(count_table,selected_genes,coefficients){
  
  # process the count data
  count_table <- apply(count_table,2,function(x) x/sum(x))   # normalization - tss
  count_table <- t(count_table)   # transformation
  
  # data frame to save the simulated phenotype
  phenotype_df <- data.frame(sample_id=rownames(count_table),linear=NA,quadratic=NA,inverse=NA,logistic=NA)
  
  # simulate the quantitative phenotype based on different relationship between phenotype and related genes
  noise <- rnorm(sample_size)
  component <- as.matrix(count_table[,selected_genes])%*%coefficients
  component_quadratic <- as.matrix(apply(count_table[,selected_genes],c(1,2),function(x) x^2))%*%coefficients
  phenotype_df[,"linear"] <- 1e2*component+noise
  phenotype_df[,"quadratic"] <- 1e4*component_quadratic+noise
  phenotype_df[,"inverse"] <- 10/component+noise
  phenotype_df[,"logistic"] <- 1e2/(1+exp(component))+noise
  
  # round two decimal places
  phenotype_df$linear <- round(phenotype_df$linear,2)
  phenotype_df$quadratic <- round(phenotype_df$quadratic,2)
  phenotype_df$inverse <- round(phenotype_df$inverse,2)
  phenotype_df$logistic <- round(phenotype_df$logistic,2)
  
  # return
  return(phenotype_df)
  
}


#======================================================================================================================#
### simulation ###
#======================================================================================================================#
if(!dir.exists("scenario1/sim_data")) dir.create("scenario1/sim_data")
for(alpha in alphas){
  sim_tabs <- list()
  sim_meta <- list()
  for(i in 1:100){
    sim_tabs[[i]] <- sim.count.func(count1=count1,count2=count2,sample_size=sample_size,library_size=library_size,
                                    num_genes=num_genes,alpha=alpha,seed=i)
    sim_meta[[i]] <- list(sim.phenotype.func(count_table=sim_tabs[[i]][[1]],selected_genes=selected_genes,coefficients=coefficients),
                          sim.phenotype.func(count_table=sim_tabs[[i]][[2]],selected_genes=selected_genes,coefficients=coefficients))
  }
  saveRDS(sim_tabs, paste0("scenario1/sim_data/sim_alpha",alpha,".rds"))
  saveRDS(sim_meta, paste0("scenario1/sim_data/sim_alpha",alpha,"_meta.rds"))
  print(paste0("alpha=",alpha))
}
# Simulated data only needs to be generated once. 
# Once it's generated, you can comment out the simulation section of the code.


#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("scenario1/norm_data")) dir.create("scenario1/norm_data")
if(!dir.exists(paste0("scenario1/norm_data/",norm_method))) dir.create(paste0("scenario1/norm_data/",norm_method))
for(alpha in alphas){
  sim_tabs <- readRDS(paste0("scenario1/sim_data/sim_alpha",alpha,".rds"))
  norm_tabs <- list()
  for(i in 1:100){
    norm_tabs[[i]] <-  norm.func(p1=sim_tabs[[i]][[1]],p2=sim_tabs[[i]][[2]],norm_method=norm_method)
  }
  saveRDS(norm_tabs, paste0("scenario1/norm_data/",norm_method,"/norm_alpha",alpha,"_",norm_method,".rds"))
  print(paste0("alpha=",alpha,",norm_method=",norm_method))
}


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("scenario1/pred_results")) dir.create("scenario1/pred_results")
if(!dir.exists(paste0("scenario1/pred_results/",norm_method))) dir.create(paste0("scenario1/pred_results/",norm_method))
for(alpha in alphas){
  for(phenotype in phenotypes){
    norm_tabs=readRDS(paste0("scenario1/norm_data/",norm_method,"/norm_alpha",alpha,"_",norm_method,".rds"))
    sim_meta=readRDS(paste0("scenario1/sim_data/sim_alpha",alpha,"_meta.rds"))
    pred_df <- data.frame(alpha=numeric(0),norm_method=character(0),pred_method=character(0),phenotype=character(0),
                          rmse=numeric(0),mape=numeric(0))
    for(i in 1:100){
      pred_res <- sim.pred.func(trn=norm_tabs[[i]][[1]],tst=norm_tabs[[i]][[2]],
                                trn_meta=sim_meta[[i]][[1]],tst_meta=sim_meta[[i]][[2]],
                                phenotype=phenotype,pred_method=pred_method)
      pred_df[i,] <- c(alpha=alpha,norm_method=norm_method,pred_method=pred_method,phenotype=phenotype,
                       rmse=round(pred_res["rmse"],3),mape=round(pred_res["mape"],3))
    }
    saveRDS(pred_df, paste0("scenario1/pred_results/",norm_method,"/norm_alpha",alpha,"_",norm_method,"_",pred_method,"_",phenotype,".rds"))
    print(paste0("alpha=",alpha,",norm_method=",norm_method,",pred_method=",pred_method,",phenotype=",phenotype))
  }
}






