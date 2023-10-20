########################################################################################################################
### simulation in scenario 3
### run the code in SYS-7089P-TR4T sever
### 2023/10/20
########################################################################################################################

setwd("/home/wangbb_gpu/normalization_quantative_phenotype/scenario3")

suppressPackageStartupMessages(library(DirichletReg))
suppressPackageStartupMessages(library(caret))
source("helper.R")


#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c("TSS","rfr")

# parameters
num_genes_overlap=c(2,4,6,8,10)             # number of overlapped genes
sample_size=100                             # sample size
library_size=1000000                        # library size
num_genes=10                                # number of phenotype related gene 
count=readRDS("data/FengQ_2015_ctrl_count.rds")               # template data
phenotypes=c("linear","quadratic","inverse","logistic")       # relationships between quantitative phenotype and phenotype-associated taxa
norm_method=command_args[1]                 # normalization method
pred_method=command_args[2]                 # prediction method

# make the names for genes in two populations
rownames(count) <- make.names(rownames(count))

# choose phenotype related gene
set.seed(1234); selected_genes <- sample(rownames(count),num_genes,replace=F)
set.seed(1234); extra_genes <- sample(rownames(count)[!rownames(count)%in%selected_genes],num_genes,replace=F)
selected_genes_df <- as.data.frame(matrix(NA,nrow=length(num_genes_overlap),ncol=num_genes,
                                          dimnames=list(paste0("overlap",num_genes_overlap),paste0("genes",1:num_genes))))
for(i in 1:length(num_genes_overlap)){
  if(num_genes_overlap[i]<num_genes) 
    selected_genes_df[i,] <- c(selected_genes[1:(num_genes_overlap[i]/2)],
                               extra_genes[1:((num_genes-num_genes_overlap[i])/2)],
                               selected_genes[(num_genes/2+1):(num_genes/2+num_genes_overlap[i]/2)],
                               extra_genes[((num_genes-2)/2+1):((num_genes-2)/2+(num_genes-num_genes_overlap[i])/2)])
  else selected_genes_df[i,] <- selected_genes
}

# coefficients for selected genes
set.seed(1234); coefficients <- c(runif(num_genes/2,3,5), runif((num_genes/2),-5,-3))
set.seed(12345); extra_coefficients <- c(runif((num_genes-2)/2,3,5), runif((num_genes-2)/2,-5,-3))
coefficients_df <- as.data.frame(matrix(NA,nrow=length(num_genes_overlap),ncol=num_genes,
                                        dimnames=list(paste0("overlap",num_genes_overlap),paste0("coef",1:num_genes))))
for(i in 1:length(num_genes_overlap)){
  if(num_genes_overlap[i]<num_genes) 
    coefficients_df[i,] <- c(coefficients[1:(num_genes_overlap[i]/2)],
                             extra_coefficients[1:((num_genes-num_genes_overlap[i])/2)],
                             coefficients[(num_genes/2+1):(num_genes/2+num_genes_overlap[i]/2)],
                             extra_coefficients[((num_genes-2)/2+1):((num_genes-2)/2+(num_genes-num_genes_overlap[i])/2)])
  else coefficients_df[i,] <- coefficients
}


#======================================================================================================================#
### functions ###
#======================================================================================================================#
### function for simulating count table
sim.count.func <- function(count,sample_size,library_size,seed){
  
  # generate the probability vectors for each sample based on the base probability vectors
  prob <- rowSums(count)/sum(rowSums(count))
  set.seed(seed); prob1s <- rdirichlet(sample_size,1e6*prob)
  set.seed(seed+100); prob2s <- rdirichlet(sample_size,1e6*prob)
  
  # simulate the counts table using multinomial distribution MN(library_size,prob)
  sim1 <- matrix(NA,nrow=length(prob),ncol=sample_size,dimnames=list(names(prob),paste0("sim1_",1:sample_size)))
  sim2 <- matrix(NA,nrow=length(prob),ncol=sample_size,dimnames=list(names(prob),paste0("sim2_",1:sample_size)))
  for(i in 1:sample_size){
    set.seed(i); sim1[,i] <- rmultinom(1, size=library_size, prob=prob1s[i,])
    set.seed(i); sim2[,i] <- rmultinom(1, size=library_size, prob=prob2s[i,])
  }
  
  # return
  return(list(sim1,sim2))
  
}


### function for simulate the quantitative phenotype
sim.phenotype.func <- function(count_table,selected_genes,coefficients){
  
  # process the count data
  #count_table=sim.count.func(count=count,sample_size=sample_size,library_size=library_size,seed=1)[[1]]
  count_table <- apply(count_table,2,function(x) x/sum(x))   # normalization - tss
  count_table <- t(count_table)   # transformation
  
  # data frame to save the simulated phenotype
  phenotype_df <- data.frame(sample_id=rownames(count_table),linear=NA,quadratic=NA,inverse=NA,logistic=NA)
  
  # simulate the quantitative phenotype based on different relationship between phenotype and related genes
  noise <- rnorm(sample_size)
  component <- as.matrix(count_table[,selected_genes])%*%coefficients
  component_quadratic <- as.matrix(apply(count_table[,selected_genes],c(1,2),function(x) x^2))%*%coefficients
  phenotype_df[,"linear"] <- 1e4*component+noise
  phenotype_df[,"quadratic"] <- 1e6*component_quadratic+noise
  phenotype_df[,"inverse"] <- 1/component+noise
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
if(!dir.exists("scenario3")) dir.create("scenario3")
if(!dir.exists("scenario3/sim_data")) dir.create("scenario3/sim_data")
# simulate the count table
sim_tabs <- list()
for(i in 1:100){
  sim_tabs[[i]] <- sim.count.func(count=count,sample_size=sample_size,library_size=library_size,seed=i)
}
saveRDS(sim_tabs,"scenario3/sim_data/sim_tabs.rds")

# simulate the phenotypes
for(j in 1:length(num_genes_overlap)){
  sim_meta <- list()
  for(i in 1:100){
    sim_meta[[i]] <- list(sim.phenotype.func(count_table=sim_tabs[[i]][[1]],
                                             selected_genes=as.character(selected_genes_df[length(num_genes_overlap),]),
                                             coefficients=as.numeric(coefficients_df[length(num_genes_overlap),])),
                          sim.phenotype.func(count_table=sim_tabs[[i]][[2]],
                                             selected_genes=as.character(selected_genes_df[j,]),
                                             coefficients=as.numeric(coefficients_df[j,])))
  }
  saveRDS(sim_meta,paste0("scenario3/sim_data/sim_meta_overlap",num_genes_overlap[j],".rds"))
  print(paste0("overlap_genes=",num_genes_overlap[j]))
}
# Simulated data only needs to be generated once. 
# Once it's generated, you can comment out the simulation section of the code.


#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("scenario3/norm_data")) dir.create("scenario3/norm_data")
if(!dir.exists(paste0("scenario3/norm_data/",norm_method))) dir.create(paste0("scenario3/norm_data/",norm_method))
sim_tabs <- readRDS("scenario3/sim_data/sim_tabs.rds")
norm_tabs <- list()
for(i in 1:100){
  norm_tabs[[i]] <-  norm.func(p1=sim_tabs[[i]][[1]],p2=sim_tabs[[i]][[2]],norm_method=norm_method)
}
saveRDS(norm_tabs, paste0("scenario3/norm_data/",norm_method,"/norm_",norm_method,".rds"))
print(paste0("norm_method=",norm_method))


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("scenario3/pred_results")) dir.create("scenario3/pred_results")
if(!dir.exists(paste0("scenario3/pred_results/",norm_method))) dir.create(paste0("scenario3/pred_results/",norm_method))
phenotypes=c("linear","quadratic","inverse","logistic")
for(j in 1:length(num_genes_overlap)){
  for(phenotype in phenotypes){
    norm_tabs=readRDS(paste0("scenario3/norm_data/",norm_method,"/norm_",norm_method,".rds"))
    sim_meta=readRDS(paste0("scenario3/sim_data/sim_meta_overlap",num_genes_overlap[j],".rds"))
    pred_df <- data.frame(overlap_geens=numeric(0),norm_method=character(0),pred_method=character(0),
                          phenotype=character(0),rmse=numeric(0),mape=numeric(0))
    for(i in 1:100){
      pred_res <- sim.pred.func(trn=norm_tabs[[i]][[1]],tst=norm_tabs[[i]][[2]],
                                trn_meta=sim_meta[[i]][[1]],tst_meta=sim_meta[[i]][[2]],
                                phenotype=phenotype,pred_method=pred_method)
      pred_df[i,] <- c(overlap_geens=num_genes_overlap[j],norm_method=norm_method,pred_method=pred_method,phenotype=phenotype,
                       rmse=round(pred_res["rmse"],3),mape=round(pred_res["mape"],3))
    }
    saveRDS(pred_df, paste0("scenario3/pred_results/",norm_method,"/norm_overlap",num_genes_overlap[j],"_",norm_method,"_",pred_method,"_",phenotype,".rds"))
    print(paste0("overlap_genes=",num_genes_overlap[j],",norm_method=",norm_method,",pred_method=",pred_method,",phenotype=",phenotype))
  }
}




