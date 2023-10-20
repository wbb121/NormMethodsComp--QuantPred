########################################################################################################################
### simulation in scenario 2
### run the code in SYS-7089P-TR4T sever
### 2023/10/20
########################################################################################################################

setwd("/home/wangbb/normalization_quantative_phenotype/scenario2")

suppressPackageStartupMessages(library(DirichletReg))
suppressPackageStartupMessages(library(MCMCpack))
suppressPackageStartupMessages(library(caret))
source("helper.R")


#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c("TSS","rfr")

# parameters
parameters <- list(c(0,1),c(500,1),c(1000,1),c(0,2),c(0,4))   # combinations of batch_mean and batch_var
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
# coefficients for selected genes
set.seed(1234); coefficients <- c(runif(num_genes-round(num_genes/2),3,5), runif(round(num_genes/2),-5,-3))


#======================================================================================================================#
### functions ###
#======================================================================================================================#
### functions related to inverse gamma distribution
# transform the mean and variance to scale and shape parameters
mv2ab <- function(m, v){
  a <- 2 + m^2/v
  b <- m * (a-1)
  return(list(alpha=a, beta=b))
}
# transform the scale and shape parameters to mean and variance
ab2mv <- function(a, b){
  m <- b / (a-1)
  v <- b^2 / ((a-1)^2*(a-2))
  return(list(mean=m, var=v))
}


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


### function for simulating batch effect
sim.batch.func <- function(count_table,meta,phenotype,hyper_pars,seed){
  
  #count_table=sim.count.func(count=count,sample_size=sample_size,library_size=library_size,seed=seed)[[1]];
  #metadata=sim.phenotype.func(count_table=count_table,selected_genes=selected_genes,coefficients=coefficients);phenotype="linear"
  
  # Simulate batch parameters from hyper-pars
  set.seed(seed); gamma <- rnorm(nrow(count_table), mean=hyper_pars$hyper_mu, sd=hyper_pars$hyper_sd)
  set.seed(seed); delta2 <- rinvgamma(nrow(count_table), shape=hyper_pars$hyper_alpha, scale=hyper_pars$hyper_beta)
  
  # Simulate batch effect
  # fit linear model to data with no batch parameters, calculate residual variance
  X <- model.matrix(~Condition, data=data.frame(Condition=meta[,phenotype]))
  beta <- solve(t(X) %*% X) %*% t(X) %*% t(count_table)
  resid <- count_table - t(X %*% beta)
  # spike-in batch variance: multiply by condition adjusted data with delta
  resid_varbatch <- resid*sqrt(delta2)
  # construct mean batch parameter design matrix for adding gamma
  X_batch <- model.matrix(~-1+Batch, data=data.frame(Batch=rep(1,ncol(count_table))))
  # new data with added batch effect
  new_count_table <- t(cbind(X, X_batch) %*% rbind(beta, gamma)) + resid_varbatch
  # round the new data, making it integers
  new_count_table <- round(new_count_table,0)
  # if the count<0, set it 0
  new_count_table[new_count_table<0] <- 0
  
  # return
  return(new_count_table)
  
}


#======================================================================================================================#
### simulation ###
#======================================================================================================================#
if(!dir.exists("scenario2/sim_data")) dir.create("scenario2/sim_data")
sim_tabs <- list()
sim_meta <- list()
for(i in 1:100){
  # simulate the count table and metadata
  sim_tabs[[i]] <- sim.count.func(count=count,sample_size=sample_size,library_size=library_size,seed=i)
  sim_meta[[i]] <- list(sim.phenotype.func(count_table=sim_tabs[[i]][[1]],selected_genes=selected_genes,coefficients=coefficients),
                        sim.phenotype.func(count_table=sim_tabs[[i]][[2]],selected_genes=selected_genes,coefficients=coefficients))
}
saveRDS(sim_tabs,"scenario2/sim_data/sim_tabs.rds")
saveRDS(sim_meta,"scenario2/sim_data/sim_meta.rds")

# simulate the batch effect
sim_tabs <- readRDS("scenario2/sim_data/sim_tabs.rds")
sim_meta <- readRDS("scenario2/sim_data/sim_meta.rds")
phenotypes <- c("linear","quadratic","inverse","logistic")
for(parameter in parameters){
  batch_mean <- parameter[1]
  batch_var <- parameter[2]
  # hyper-parameters for simulating the batch effect
  # the mean of gene g caused by batch: gamma_g~N(mu,sd2)
  # the variance of gene g caused by batch: delta_g~InvGamma(alpha,beta)
  hyper_pars <- list(hyper_mu=batch_mean, hyper_sd=sqrt(0.01),
                     hyper_alpha=mv2ab(m=batch_var,v=0.01)$alpha,
                     hyper_beta=mv2ab(m=batch_var,v=0.01)$beta)
  # simulate the batch effect
  for(phenotype in phenotypes){
    sim_tabs_adjust <- list()
    for(i in 1:100){
      sim_tabs_adjust_trn <- sim.batch.func(count_table=sim_tabs[[i]][[1]],meta=sim_meta[[i]][[1]],
                                            phenotype=phenotype,hyper_pars=hyper_pars,seed=i)
      sim_tabs_adjust_tst <- sim_tabs[[i]][[2]]
      sim_tabs_adjust[[i]] <- list(sim_tabs_adjust_trn,sim_tabs_adjust_tst)
    }
    saveRDS(sim_tabs_adjust,paste0("scenario2/sim_data/sim_tabs_mean",batch_mean,"_var",batch_var,"_",phenotype,".rds"))
    print(paste0("batch_mean=",batch_mean,",batch_var=",batch_var,",phenotype=",phenotype))
  }
}
# Simulated data only needs to be generated once. 
# Once it's generated, you can comment out the simulation section of the code.


#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("scenario2/norm_data")) dir.create("scenario2/norm_data")
if(!dir.exists(paste0("scenario2/norm_data/",norm_method))) dir.create(paste0("scenario2/norm_data/",norm_method))
for(parameter in parameters){
  batch_mean <- parameter[1]
  batch_var <- parameter[2]
  for(phenotype in phenotypes){
    sim_tabs <- readRDS(paste0("scenario2/sim_data/sim_tabs_mean",batch_mean,"_var",batch_var,"_",phenotype,".rds"))
    norm_tabs <- list()
    for(i in 1:100){
      norm_tabs[[i]] <-  norm.func(p1=sim_tabs[[i]][[1]],p2=sim_tabs[[i]][[2]],norm_method=norm_method)
    }
    saveRDS(norm_tabs, paste0("scenario2/norm_data/",norm_method,"/norm_mean",batch_mean,"_var",batch_var,"_",phenotype,".rds"))
    print(paste0("batch_mean=",batch_mean,",batch_var=",batch_var,",phenotype=",phenotype,",norm_method=",norm_method))
  }
}


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("scenario2/pred_results")) dir.create("scenario2/pred_results")
if(!dir.exists(paste0("scenario2/pred_results/",norm_method))) dir.create(paste0("scenario2/pred_results/",norm_method))
for(parameter in parameters){
  batch_mean <- parameter[1]
  batch_var <- parameter[2]
  for(phenotype in phenotypes){
    norm_tabs <- readRDS(paste0("scenario2/norm_data/",norm_method,"/norm_mean",batch_mean,"_var",batch_var,"_",phenotype,".rds"))
    sim_meta <- readRDS("scenario2/sim_data/sim_meta.rds")
    pred_df <- data.frame(batch_mean=numeric(0),batch_var=numeric(0),norm_method=character(0),pred_method=character(0),
                          phenotype=character(0),rmse=numeric(0),mape=numeric(0))
    for(i in 1:100){
      pred_res <- sim.pred.func(trn=norm_tabs[[i]][[1]],tst=norm_tabs[[i]][[2]],
                                trn_meta=sim_meta[[i]][[1]],tst_meta=sim_meta[[i]][[2]],
                                phenotype=phenotype,pred_method=pred_method)
      pred_df[i,] <- c(batch_mean=batch_mean,batch_var=batch_var,norm_method=norm_method,pred_method=pred_method,
                       phenotype=phenotype,rmse=round(pred_res["rmse"],3),mape=round(pred_res["mape"],3))
    }
    saveRDS(pred_df,paste0("scenario2/pred_results/",norm_method,"/pred_",pred_method,"_mean",batch_mean,"_var",batch_var,"_",phenotype,".rds"))
    print(paste0("batch_mean=",batch_mean,",batch_var=",batch_var,",phenotype=",phenotype,",norm_method=",norm_method,",pred_method=",pred_method))
  }
}







