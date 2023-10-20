########################################################################################################################
### figures in the draft
### 2023/10/20
########################################################################################################################

setwd("/Users/beibeiwang/Desktop/Normalization_quantitative/")

### packages
all(sapply(c("vegan","ape","reshape2","ggplot2","ggtree","viridis","RColorBrewer","ggsci","aplot","ggpubr","ggrepel",
             "grid","ggplotify","patchwork","dplyr"),require, character.only=TRUE))


#======================================================================================================================#
### part 1: basic analysis using datasets in curatedMetagenomicData ###
#======================================================================================================================#
### data 
meta <- readRDS("data/meta.rds")
count <- readRDS("data/count.rds")
norm_methods <- readRDS("data/norm_methods.rds")


#======================================================================================================================#
### supplementary figure s1: BMI
boxplot_bmi <- ggplot(meta, aes(x=study_name,y=BMI,fill=study_name))+
  geom_boxplot(alpha=0.7)+
  geom_hline(yintercept=mean(meta$BMI),linetype="dashed",color="grey")+
  annotate(geom="curve",x=28.5,y=12,xend=31.2,yend=mean(meta$BMI),curvature=.3,arrow=arrow(length=unit(2,"mm")))+
  annotate(geom="text",x=26,y=11,label=paste0("average BMI = ",round(mean(meta$BMI),2)))+
  stat_compare_means(method="wilcox.test",ref.group=".all.",label="p.signif")+
  stat_summary(fun="mean", geom="point", shape=20, size=2.5, color="red", fill="red",alpha=0.7)+
  scale_fill_igv("default")+
  labs(x="",y="BMI",title="")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1))+
  theme(legend.position="none")
boxplot_bmi
ggsave("figures/S1_fig.tiff",width=10,height=6,dpi=300)


#======================================================================================================================#
### supplementary figure s2: shannon index
count_norm <- as.data.frame(apply(count,2,function(x) x/sum(x)))
shannon_div <- diversity(t(count_norm),index='shannon')
shannon_df <- data.frame(sample_id=names(shannon_div),shannon=shannon_div)
shannon_df <- merge(shannon_df,meta,by="sample_id")
boxplot_shannon <- ggplot(shannon_df, aes(x=study_name,y=shannon,fill=study_name))+
  geom_boxplot(alpha=0.7)+
  geom_hline(yintercept=mean(shannon_df$shannon),linetype="dashed",color="grey")+
  annotate(geom="curve",x=29,y=0.5,xend=31.2,yend=mean(shannon_df$shannon),curvature=0.2,arrow=arrow(length=unit(2,"mm")))+
  annotate(geom="text",x=26,y=0.3,label=paste0("average Shannon Index= ",round(mean(shannon_df$shannon),2)))+
  stat_compare_means(method="wilcox.test",ref.group=".all.",label="p.signif")+
  stat_summary(fun="mean", geom="point", shape=20, size=2.5, color="red", fill="red",alpha=0.7)+
  scale_fill_igv("default")+
  labs(x="",y="Shannon Index",title="")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1))+
  theme(legend.position="none")
boxplot_shannon
ggsave("figures/S2_fig.tiff",width=10,height=6,dpi=300)


#======================================================================================================================#
### figure 1: basic analysis 
### pcoa plot
count_norm <- as.data.frame(apply(count,2,function(x) x/sum(x)))
bray_div <- vegdist(t(count_norm),method="bray")
pcoa_lst <- pcoa(bray_div)
pcoa_df <- data.frame(sample_id=rownames(pcoa_lst$vectors),axis1=pcoa_lst$vectors[,1],axis2=pcoa_lst$vectors[,2])
pcoa_df <- merge(pcoa_df,meta,by="sample_id")
pcoa_summ_df <- pcoa_df %>% group_by(study_name) %>% summarise(mean_axis1=mean(axis1),mean_axis2=mean(axis2),sample_size=length(sample_id))
scatter_pcoa <- ggplot(pcoa_summ_df,aes(x=mean_axis1,y=mean_axis2,color=study_name))+
  geom_point(aes(size=sample_size),alpha=0.7)+
  geom_text_repel(aes(label=study_name))+
  scale_color_manual(values=pal_igv("default")(length(unique(pcoa_df$study_name))))+
  labs(x=paste0("PCoA1(",round(pcoa_lst$values[,"Relative_eig"][1],3)*100,"%)"), y=paste0("PCoA2(",round(pcoa_lst$values[,"Relative_eig"][2],3)*100,"%)"),title="PCoA")+
  theme_bw()+theme(legend.position="none")
scatter_pcoa


### permanova analysis
# To avoid problems related to variable ordering, total variance explained by each variable was evaluated independently 
# of other variables, and thus should be regarded as total variance explainable by that variable (reference: HMP2 paper)
# considered covariates: study, country, dna extraction kit, age, gender
# the samples with NA values were removed from analysis
permanova_df <- data.frame(variable=character(0),r2=numeric(0),pvalue=numeric(0))
covariates <- c("study_name","country","DNA_extraction_kit","sequencing_platform","age","gender","BMI")
permanova.func <- function(distance,meta_df,variable){
  meta_df <- meta_df[!is.na(meta_df[,variable]),]
  distance <- as.dist(as.matrix(distance)[meta_df$sample_id,meta_df$sample_id])
  formulas <- as.formula(paste0("distance~",variable))
  permanova <- adonis2(formulas,data=meta_df,permutations=1000,parallel=10)
  c(variable,round(permanova[variable,"R2"],3),permanova[variable,"Pr(>F)"])
}
for(i in 1:length(covariates)){
  permanova_df[i,] <- permanova.func(distance=bray_div,meta_df=meta,variable=covariates[i])
  print(covariates[i])
}
saveRDS(permanova_df,"results/permanova_res.rds")
# as the p value for considered covariates all equal to 0.001, so annotated the p value *** at the top of the bar rather than a separated bar plot
permanova_df$r2 <- as.numeric(permanova_df$r2)
permanova_df$labels <- rep("***",nrow(permanova_df))
permanova_df$variable <- factor(permanova_df$variable,levels=permanova_df$variable[order(permanova_df$r2,decreasing=T)])
bar_permanova <- ggplot(permanova_df, aes(x=variable,y=r2,fill=variable))+
  geom_bar(stat="identity",width=0.8)+
  geom_text(aes(label=labels),color="black",size=4,hjust=0,vjust=1)+
  scale_x_discrete(breaks=c("study_name","country","DNA_extraction_kit","sequencing_platform","age","gender","BMI"),
                   labels=c("Datasets","Country","DNA-Exk","Seq-Plat","Age","Gender","BMI"))+
  scale_fill_manual(values=rep("steelblue",nrow(permanova_df)))+
  labs(x="Variable",y="PERMANOVA R2",title="PERMANOVA")+
  theme_bw()+theme(legend.position="none")+coord_flip()
bar_permanova


### heatmapes for rmse using TSS normalized abundances
res_summ_order_df <- readRDS("results/curatedMetagenomicData_res_summ_order.rds")
tss_pred_df <- res_summ_order_df[res_summ_order_df$norm_method=="TSS",]
tss_pred_df$median_rmse <- round(tss_pred_df$median_rmse,2)
heatmap_tss_pred_rmse <- ggplot(tss_pred_df, aes(x=study2,y=study1,fill=median_rmse))+ 
  geom_tile()+
  scale_fill_distiller(palette="RdBu",direction=-1,na.value="transparent")+
  geom_text(aes(study2,study1,label=median_rmse),color="black",size=4)+
  labs(x="Dataset 1", y="Dataset 2", fill="Median RMSE")+
  labs(title="Median RMSE")+
  theme_bw()+theme(legend.position="bottom")+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(panel.grid.minor=element_blank())+
  theme(panel.border=element_blank())
heatmap_tss_pred_rmse


### combine the plots
ggarrange(ggarrange(scatter_pcoa,bar_permanova,ncol=2,nrow=1,widths=c(3,1),labels=c("A","B")),
          heatmap_tss_pred_rmse,ncol=1,nrow=2,heights=c(1,1.5),labels=c("","C"))
ggsave("figures/Fig1.tiff",width=13,height=13,dpi=300)


#======================================================================================================================#
### part 2: simulation analysis ###
#======================================================================================================================#
### data
scenario1_summ_df <- readRDS("data/scenario1_res_summ.rds")
scenario2_summ_df <- readRDS("data/scenario2_res_summ.rds")
scenario3_summ_df <- readRDS("data/scenario3_res_summ.rds")
norm_methods <- readRDS("data/norm_methods.rds")


#======================================================================================================================#
### figure 2: heatmaps for rmse in scenario1
scenario1_summ_df <- merge(scenario1_summ_df,norm_methods,by="norm_method")
scenario1_summ_df$alpha <- factor(scenario1_summ_df$alpha,levels=c(0,0.2,0.4,0.6,0.8,1))
scenario1_summ_df$median_rmse <- round(scenario1_summ_df$median_rmse,2)
scenario1.heatmap.rmse.func <- function(res_summ_df,title){
  heatmaps <-  ggplot(res_summ_df, aes(x=alpha, y=annotated_method, fill=median_rmse))+ 
    geom_tile() +
    geom_text(aes(x=alpha, y=annotated_method, label=median_rmse), color="black", size=4)+
    labs(x="population effect",y="method",fill="Median RMSE",title=title)+
    scale_fill_distiller(palette="Reds",direction=-1)+
    theme_bw()+
    theme(legend.position="bottom")+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
scenario1_rmse_heatmaps <- list()
scenario1_rmse_heatmaps[[1]] <- scenario1.heatmap.rmse.func(res_summ_df=scenario1_summ_df[scenario1_summ_df$phenotype=="linear",],title="Linear")
scenario1_rmse_heatmaps[[2]] <- scenario1.heatmap.rmse.func(res_summ_df=scenario1_summ_df[scenario1_summ_df$phenotype=="quadratic",],title="Quadratic")
scenario1_rmse_heatmaps[[3]] <- scenario1.heatmap.rmse.func(res_summ_df=scenario1_summ_df[scenario1_summ_df$phenotype=="inverse",],title="Inverse")
scenario1_rmse_heatmaps[[4]] <- scenario1.heatmap.rmse.func(res_summ_df=scenario1_summ_df[scenario1_summ_df$phenotype=="logistic",],"Logistic")
# combine the plots
ggarrange(plotlist=scenario1_rmse_heatmaps,labels=LETTERS[1:4],nrow=2,ncol=2)
ggsave("figures/Fig2.tiff",width=10,height=10,dpi=300)


#======================================================================================================================#
### figure 3: heatmaps for rmse in scenario2
scenario2_summ_df <- merge(scenario2_summ_df,norm_methods,by="norm_method")
scenario2_summ_df$group <- paste0("m",scenario2_summ_df$batch_mean,"v",scenario2_summ_df$batch_var)
scenario2_summ_df$group <- factor(scenario2_summ_df$group,levels=c("m0v1","m0v2","m0v4","m500v1","m1000v1"))
scenario2_summ_df$median_rmse <- round(scenario2_summ_df$median_rmse,2)
scenario2.heatmap.rmse.func <- function(res_summ_df,title){
  heatmaps <-  ggplot(res_summ_df, aes(x=group, y=annotated_method, fill=median_rmse))+ 
    geom_tile() +
    geom_text(aes(x=group, y=annotated_method, label=median_rmse), color="black", size=4)+
    labs(x="factor",y="method",fill="Median RMSE",title=title)+
    scale_fill_distiller(palette="Reds",direction=-1)+
    #scale_x_discrete(position="top")+
    theme_bw()+theme(legend.position="bottom")+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
scenario2_rmse_heatmaps <- list()
scenario2_rmse_heatmaps[[1]] <- scenario2.heatmap.rmse.func(res_summ_df=scenario2_summ_df[scenario2_summ_df$phenotype=="linear",],title="Linear")
scenario2_rmse_heatmaps[[2]] <- scenario2.heatmap.rmse.func(res_summ_df=scenario2_summ_df[scenario2_summ_df$phenotype=="quadratic",],title="Quadratic")
scenario2_rmse_heatmaps[[3]] <- scenario2.heatmap.rmse.func(res_summ_df=scenario2_summ_df[scenario2_summ_df$phenotype=="inverse",],title="Inverse")
scenario2_rmse_heatmaps[[4]] <- scenario2.heatmap.rmse.func(res_summ_df=scenario2_summ_df[scenario2_summ_df$phenotype=="logistic",],"Logistic")
# combine the plots
ggarrange(plotlist=scenario2_rmse_heatmaps,labels=LETTERS[1:4],nrow=2,ncol=2)
ggsave("figures/Fig3.tiff",width=10,height=10,dpi=300)


#======================================================================================================================#
### figure 4: heatmaps for rmse in scenario3
scenario3_summ_df <- merge(scenario3_summ_df,norm_methods,by="norm_method")
scenario3_summ_df$overlap <- factor(scenario3_summ_df$overlap,levels=c(2,4,6,8,10))
scenario3_summ_df$median_rmse <- round(scenario3_summ_df$median_rmse,2)
scenario3.heatmap.rmse.func <- function(res_summ_df,title){
  heatmaps <-  ggplot(res_summ_df, aes(x=overlap, y=annotated_method, fill=median_rmse))+ 
    geom_tile() +
    geom_text(aes(x=overlap, y=annotated_method, label=median_rmse), color="black", size=4)+
    labs(x="overlap",y="method",fill="Median RMSE",title=title)+
    scale_fill_distiller(palette="Reds",direction=-1)+
    #scale_x_discrete(position="top")+
    theme_bw()+theme(legend.position="bottom")+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
scenario3_rmse_heatmaps <- list()
scenario3_rmse_heatmaps[[1]] <- scenario3.heatmap.rmse.func(res_summ_df=scenario3_summ_df[scenario3_summ_df$phenotype=="linear",],title="Linear")
scenario3_rmse_heatmaps[[2]] <- scenario3.heatmap.rmse.func(res_summ_df=scenario3_summ_df[scenario3_summ_df$phenotype=="quadratic",],title="Quadratic")
scenario3_rmse_heatmaps[[3]] <- scenario3.heatmap.rmse.func(res_summ_df=scenario3_summ_df[scenario3_summ_df$phenotype=="inverse",],title="Inverse")
scenario3_rmse_heatmaps[[4]] <- scenario3.heatmap.rmse.func(res_summ_df=scenario3_summ_df[scenario3_summ_df$phenotype=="logistic",],"Logistic")+
  scale_fill_distiller(palette="Reds",limits=c(0.975,1.03),breaks=c(0.98,1.00,1.02),labels=c("0.98", "1.00", "1.02"))
# combine the plots
ggarrange(plotlist=scenario3_rmse_heatmaps,labels=LETTERS[1:4],nrow=2,ncol=2)
ggsave("figures/Fig4.tiff",width=10,height=10,dpi=300)


#======================================================================================================================#
### part 3: real data prediction ###
#======================================================================================================================#
### data
norm_methods <- readRDS("data/norm_methods.rds")
res_summ_order_df <- readRDS("data/curatedMetagenomicData_res_summ_order.rds")
res_summ_order_df <- merge(res_summ_order_df,norm_methods,by="norm_method")


### supplementary figure s3: boxplot for median RMSE of certain test data
studies <- unique(res_summ_order_df$study1)
studies <- studies[order(studies)]
box_pred_rmse <- list()
for(study in studies){
  box_pred_rmse[[study]] <- ggplot(res_summ_order_df[res_summ_order_df$study2==study,], aes(x=annotated_method,y=median_rmse,fill=class))+
    geom_boxplot(alpha=0.7)+
    stat_compare_means(method="wilcox.test",ref.group=".all.",label="p.signif",vjust=0.5)+
    stat_summary(fun="mean", geom="point", shape=20, size=2.5, color="red", fill="red",alpha=0.7)+
    scale_fill_jco()+
    labs(x="",y="Median RMSE",fill="Group of normalization methods",title=paste0("test: ",study))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
    theme(legend.position="bottom")
  #if(study%in%studies[1:(length(studies)-4)]) box_pred_rmse[[study]] <- box_pred_rmse[[study]]+labs(x="")+theme(axis.text.x = element_blank())
}
ggarrange(plotlist=box_pred_rmse,nrow=8,ncol=4,common.legend=T,labels=paste0(rep(LETTERS[1:8],each=4),rep(1:4,8))[1:length(studies)])
ggsave("figures/S3_fig.tiff",width=18,height=20,dpi=300)


### figure 5: boxplot for the ranks of normalization method according to median rmse
box_ranks_rmse <- ggplot(res_summ_order_df, aes(x=annotated_method,y=rank_median_rmse,fill=class))+
  #geom_boxplot(alpha=0.7)+
  geom_violin(width=0.9,trim=FALSE,alpha=0.7)+
  geom_boxplot(width=0.2,fill="white",position=position_dodge(width=0.9))+
  scale_fill_jco()+
  labs(x="",y="Rank",fill="Group of normalization methods",title="Median RMSE Ranks") +
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(legend.position="bottom")
box_ranks_rmse+labs(title="")
ggsave("figures/Fig5.tiff",width=10,height=4,dpi=300)





