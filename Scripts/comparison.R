library(dplyr)
library(ggplot2);library(ggpubr);library(pheatmap); library(gplots)
library("pheatmap")
library("RColorBrewer")
library(tidyr)
library(data.table)
library(monet)
library(ComplexHeatmap)

y_pr_path <- here::here("data/processed/ge_full_Sp1WT.txt") #RNA-seq Wild Type file path
a_pr_path <- here::here("data/processed/ATAC_WT-time.txt") #ATAC-seq Wild Type file path

## PROCESSED DATA 
y_pr<- fread(y_pr_path)
y_pr<- monet:::geneDtMelt(y_pr) %>% monet:::geneStanderdise() %>% dcast(gene ~ variable)
genes<-y_pr[,1]
y_pr<- y_pr[,-1]
y_pr$total_exp <- rowSums(y_pr) #add a total expression column
y_pr<- cbind(genes, y_pr)
y_pr<- data.frame(y_pr)

a_pr <- read.table(a_pr_path, header=TRUE, sep="\t")
names(a_pr) <- c("Gene","HB","HE1","HE2","HP")

## SIMULATED DATA 
y_sim <- data.table(t(monet_sim$y))
y_sim$gene <- seq(1:nrow(y_sim))
a_sim <- data.table(t(monet_sim$a))
names(a_sim) <- c("HB","HE1","HE2","HP")

y0_sim<-  data.table(monet_sim$y0)
y_sim<-cbind(y0_sim, y_sim) # add time point 0 to gene expression data
names(y_sim)[1] <- "ES_Sp1WT"

# standardize all the simulated data
y_sim <- cbind(y_sim[, "gene"], asinh(y_sim[, -"gene"])) %>% monet:::geneDtMelt() %>% monet:::geneStanderdise() %>% dcast(gene ~ variable)
y_sim<- y_sim[,-"gene"]
y_sim$total_exp <- rowSums(y_sim) #add a total expression column

####### RNASEQ comparison #######

# Densities of the total expression distributions
pdf("results/plots/Comparison_plots/2000g/density_totalexp_sim_30TF_8000it_2000wup.pdf") 
ggplot(y_sim, aes(x=total_exp)) + 
  geom_density(color="black", fill="#D2D2D2")+ 
  geom_vline(aes(xintercept=mean(total_exp)), color="black", linetype="dashed", size=1)+
  labs(title="Distribution of the total expression of genes of the simulation",x="Total expression", y = "Count")+
  theme_classic()
dev.off()

pdf("results/plots/Comparison_plots/density_totalexp_pr.pdf") 
ggplot(y_pr, aes(x=total_exp)) + 
  geom_density(color="black", fill="#ABCCE2")+ 
  geom_vline(aes(xintercept=mean(total_exp)), color="black", linetype="dashed", size=1)+
  labs(title="Distribution of the total expression of genes of the processed data",x="Total expression", y = "Count")+
  theme_classic()
dev.off()

# Boxplots
  # put data into long format
names(y_sim)<- c("ES","HB","HE1","HE2","HP", "total_exp")
y_sim_long <- y_sim %>% pivot_longer(cols= c("ES","HB","HE1","HE2","HP"),
                                     names_to='Time_point',
                                     values_to='Expression')

pdf("results/plots/Comparison_plots/2000g/boxplots_sim_30TF_8000it_2000wup.pdf")
ggplot(y_sim_long, aes(x=Time_point, y=Expression, fill=Time_point)) + geom_boxplot()+ theme_classic()+ labs(title="Plot of expression per time point of simulated data",fill="Time point", y = "Expression")+ scale_fill_brewer(palette="Dark2",labels=c('ES', 'HB', 'HE1', 'HE2','HP'))
dev.off()

names(y_pr)<- c("gene","ES","HB","HE1","HE2","HP", "total_exp")
y_pr_long <- y_pr %>% pivot_longer(cols=c("ES","HB","HE1","HE2","HP"),
                                   names_to='Time_point',
                                   values_to='Expression')

pdf("results/plots/Comparison_plots/boxplots_pr.pdf")
ggplot(y_pr_long, aes(x=Time_point, y=Expression, fill=Time_point)) + geom_boxplot()+ theme_classic()+ labs(title="Plot of expression per time point of processed real data",fill="Time point", y = "Expression")+scale_x_discrete(limits=c("ES","HB","HE1","HE2","HP"))+ scale_fill_brewer(palette="Dark2",labels=c('ES', 'HB', 'HE1', 'HE2','HP'))
dev.off()

# The two boxplots in the same figure
y_sim_long$Data <- rep("Simulated",nrow(y_sim_long))
y_pr_long$Data <- rep("Observed",nrow(y_pr_long))
data_long<- rbind(y_sim_long, y_pr_long[,2:5])
data_long$Time_point[data_long$Time_point=="ES_Sp1WT"] <- "ES"
data_long$Time_point[data_long$Time_point=="Flk_Sp1WT"] <- "HB"
data_long$Time_point[data_long$Time_point=="HE1_Sp1WT"] <- "HE1"
data_long$Time_point[data_long$Time_point=="HE2_Sp1WT"] <- "HE2"
data_long$Time_point[data_long$Time_point=="HP_Sp1WT"] <- "HP"

pdf("results/plots/Comparison_plots/2000g/boxplots_30TF_8000it_2000wup.pdf") 
ggplot(data_long, aes(x=Time_point, y=Expression, fill=interaction(Time_point,Data))) + geom_boxplot()+ theme_classic()+ labs(title="Expression distribution per time point of simulated and observed data",fill="Data", y = "Expression") + scale_fill_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#91dec7", "#f6ae76", "#bfbcec","#fe97cb" ,"#bde590" ))
dev.off()

## Densities per time point
y_sim<- data.frame(y_sim)
pdf("results/plots/Comparison_plots/2000g/densities_timepoints_30TF_8000it_2000wup.pdf") 
list_plots=list()
legend=c("ES","HB","HE1","HE2","HP")
colors=c("#AFE183", "#6FDBA5", "#6FD7DB", "#6FB0DB", "#6485C6")
for (n in 1:5){
  x<- data.frame(timepoint= y_sim[, n], Data=rep("Simulated", nrow(y_sim)))
  y<-data.frame(timepoint=y_pr[,n+1], Data=rep("Observed",nrow(y_pr)))
  data <-rbind(x,y)
  p <- data %>%
    ggplot( aes(x=timepoint, fill=Data)) +
    geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c(colors[n], "#3a3a3a")) +
    theme_classic() + xlab(legend[n])
  
  
  list_plots[[paste0("run ", n)]] =p
}
gridExtra::grid.arrange(grobs = list_plots, ncol=2,top="Density distributions of gene expression at different time points", common.legend = TRUE, legend="bottom")
dev.off()

# Plot of 10 random genes trajectories
# Plot of 10 random genes trajectories
# Take a random susbet of 10 genes of the simulated data set
sim_10 <- subset(y_sim, rownames(y_sim) %in% sample(1:nrow(y_sim), 10, replace=F)) 
sim_10<-sim_10[,-1]
names(sim_10)<- c("1","2","3","4","5")
sim_10$gene <- rownames(sim_10)
# put data in long format
sim_10_long <- sim_10 %>% pivot_longer(cols=c("1","2","3","4","5"),names_to='Time_point',values_to='Expression')
sim_10_long$Time_point = as.numeric(sim_10_long$Time_point)

# Take a random susbet of 10 genes of the processed data set
pr_10  <- subset(y_pr, rownames(y_pr) %in% sample(1:nrow(y_pr), 10, replace=F))[,-7]
names(pr_10)<- c("gene","1","2","3","4","5")
# put data in long format
pr_10_long <- pr_10 %>% pivot_longer(cols=c("1","2","3","4","5"),names_to='Time_point',values_to='Expression')
pr_10_long$Time_point = as.numeric(pr_10_long$Time_point)

# Trajectories in grid form
pdf("results/plots/Comparison_plots/2000g/trajectories_grid_sim_30TF_8000it_2000wup.pdf",width = 8) 
ggplot(sim_10_long, aes(x=Time_point, y= Expression, color=as.factor(gene), group=gene))+ geom_line()+theme_classic() + theme(legend.position = "none") + labs(title = "Trajectories of 10 genes of simulated data ", y = "Expression", x = "Time point")+ scale_x_discrete(limits=c("ES","HB","HE1","HE2","HP")) + facet_wrap(~gene, ncol=2)
dev.off()

pdf("results/plots/Comparison_plots/trajectories_grid_pr.pdf",width = 6) 
ggplot(pr_10_long, aes(x=Time_point, y= Expression, color=as.factor(gene), group=gene))+ geom_line(linewidth=1.2)+theme_classic() + theme(legend.position = "none") + labs(title = "Trajectories of 10 genes of observed data ", y = "Expression", x = "Time point")+ scale_x_discrete(limits=c("ES","HB","HE1","HE2","HP"))+ facet_wrap(~gene, ncol=2)
dev.off()

p3 <- ggplot(sim_10_long, aes(x=Time_point, y= Expression, color=as.factor(gene), group=gene))+ geom_line(linewidth=1)+theme_classic() + theme(legend.position = "none") + labs(title = "Trajectories of 10 genes of simulated data ", y = "Expression", x = "Time point")+ scale_x_discrete(limits=c("ES","HB","HE1","HE2","HP")) + facet_wrap(~gene, ncol=2)
p4 <- ggplot(pr_10_long, aes(x=Time_point, y= Expression, color=as.factor(gene), group=gene))+ geom_line(linewidth=1)+theme_classic() + theme(legend.position = "none") + labs(title = "Trajectories of 10 genes of observed data ", y = "Expression", x = "Time point")+ scale_x_discrete(limits=c("ES","HB","HE1","HE2","HP"))+ facet_wrap(~gene, ncol=2)
arrange <- ggarrange(p4, p3, ncol = 1, nrow = 2)
ggsave(path="results/plots/Comparison_plots/2000g", filename="trajectories_grid_pr_sim_30TF_8000wup_2000it.png", arrange,width = 5.5)

####### ATACSEQ comparison #######
summary(a_sim)
summary(a_pr)

a_sim_matrix <- as.matrix(a_sim)
a_pr_matrix <- as.matrix(a_pr[,2:5])
rownames(a_pr_matrix)<- a_pr[,1]

pdf("results/plots/Comparison_plots/2000g/Heatmap_sim_30TF_8000it_2000wup.pdf") 
Heatmap(a_sim_matrix, name="Chromatin state", col=c("#FF7465", "#84E875"),
column_title = "Chromatin state of the simulated data", column_title_gp = gpar(fontsize = 15, fontface = "bold"))
dev.off()

rownames(a_pr_matrix) = NULL
pdf("results/plots/Comparison_plots/2000g/Heatmap_pr.pdf") 
Heatmap(a_pr_matrix, name="Chromatin state", col=c("#FF7465", "#84E875"),
column_title = "Chromatin state of the observed data", column_title_gp = gpar(fontsize = 15, fontface = "bold"))
dev.off()

colMeans(a_sim_matrix)
colMeans(a_pr_matrix)
