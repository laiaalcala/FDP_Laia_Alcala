library(tidyverse); library(dplyr)
library(monet); library(rstan)
library(stringr)
library(ggplot2); library(ggpubr); library(viridis);library(gridExtra)
library(data.table)
library(matrixcalc) # for the hardmard product
library(ComplexHeatmap); library(circlize)

no_g <- 1000 # Number of genes we want to work with
no_tf <- 20 # Number of transcription factors
de_b <- 2 
alpha <- 1
iter <- c(12000) # Number of iterations
wup <- c(1500) # warm-up
test <- 1:10

SD <- 2

#### FOR WILD TYPE #####
gn_type <- "WT"

data_path <- here::here("data/processed")
#gene_path <-here::here("data/processed/ge_full_Sp1WT_top.txt")
gene_path <- "~/cazierj-msc-bioinf-dl/lca259/20TF_0.5_NINK/ge_full_Sp1WT_top_nink_0.5.txt"
#atac_path <-here::here("data/processed/ATAC_WT-time_top.txt")
atac_path <-"~/cazierj-msc-bioinf-dl/lca259/20TF_0.5_NINK/ATAC_Sp1WT_time_top_nink_0.5.txt"

# Function to initialize the monetData class. 
# Inputs: the RNA-seq and ATAC-seq data paths, the column of gene names in the RNA-seq data, the transformation function to be applied and the column of gene names in the ATAC-seq data 
dat_monet <-monet:::monetData_init(gene_path, atac_path, "gene", "asinh", "Gene")

# Function to verify that RNA-seq and ATAC-seq data contain the same genes. It creates a subset with the matching genes.
# Inputs: the monetData class and the number of genes we want to obtain in the subset.
monet_data <-monet::monet_data_filter(dat_monet, no_genes = no_g)

# Function to prepare input monet class for Rstan model
# Inputs: the monet data class and a named list containing the number of TF and b.The output is a list for monet stan model.
monet_l <- monet::monet_input_prep(monet_data,list(no_tf = no_tf, de_b = de_b))

# FIXME: Remove optim results and add functions to extract sampling information
# Fit the monet stan model to the preprocessed data
optim_out <- stan_monet.fit(monet_l, iter = iter,warmup = wup, cores = 4)

save(monet_data, optim_out,
     file = stringr::str_c("results/tables/obs_laia_1000/single-optim_",
                           no_g, "_", no_tf, "-", de_b, "-",
                           alpha,"-", gn_type,"-", SD, ".Rdata"))
load(stringr::str_c(here::here("results/tables/obs_laia_1000/single-optim_"),
                    no_g, "_", no_tf, "-", de_b, "-", alpha,"-", gn_type,"-", SD, ".Rdata"))

# MCMC
pdf("results/plots/Obs_data_plots/MCMC_trace_logpost_WT_0.5.pdf", width=25) 
traceplot(optim_out, "lp__")+
  labs(title = "MCMC log posterior traceplot",subtitle ="1000 genes - 20TFs - 12000 iterations - 1500 w.up",x = "Iterations", y= "Log posterior")+
  theme(axis.text = element_text(size = 12),
    axis.title.y = element_text(size = rel(3)),
    axis.title.x = element_text(size = rel(3)),
    plot.title = element_text(size = rel(3.2)),
    plot.subtitle=element_text(size=rel(3)),
    legend.text=element_text(size=rel(2.8)),
    legend.title=element_text(size=rel(2.5)),
    legend.key.size = unit(2, 'cm'))
dev.off()

# Extract posterior means
w_all_WT <- get_posterior_mean(optim_out, "w")
b_all_WT <- get_posterior_mean(optim_out, "b")
x_all_WT <- get_posterior_mean(optim_out, "x")

# Extract posterior across all chains (last column)
w_est_WT <- w_all_WT[, ncol(w_all_WT)]
b_est_WT <- b_all_WT[, ncol(b_all_WT)]
x_est_WT <- x_all_WT[, ncol(x_all_WT)]

bgk_WT <- as.data.frame(b_est_WT) 
xk_WT <- as.data.frame(x_est_WT)

xk_matrix_WT <- matrix(NA,ncol=20, nrow=4 )
c<-0
for (i in seq(nrow(xk_WT)/20)){
  c <- c+1
  xk_matrix_WT[i,] <- xk_WT[c:(c+19),]
  c <- c+19
}

bgk_matrix_WT <- matrix(NA,ncol=20, nrow=1000 )
c<-0
for (i in seq(nrow(bgk_WT)/20)){
  c <- c+1
  bgk_matrix_WT[i,] <- bgk_WT[c:(c+19),]
  c <- c+19
}

xk_matrix_WT <- t(xk_matrix_WT) # Dimension: TF(10) x Time points (4)
bgk_matrix_WT <- t(bgk_matrix_WT); colnames(bgk_matrix_WT)= monet_data@gene_filt_list  # Dimension: TF (10) x Genes (1000) 
wk_matrix_WT <- as.matrix(diag(w_est_WT)) # Dimension: TF (10) x TF (10)

genes_WT <- colnames(bgk_matrix_WT)

# Heatmaps of B and X
names(bgk_matrix_WT) <- NULL
rownames(bgk_matrix_WT) <- seq(1:20)
col_fun = colorRamp2(c(min(bgk_matrix_WT),min(bgk_matrix_WT) / 10, min(bgk_matrix_WT) / 100, 0,max(bgk_matrix_WT) / 10, max(bgk_matrix_WT) / 100, max(bgk_matrix_WT)),
                     c( "darkblue", "cadetblue3", "lightcyan", "white", "lightcoral", "orange", "red"))
pdf("results/plots/Obs_data_plots/bgk_PE_WT_20TF_1000g_2.pdf", width=12) 
h_bgk <- Heatmap(bgk_matrix_WT, name="Gene expression",col=col_fun, cluster_rows = FALSE,row_names_side = "left", show_column_names=FALSE, show_row_dend = FALSE,column_dend_height = unit(3, "cm"), 
                 row_dend_width = unit(3, "cm"),heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 14)))
draw(h_bgk, heatmap_legend_side = "bottom")
dev.off()

x_df_WT <- as.data.frame(xk_matrix_WT)
names(x_df_WT) <- c("1","2","3","4")
rownames(x_df_WT) <- seq(1:20)

pdf("results/plots/Obs_data_plots/xk_PE_WT_20TF_1000g_2.pdf", width=3, height=6 )
h_xk <- Heatmap(as.matrix(x_df_WT), name="WT TF interaction",cluster_rows = FALSE,cluster_columns  = FALSE, row_title= "TF" ,show_row_dend = FALSE, show_column_dend = FALSE, row_names_side = "left",col=c("purple", "white"),heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 12)))
draw(h_xk, heatmap_legend_side = "bottom")
dev.off()

# Weights DP
w_est_WT_df <- as.data.frame(w_est_WT)
w_est_WT_df$TF <- seq(1:20)

pdf("results/plots/Obs_data_plots/wk_PE_WT_20TF_1000g_2.pdf", width=2, height=4) 
ggplot ( w_est_WT_df, aes(as.factor(TF),w_est_WT)) +geom_bar(stat = "identity", fill="lightsalmon") +ylim(c(0, 0.5)) +labs ( y= "DP weight", x="TF")+ theme_classic()+coord_flip()
dev.off()

# yg computation
rnaseq_WT <- monet_data@gene_exp$E14_flk_fpkm #y0
atacseq_WT <- as.data.frame(monet_data@atac_seq)
atacseq_WT <- atacseq_WT[atacseq_WT$gene %in% genes_WT,]
atacseq_WT <- as.matrix(atacseq_WT[,-1])
class(atacseq_WT) <- "numeric"
yg_WT <- as.data.frame(monet_sim_yg(rnaseq_WT, atacseq_WT, bgk_matrix_WT, xk_matrix_WT,wk_matrix_WT))


#### FOR SP1 KNOCKOUT #####
gn_type <- "Sp1hyp"

data_path <- here::here("data/processed")
#chip_path <-here::here("data/processed/ES_ChIP_binary_table.txt")
#gene_path <-here::here("data/processed/ge_full_Sp1hyp_top.txt")
gene_path <- "~/cazierj-msc-bioinf-dl/lca259/20TF_0.5_NINK/ge_full_Sp1hyp_top_nink_0.5.txt"
#atac_path <-here::here("data/processed/ATAC_Sp1hyp-time_top.txt")
atac_path <- "~/cazierj-msc-bioinf-dl/lca259/20TF_0.5_NINK/ATAC_Sp1hyp_time_top_nink_0.5.txt"

dat_monet <-monet:::monetData_init(gene_path, atac_path, "gene", "asinh", "Gene")
monet_data <-monet::monet_data_filter(dat_monet, no_genes = no_g)
monet_l <- monet::monet_input_prep(monet_data,list(no_tf = no_tf, de_b = de_b))

optim_out_Sp1hyp <- stan_monet.fit(monet_l, iter = iter,warmup = wup, cores = 4)

save(monet_data, optim_out_Sp1hyp,
     file = stringr::str_c(here::here("results/tables/obs_laia_1000/single-optim_"),
                           no_g, "_", no_tf, "-", de_b, "-",
                           alpha,"-", gn_type,"-", SD, ".Rdata"))

load(stringr::str_c(here::here("results/tables/obs_laia_1000/single-optim_"),
                    no_g, "_", no_tf, "-", de_b, "-", alpha,"-", gn_type,"-", SD, ".Rdata"))

load("~/cazierj-msc-bioinf-dl/laia-data/single-optim_1000_10-2-1-Sp1hyp.Rdata")

# MCMC
pdf("results/plots/Obs_data_plots/MCMC_trace_logpost_Sp1hyp_0.5.pdf", width=25) 
traceplot(optim_out, "lp__")+
  labs(title = "MCMC log posterior traceplot",subtitle ="1000 genes - 20TFs - 12000 iterations - 1500 w.up",x = "Iterations", y= "Log posterior")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        plot.title = element_text(size = rel(3.2)),
        plot.subtitle=element_text(size=rel(3)),
        legend.text=element_text(size=rel(2.8)),
        legend.title=element_text(size=rel(2.5)),
        legend.key.size = unit(2, 'cm'))
dev.off()

# Extract posterior means
w_all_SP1 <- get_posterior_mean(optim_out_Sp1hyp, "w")
b_all_SP1 <- get_posterior_mean(optim_out_Sp1hyp, "b")
x_all_SP1 <- get_posterior_mean(optim_out_Sp1hyp, "x")

# Extract posterior across all chains (last column)
w_est_SP1<- w_all_SP1[, ncol(w_all_SP1)]
b_est_SP1<- b_all_SP1[, ncol(b_all_SP1)]
x_est_SP1 <- x_all_SP1[, ncol(x_all_SP1)]

bgk_SP1 <- as.data.frame(b_est_SP1) 
xk_SP1 <- as.data.frame(x_est_SP1)

xk_matrix_SP1 <- matrix(NA,ncol=20, nrow=4 )
c<-0
for (i in seq(nrow(xk_SP1)/20)){
  c <- c+1
  xk_matrix_SP1[i,] <- xk_SP1[c:(c+19),]
  c <- c+19
}

bgk_matrix_SP1 <- matrix(NA,ncol=20, nrow=1000 )
c<-0
for (i in seq(nrow(bgk_SP1)/20)){
  c <- c+1
  bgk_matrix_SP1[i,] <- bgk_SP1[c:(c+19),]
  c <- c+19
}

xk_matrix_SP1 <- t(xk_matrix_SP1) # Dimension: TF(10) x Time points (4)
bgk_matrix_SP1 <- t(bgk_matrix_SP1)  # Dimension: TF (10) x Genes (1000) 
wk_matrix_SP1 <- as.matrix(diag(w_est_SP1)) # Dimension: TF (10) x TF (10)

genes_SP1 <- monet_data@gene_filt_list

# Heatmaps of B and X
names(bgk_matrix_SP1) <- NULL
rownames(bgk_matrix_SP1) <- seq(1:20)
col_fun = colorRamp2(c(min(bgk_matrix_SP1),min(bgk_matrix_SP1) / 10, min(bgk_matrix_SP1) / 100, 0,max(bgk_matrix_SP1) / 10, max(bgk_matrix_SP1) / 100, max(bgk_matrix_SP1)),
                     c( "darkblue", "cadetblue3", "lightcyan", "white", "lightcoral", "orange", "red"))
pdf("results/plots/Obs_data_plots/bgk_PE_SP1_20TF_1000g_2.pdf", width=12) 
h_bgk <- Heatmap(bgk_matrix_SP1, name="Gene expression",col=col_fun, cluster_rows = FALSE,row_names_side = "left", show_column_names=FALSE, show_row_dend = FALSE,column_dend_height = unit(3, "cm"), 
                 row_dend_width = unit(3, "cm"),heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 14)))
draw(h_bgk, heatmap_legend_side = "bottom")
dev.off()

x_df_SP1 <- as.data.frame(xk_matrix_SP1)
names(x_df_SP1) <- c("1","2","3","4")
rownames(x_df_SP1) <- seq(1:20)

pdf("results/plots/Obs_data_plots/xk_PE_SP1_20TF_1000g_2.pdf", width=3, height=6 )
h_xk <- Heatmap(as.matrix(x_df_SP1), name="SP1 knockout TF interaction",cluster_rows = FALSE,cluster_columns  = FALSE, row_title= "TF" ,show_row_dend = FALSE, show_column_dend = FALSE, row_names_side = "left",col=c("purple", "white"),heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 12)))
draw(h_xk, heatmap_legend_side = "bottom")
dev.off()

# Weights DP
w_est_SP1_df <- as.data.frame(w_est_SP1)
w_est_SP1_df$TF <- seq(1:20)

pdf("results/plots/Obs_data_plots/wk_PE_SP1_20TF_1000g_2.pdf", width=2, height=4) 
ggplot ( w_est_SP1_df, aes(as.factor(TF),w_est_SP1)) +geom_bar(stat = "identity", fill="lightsalmon") +ylim(c(0, 0.5)) +labs ( y= "DP weight", x="TF")+ theme_classic()+coord_flip()
dev.off()

# yg computation
rnaseq_SP1 <- monet_data@gene_exp$Sp1hyp_Flk_fpkm #y0
atacseq_SP1 <- as.data.frame(monet_data@atac_seq)
atacseq_SP1 <- atacseq_SP1[atacseq_SP1$gene %in% genes_SP1,]
atacseq_SP1 <- as.matrix(atacseq_SP1[,-1])
class(atacseq_SP1) <- "numeric"
yg_SP1 <- as.data.frame(monet_sim_yg(rnaseq_SP1, atacseq_SP1, bgk_matrix_SP1, xk_matrix_SP1,wk_matrix_SP1))

## Traceplot of key genes ##
# Add a column for gene names
yg_WT$gene <- genes_WT 
yg_SP1$gene <- genes_SP1

#Look matching genes
matching_genes <- intersect(genes_WT,genes_SP1)
matching_genes

# Select key genes that are in both data sets
key <- c( "Hist1h1a", "Cdh5", "Slk" , "Hist1h3f", "Cdca7","Crip1"  )
# Take y0 of key genes from the data
rnaseq_SP1 <- as.data.frame(rnaseq_SP1); rnaseq_SP1$gene <- genes_SP1
y0_key_SP1<- rnaseq_SP1[rnaseq_SP1$gene %in% key,]; y0_key_SP1<- y0_key_SP1[,-2]
rnaseq_WT <- as.data.frame(rnaseq_WT); rnaseq_WT$gene <- genes_WT
y0_key_WT<- rnaseq_WT[rnaseq_WT$gene %in% key,]; y0_key_WT<- y0_key_WT[,-2]


yg_key_WT <- cbind(y0_key_WT,yg_WT[yg_WT$gene %in% key, ]) 
yg_key_SP1 <- cbind(y0_key_SP1, yg_SP1[yg_SP1$gene %in% key, ])

# Add a column for type of data
yg_key_WT$obs <- "Wild type"
yg_key_SP1$obs <- "Sp1 perturbation"

names(yg_key_WT)<- c("ES","HB","HE1","HE2","HP","Gene", "Data")
names(yg_key_SP1)<- c("ES","HB","HE1","HE2","HP","Gene", "Data")

yg_key<- rbind(yg_key_WT, yg_key_SP1)

# Put the data in long format
yg_key_long <- yg_key %>% pivot_longer(cols=c("ES","HB","HE1","HE2","HP"), names_to='Time_point',values_to='Expression') 

pdf("results/plots/Obs_data_plots/key_trajectories.pdf",width = 8) 
ggplot(yg_key_long, aes(x=Time_point, y= Expression, color=Data, group=Data))+ geom_line()+
  geom_point()+theme_minimal()  + labs(title = "Gene expression trajectories in hematopoiesis", y = "Expression", x = "Time point")+ scale_x_discrete(limits=c("ES","HB","HE1","HE2","HP")) + facet_wrap(~Gene, ncol=2)+scale_color_manual(values=c("#54CE46", "#A749D5"))
dev.off()

##### Difference x_k heatmaps
x_dif_0.5 <- x_df_WT - x_df_SP1

pdf("results/plots/Obs_data_plots/xk_PE_0.5_dif.pdf", width=3, height=6 )
col_fun = colorRamp2(c(min(x_dif_0.5), 0, max(x_dif_0.5)),
                     c( "red", "yellow","green"))
h_xk <- Heatmap(as.matrix(x_dif_0.5), name="TF interaction difference",cluster_rows = FALSE,cluster_columns  = FALSE, row_title= "TF" ,show_row_dend = FALSE, show_column_dend = FALSE, row_names_side = "left",col=col_fun,heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 12)))
draw(h_xk, heatmap_legend_side = "bottom")
dev.off()

x_dif_2 <- x_df_WT - x_df_SP1
pdf("results/plots/Obs_data_plots/xk_PE_2_dif.pdf", width=3, height=6 )
col_fun = colorRamp2(c(min(x_dif_2), 0, max(x_dif_2)),
                     c( "red", "yellow","green"))
h_xk <- Heatmap(as.matrix(x_dif_2), name="TF interaction difference",cluster_rows = FALSE,cluster_columns  = FALSE, row_title= "TF" ,show_row_dend = FALSE, show_column_dend = FALSE, row_names_side = "left",col=col_fun,heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 12)))
draw(h_xk, heatmap_legend_side = "bottom")
dev.off()
