# edit: 16 May 2023
library(monet)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggplot2); library(ggpubr); library(viridis);library(gridExtra)
library(data.table)
library(matrixcalc) # for the hardmard product
library(ComplexHeatmap); library(circlize)
library(rstan)

# ## MCMC sampling output of the STAN model fitting
load("~/Downloads/monet-applied/results/tables/sim_laia_2000g/sim_genes_2000g30TF-a_1_8000it_2000wp.Rdata")

# Create a new monet sample object
monet_res <- new("monetSampling", monet_smpl)

# remove large object for efficiency
rm(monet_smpl)
str(monet_res, 2)

###  Traceplots MCMC
pdf("results/plots/Datafitting_plots/2000g/MCMC_trace_logpost_8000it_2000wup_30TF.pdf", width=25) 
traceplot(monet_res@stanfit[[1]], pars="lp__", size = 0.5, alpha = 0.6) +
labs(title = "MCMC log posterior traceplot",subtitle ="2000 genes - 30TFs - 8000 iterations - 2000 w.up",x = "Iterations", y= "Log posterior")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_text(size = rel(3.2)),
        axis.title.x = element_text(size = rel(3.2)),
        plot.title = element_text(size = rel(3.2)),
        plot.subtitle=element_text(size=rel(3.6)),
        legend.text=element_text(size=rel(2.8)),
        legend.title=element_text(size=rel(2.5)),
        legend.key.size = unit(2, 'cm'))
dev.off()

pdf("results/plots/Datafitting_plots/2000g/MCMC_trace_w__30TF_8000it_2000wup.pdf", width=25) 
traceplot(monet_res@stanfit[[1]], pars="w", size = 0.5, alpha = 0.6) 
dev.off()

pdf("results/plots/Datafitting_plots/2000g/MCMC_trace_w__30TF_8000it_2000wup_5.pdf", width=25) 
traceplot(monet_res@stanfit[[1]], pars=c("w[6]","w[8]","w[19]", "w[21]","[26]","w[30]") , size = 0.5, alpha = 0.6) 
dev.off()

### Densities of the parameter values with a line in the posterior mean
pars=c("x[2,6]", "x[3,29]" ,"x[4,8]","x[4,26]","b[1,1]", "b[241,3]","b[675,6]", "b[836,19]","w[6]", "w[10]", "w[19]", "w[22]")

pdf("results/plots/Datafitting_plots/2000g/dens_stanfit_30TF_8000it_2000wup_2.pdf")
stan_dens(monet_res@stanfit[[1]],ncol=3,fill="#CAEEF5", pars=pars)+ 
  ggtitle("Density Estimates") + theme_classic()
dev.off()

list_plots=list()
pdf("results/plots/Datafitting_plots/2000g/dens_stanfit_grid.pdf")
for (e in pars){
  all <- rstan::get_posterior_mean(monet_res@stanfit[[1]], par = e)
  est <- all[, ncol(all)]
  print(est)
  d <- stan_dens(monet_res@stanfit[[1]],ncol=3,fill="#CAEEF5", pars=e)+ geom_vline(aes(xintercept = est))+theme_classic() +ylab("")
  list_plots[[e]] =d
}
gridExtra::grid.arrange(grobs = list_plots, ncol=4,top="Densities of the posterior point estimates", common.legend = TRUE, legend="bottom")
dev.off()

# p12 <- stan_dens(monet_res@stanfit[[1]],ncol=3,fill="#CAEEF5", pars="w[22]")+ geom_vline(aes(xintercept = 0.002808838))+theme_classic() +ylab("")   
# 
# pdf("results/plots/Datafitting_plots/2000g/dens_stanfit_grid.pdf", width=14)
# grid.arrange(p1, p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, nrow = 3, top= "Densities of the posterior point estimates")
# dev.off()

#######

# Extract posterior means
w_all <- rstan::get_posterior_mean(monet_res@stanfit[[1]], par = "w")
b_all <- rstan::get_posterior_mean(monet_res@stanfit[[1]], par = "b")
x_all <- rstan::get_posterior_mean(monet_res@stanfit[[1]], par = "x")

# Extract posterior across all chains (last column)
w_est <- w_all[, ncol(w_all)]
b_est <- b_all[, ncol(b_all)]
x_est <- x_all[, ncol(x_all)]

wk <- as.matrix(diag(w_est)) # Dimension: TF (10) x TF (10)
bgk <- as.data.frame(b_est) 
xk <- as.data.frame(x_est)

xk_matrix <- matrix(NA,ncol=30, nrow=4 )

c<-0
for (i in seq(nrow(xk)/30)){
  c <- c+1
  xk_matrix[i,] <- xk[c:(c+29),]
  c <- c+29
}

bgk_matrix <- matrix(NA,ncol=30, nrow=2000 )
c<-0
for (i in seq(nrow(bgk)/30)){
  c <- c+1
  bgk_matrix[i,] <- bgk[c:(c+29),]
  c <- c+29
}

xk_matrix <- t(xk_matrix) # Dimension: TF(10) x Time points (4)
bgk_matrix <- t(bgk_matrix)  # Dimension: TF (10) x Genes (2000) 

#### Plots ####
  # Look at the density distributions of the parameters
  # Density plot for wk
pdf("results/plots/Datafitting_plots/2000g/wk_density_MCMC_30TF_8000it_2000wup.pdf") 
ggplot(as.data.frame(w_est), aes(x=w_est)) + 
  geom_density(color="black", fill="lightsalmon")+ 
  geom_vline(aes(xintercept=mean(w_est)), color="black", linetype="dashed", size=1)+
  labs(title="Density distribution of w of MCMC sampling of the simulated data",x="wk", y = "Count")+
  theme_classic() + xlim(0,0.75)
dev.off()

  # Density plot for bgk
pdf("results/plots/Datafitting_plots/2000g/bgk_density_MCMC_30TF_8000it_2000wup.pdf") 
ggplot(bgk, aes(x=b_est)) + 
  geom_density(color="black", fill="#80CBC4")+ 
  geom_vline(aes(xintercept=mean(b_est)), color="black", linetype="dashed", size=1)+xlim(-5,5)+
  labs(title="Density distribution of b of MCMC sampling of the simulated data",x="bgk", y = "Count")+
  theme_classic()
dev.off()

  # Density plot for xk
pdf("results/plots/Datafitting_plots/2000g/xk_density_MCMC_30TF_8000it_2000wup.pdf") 
ggplot(xk, aes(x=x_est)) + 
  geom_density(color="black", fill="#9990BA")+ 
  geom_vline(aes(xintercept=mean(x_est)), color="black", linetype="dashed", size=1)+xlim(0,10)+
  labs(title="Density distribution of x of MCMC sampling of the simulated data",x="xk", y = "Count")+
  theme_classic()
dev.off()

  # Weights DP
w_est_df <- as.data.frame(w_est)
w_est_df$TF <- seq(1:30)

pdf("results/plots/Datafitting_plots/2000g/wk_PE_30TF_8000it_2000wup.pdf", width=2, height=4) 
ggplot ( w_est_df, aes(as.factor(TF),w_est)) +geom_bar(stat = "identity", fill="lightsalmon") +ylim(c(0, 0.5)) +labs ( y= "DP weight", x="TF")+ theme_classic()+coord_flip()
dev.off()

  #Heatmaps of B and X
names(bgk_matrix) <- NULL
rownames(bgk_matrix) <- seq(1:30)
col_fun = colorRamp2(c(-4, 0, 4), c("green", "white", "blue"))
pdf("results/plots/Datafitting_plots/2000g/bgk_PE_30TF_8000it_2000wup.pdf", width=12) 
h_bgk <- Heatmap(bgk_matrix, name="Gene expression",col=col_fun, cluster_rows = FALSE,row_names_side = "left",show_column_dend = FALSE, show_row_dend = FALSE,heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 14)))
draw(h_bgk, heatmap_legend_side = "bottom")
dev.off()
##
col_fun = colorRamp2(c(min(bgk_matrix),min(bgk_matrix) / 10, min(bgk_matrix) / 100, 0,max(bgk_matrix) / 10, max(bgk_matrix) / 100, max(bgk_matrix)),
                     c( "darkblue", "cadetblue3", "lightcyan", "white", "lightcoral", "orange", "red"))
pdf("results/plots/Datafitting_plots/2000g/bgk_PE_30TF_8000it_2000wup_2.pdf", width=12) 
h_bgk <- Heatmap(bgk_matrix, name="Gene expression",col=col_fun, cluster_rows = FALSE,row_names_side = "left", show_row_dend = FALSE,show_row_names=TRUE,show_column_dend=FALSE,heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 14)))
draw(h_bgk, heatmap_legend_side = "bottom")
dev.off()
##

x_est_df <- as.data.frame(xk_matrix)
names(x_est_df) <- c("1","2","3","4")
rownames(x_est_df) <- seq(1:30)

pdf("results/plots/Datafitting_plots/2000g/xk_PE_30TF_8000it_2000wup.pdf", width=3, height=6 )
h_xk <- Heatmap(as.matrix(x_est_df), name="TF interaction",cluster_rows = FALSE,cluster_columns  = FALSE, row_title= "TF" ,show_row_dend = FALSE, show_column_dend = FALSE, row_names_side = "left",col=c("purple", "white"),heatmap_legend_param = list(legend_direction = "horizontal",labels_gp = gpar(fontsize = 12)))
draw(h_xk, heatmap_legend_side = "bottom")
dev.off()


######### Compute yg ##########
# Function to compute y_g 
# Inputs: y_0 and a_g from the simulated data set
#         b_g_k, x_k and w_k from the fitting of the monet STAN model
monet_sim_yg <- function(y0, ag, bgk, xk,wk){
  zg <- t(wk %*% bgk) %*% xk
  yg <- y0 + hadamard.prod(ag, zg)
  return(yg)
  #save(yg,file = str_c("results/tables/sim_laia_2000g/yg_2000g_10TF.data"))
}

no_tf=30 # number of TF
no_tpt=4 # number of time points (without y0)

load("~/monet-applied/results/tables/sim_laia_2000g/sim_data-8000it_2000wp_30TF.Rdata")
y0<- monet_sim$y0 # take y_0 from the simulated dataset
ag<- t(monet_sim$a) # take a_g from the simulated dataset

yg <- monet_sim_yg(y0,ag,bgk_matrix,xk_matrix, wk )
yg <- data.table(yg)
yg$gene <- seq(1,nrow(yg))

yg_stand <- monet:::geneDtMelt(yg) %>% monet:::geneStanderdise() %>% dcast(gene ~ variable)
 
  #keep only the genes that are in both datasets --> 1024 genes
y_sim_filt <- data_frame()
yg_stand_filt <- data_frame()
for (x in y_sim$gene){
  if (x %in% c(yg_stand$gene)){
    y_sim_filt<-rbind(y_sim_filt, y_sim[y_sim$gene==x, ])
    yg_stand_filt<-rbind(yg_stand_filt, yg_stand[yg_stand$gene==x, ])
    }
}

  #compute the sum of differences of gene expression of the matching genes of the simulated and estimated
y_sim_filt_2 <- as.data.frame(y_sim_filt[,3:6])
yg_stand_filt_2 <- as.data.frame(yg_stand_filt[,2:5])
dif =c()
for (x in seq(1:nrow(y_sim_filt_2))){
  d<-0
  for (y in seq(1:ncol(y_sim_filt_2))){
    d<- d+ abs(y_sim_filt_2[x,y]-yg_stand_filt_2[x,y])
    #print(d)
  }
  dif <-append(dif,d)
}

  # take the 10 genes with the smallest difference
indexes <- seq(1:length(dif))
min_indexes <- c()
for (i in seq(1:10)){
  min_indexes <- append(min_indexes,which.min(dif))
  dif <- dif[-which.min(dif)]
}

min_genes <- yg_stand_filt[min_indexes, gene]

  # Plot of the same 10 trajectories
yg_10 <- yg_stand_filt[yg_stand_filt$gene %in% min_genes, ] # take the selected genes
names(yg_10)<- c("gene","1","2","3","4")

# yg_10_long <- yg_10 %>% pivot_longer(cols=c("1","2","3","4"), names_to='Time_point',values_to='Expression') 
# yg_10_long$Time_point = as.numeric(yg_10_long$Time_point)

sim_10 <- y_sim[y_sim$gene %in% min_genes, ] # take the selected genes
sim_10<-sim_10[,-2] # remove timepoint 0
names(sim_10)<- c("gene","1","2","3","4")

yg_10$data<- "Estimated"
sim_10$data<- "Simulated"
traj_10 <- rbind(yg_10, sim_10, fill=TRUE)

traj_10_long <- traj_10 %>% pivot_longer(cols=c("1","2","3","4"), names_to='Time_point', values_to='Expression')
traj_10_long$Time_point = as.numeric(traj_10_long$Time_point)

pdf("results/plots/Datafitting_plots/2000g/10_trajectories.pdf",width = 8) 
ggplot(traj_10_long, aes(x=Time_point, y= Expression, color=data, group=data))+ geom_line()+theme_classic()  + labs(title = "Gene expression trajectories", y = "Expression", x = "Time point")+ scale_x_discrete(limits=c("HB","HE1","HE2","HP")) + facet_wrap(~gene, ncol=2)+scale_color_manual(values=c("#808080", "#CC0000"))
dev.off()
