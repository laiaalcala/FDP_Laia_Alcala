# edit: 4 Jul 2020
# run: 4 Jul
library(tidyverse); library(monet)

#Create a data frame with the number of genes, TF, iterations, the warmup and alpha variables
par_df <- expand_grid(
  "gen" = 1000, #number of genes
  "tf" = 10, #number of transcription factors
  "iter" = 6000, #number of iterations
  "wup" = 2000, #warmup
  "alpha" = 1 )

gn_type <- "WT" #gene type assigned as Wild Type

data_path <- here::here("data/processed") #real processed data path

#Save the paths of the wild type processed data of the different data types in different variables
chip_path <- list.files(data_path, "ES_ChIP_binary", full.names = T) #ChIP-seq binary table file path
gene_path <- list.files(data_path,pattern = str_c("^ge_full.*", gn_type),full.names = T) #RNA-seq Wild Type file path
atac_path <-list.files(data_path, str_c("^ATAC.*", gn_type, ".*time.*"), full.names = T) #ATAC-seq Wild Type file path

#Create a new column containing a vector of the size of the number of TF (10). 
#Will use it to do the simulation estimating different number of TF: from 10 to 20
par_df$tf_est <- list(c(0, 10))

lapply(seq_len(nrow(par_df)), function(ip) { #Apply the ip function over the sequence of the row numbers of the data frame. 
#With this the function is applied in each row of the df at the same time which contain different variables for different simulations (in this case we only have one row).
  
  par_v <- par_df[ip, ] #take the row

  ntf <- par_v$tf #number of transcription factors of the current row
  itr <- par_v$iter #number of iterations of the current row
  wup <- par_v$wup #warmup of the current row
  ng <- par_v$gen #number of genes of the current row
  ap <- par_v$alpha #alpha of the current row

  # Set path for results
  gn_path <- here::here(str_c("results/tables/sim_paperv2", ng, "/"))
  #Create the directory of the path in case it does not exist and print a message when it has been created
  if (!dir.exists(gn_path)) {
    dir.create(gn_path)
    cat("New gene folder created \n")
  }

  #Data simulation using monet package. Inputs: character paths of the data types WT files, number of genes and number of TF
  monet_sim <- monet_sim(chip_path, atac_path, gene_path, n_g = ng, n_tf = ntf)

  #save the simulation in another file
  save(monet_sim, file = str_c(gn_path, "sim_data-120720_", ip, ".Rdata"))
  
  #Prepare simulation data. This function puts the simulated data in the format of interest for Stan
  monet_l <- prep_sim(monet_sim, list(de_b = 2, alpha = ap))

  i_test <- par_v$tf_est[[ip]] #take the vector of the TF estimations of the current row

  lapply(seq_along(i_test), function(i_est) { #Apply the i_est function over the SEQUENCE of the tf_est of the current row in a paralellized way.
    #This function.....
    tf_est <- i_test[i_est] #take the current sequence number, which belong to the number of transcription factors we want to add
    monet_l$no_tf <- monet_l$no_tf + tf_est #update the number of transcription factors on the simulated data

    cat("\n...run for ", ip, "-th row completed\n\n")
    
    #Fit the Monet model to the prepared simulated data and save the output in another variable.
    monet_smpl <- stan_monet.fit(standata = monet_l, chains = 4, iter = itr,refresh = 500, warmup = wup, 
                                 algorithm = "sampling", cores = 4)

    if (tf_est == 0) { #No estimated TF
      save(monet_smpl,file = str_c(gn_path, "sim_genes_", ng, "tf", ntf,"-a_", ap, ".Rdata"))
    
      } else {
      save(monet_smpl,file = str_c(gn_path, "tf_sim_genes_", ng, "tf", ntf,"-a_", ap, ".Rdata"))
    }

    rm(tmp, monet_l, monet_sim)
    gc()
  })
})
