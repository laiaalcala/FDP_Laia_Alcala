#' Simulation data based on monet
#'
#' @export
#'
#' @importFrom data.table fread
#' @importFrom rmutil rlaplace
#' @importFrom truncnorm rtruncnorm
monet_sim <- function(chip_path, atac_path, gene_path,
                      n_g, n_tf = 9, nt = 4, seed = as.numeric(Sys.time())) {
    
  #Set the seed of random number generator to ensure that the results are reproducible
    set.seed(seed)

    # Read ChipSeq data
    chip_d <- fread(chip_path) # [sample(.N, n_g)]

    # Sample n_g from ChipSeq data and convert to matrix
    cd_m <- cbind(
        chip_d[sample(.N, n_g)][, -"Gene"],
        chip_d[sample(.N, n_g)][, -"Gene"][, 1:1]
    ) %>%
        as.matrix()

    # Sample b_{g,k}from a a Laplace distribution with mean equal to 0 and standard deviation equal to 2
    r_lp <- rlaplace(n_g * n_tf, 0, 1) ####!!!!!!!!!!!!!
    #Values near 0.01 are changed to 0 
    r_lp[abs(r_lp) < 10^-2] <- 0
    # Add more zeroes
    i_set <- sample(1:(n_g * n_tf), tail(n_g * n_tf / 1.5), replace = TRUE)
    r_lp[i_set] <- 0
    # Convert the sample to matrix
    b_sim <- r_lp %>% matrix(nrow = n_g, ncol = n_tf)

    # Sample x_k from a truncated normal distribution with lower bound and mean equal to 0 and standard deviation equal to 2
    r_tn <- truncnorm::rtruncnorm(n_tf * nt, a = 0, mean = 0, sd = 1) ####!!!!!!!!!!!!!
    # Convert the sample to matrix
    x_sim <- r_tn %>% matrix(nrow = n_tf, ncol = nt)

    # Generate z_g by multiplying both simulations b_{g,k} and x_k.
    z <- b_sim %*% x_sim

    # Read ATACSeq data
    atac_wt <- fread(atac_path)

    # Read RNASeq data
    gene_dt <- fread(gene_path)

    gene_dt <-
        #Transform to its inverse hyperbolic sine
        cbind(gene_dt[, "gene"], asinh(gene_dt[, -"gene"])) %>% 
        # Reshape into long format and prepare for filtering where genes with lack of variation or very low expressed are removed
        geneDtMelt() %>%
        # Standardize the gene expression 
        geneStanderdise() %>%
        # Reshapes again the data into wide format
        dcast(gene ~ variable)

    # Samples from the first cell stage of expression data
    y0 <- gene_dt[sample(.N, n_g)][, 2][[1]]
    
    # Sample a_g from ATACSeq data in a matrix form
    atac_m <- as.matrix(atac_wt[, -"Gene"][sample(.N, n_g)])

    # Obtain y_g by  combining y_0, a_g, z_g and a normally distributed noise 
    y_sim <- y0 + atac_m * z
    y_sim <- y_sim + matrix(rnorm(length(y_sim), 0, 1), nrow(y_sim), ncol(y_sim))

    # Generate a list containing all the results
    sim_l <-
        list(
            no_gns = n_g,
            no_tf = n_tf,
            no_tpt = nt,
            y = t(y_sim),
            b_sim = b_sim,
            x_sim = x_sim,
            y0 = y0,
            a = t(atac_m)
        )
    
    # Return the list
    return(sim_l)
}

#' Prepare simulation data for monet
#'
#' @export
prep_sim <- function(sim_l, par = list()) {
    par_def <- list(alpha = 4, de_b = 2, sigma = 1, sigma_x = 4)

    par <- c(par_def[setdiff(names(par_def), names(par))], par)


    monet_dat <-
        c(
            sim_l[c("no_gns", "no_tf", "no_tpt", "y", "y0", "a")],
            par
        )

    return(monet_dat)
}
