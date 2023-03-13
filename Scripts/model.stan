data {
  // dimensions
  int<lower=1> no_gns; // # genes
  int<lower=1> no_tf; // # TFs
  int<lower=1> no_tpt; // # time points
  real<lower=0> de_b; //Laplace parameter
  real<lower=0> sigma; //Normal parameter
  real<lower=0> sigma_x; //NormalTruncated parameter
  // response variable
  vector[no_gns] y0; //Gene expression vector of genes at time 0
  vector[no_gns] y[no_tpt]; // response variable
  vector[no_gns] a[no_tpt]; // masking matrix (e.g. ATACseq)
  real<lower=0> alpha;
}

parameters {
  // restrict the model to take the choice of
  vector<lower=0>[no_tf] x[no_tpt]; //Shared latent transcription factor activities
  matrix[no_gns, no_tf] b; //Gene-specific transcription factor interaction network with genes as rows and TFs as columns
  simplex[no_tf] w; //Infinite vector of weights. It is a vector with non-negative values whose entries sum to 1
}

transformed parameters {
  vector[no_tpt] tmp[no_gns]; //Gene-specific latent trajectory
  vector[no_tpt] mu[no_gns]; //Gene expression time-course measurement
  //iterate over time points
  for (i_tpt in 1:no_tpt) {
    //model equations
    tmp[i_tpt] = b * (w .* x[i_tpt]);
    mu[i_tpt] = y0 + a[i_tpt] .* tmp[i_tpt];
  }
}

model {
  //PRIORS
  //iterate over time points and TFs
  for (i_tpt in 1:no_tpt) {
    for (i_tf in 1:no_tf) {
      x[i_tpt, i_tf] ~ normal(0, sigma_x)T[0, ]; //Truncated normal distribution prior
    }
  }
 
  w ~ dirichlet(rep_vector(alpha, no_tf)); //Dirichlet prior
  
  //iterate over TFs
  for (i_tf in 1:no_tf) {
    b[, i_tf] ~ double_exponential(0, de_b); //Double exponential (Laplace) distribution prior
  }
  
  //LIKELIHOOD
  //iterate over time points
  for (i_tpt in 1:no_tpt) {
    y[i_tpt] ~ normal(mu[i_tpt], sigma);
  }
}
