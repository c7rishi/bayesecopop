// Mixture model as approximate Dirichlet process
// Returns posterior probabilities of allocation for each point

data {
  int n_groups; //number of components
  int n_data; //number of data points
  real y[n_data, 2]; //data
  real<lower=0> alpha0;
}

parameters {
  real mu[n_groups];
  real<lower = 0> sigma[n_groups];

  real nu[n_groups];
  real<lower = 0> kappa[n_groups];

  real<lower=0,upper=1> v_raw[n_groups-1]; // stick breaking v
}

transformed parameters{
  simplex[n_groups] pmix;
  vector[n_groups] pmix_raw;
  real<lower=0,upper=1> v[n_groups];

  v[1:(n_groups-1)] = v_raw[:];
  v[n_groups] = 1;

  pmix_raw[1] = v[1];
  // stick-break process based on The BUGS book Chapter 11 (p.294)
  for(j in 2:n_groups){
      pmix_raw[j] = v[j]*(1-v[j-1]) * pmix_raw[j-1]/v[j-1];
  }
  pmix = pmix_raw/sum(pmix_raw[:]); // to make a simplex.
}

model {
  real contributions[n_groups];

  // priors
  mu ~ cauchy(0, 10);
  sigma ~ cauchy(0, 10);

  nu ~ cauchy(0, 10);
  kappa ~ cauchy(0, 10);

  v_raw ~ beta(1, alpha0);


  // likelihood
  for(i in 1:n_data) {
    for(k in 1:n_groups) {
         if (kappa[k] < 100) {
        contributions[k] = log(pmix[k]) +
          normal_lpdf(y[i, 1] | mu[k], sigma[k]) +
          von_mises_lpdf(y[i, 2] | nu[k], kappa[k]);
      } else {
        // for stability, when kappa is large,
        // use normal distribution
        contributions[k] = log(pmix[k]) +
          normal_lpdf(y[i, 1] | mu[k], sigma[k]) +
          normal_lpdf(y[i, 2] | nu[k], 1/sqrt(kappa[k]));
      }
    }
    target += log_sum_exp(contributions);
    // log jacobian unnecessary since just a scale transform
  }
}


generated quantities {
  real log_lik[n_data];
  vector[n_groups] contributions;
  simplex[n_groups] post_prob[n_data];
  int allocation[n_data];

  for(i in 1:n_data) {
    for(k in 1:n_groups) {
      if (kappa[k] < 100) {
        contributions[k] = log(pmix[k]) +
          normal_lpdf(y[i, 1] | mu[k], sigma[k]) +
          von_mises_lpdf(y[i, 2] | nu[k], kappa[k]);
      } else { // for stability
        contributions[k] = log(pmix[k]) +
          normal_lpdf(y[i, 1] | mu[k], sigma[k]) +
          normal_lpdf(y[i, 2] | nu[k], 1/sqrt(kappa[k]));
      }
    }
    log_lik[i] = log_sum_exp(contributions);
    post_prob[i] = softmax(contributions);
    allocation[i] = categorical_rng(post_prob[i]);
  }
}

