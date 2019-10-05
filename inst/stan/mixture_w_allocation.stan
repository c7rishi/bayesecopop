// Finite mixture model (Overfitting by superimpose)
// Return posterior probabilities
// of allocation for each point

data {
  int n_groups; //number of components
  int n_data; //number of data points
  real y[n_data, 2]; //data
  vector[n_groups] alpha;
}

parameters {
  real mu[n_groups];
  real<lower = 0> sigma[n_groups];

  real nu[n_groups];
  real<lower = 0> kappa[n_groups];

  simplex[n_groups] pmix;
}

model {
  real contributions[n_groups];

  // priors
  mu ~ cauchy(0, 10);
  sigma ~ cauchy(0, 10);

  nu ~ cauchy(0, 10);
  kappa ~ cauchy(0, 10);

  pmix ~ dirichlet(alpha);


  // likelihood
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
    target += log_sum_exp(contributions);
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

