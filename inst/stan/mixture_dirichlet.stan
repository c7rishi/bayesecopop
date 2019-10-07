// Mixture model as approximate Dirichlet process
// Does not return posterior probabilities of allocation for each point

data {
  int n_groups; //number of components
  int n_data; //number of data points
  real y[n_data, 2]; //data
  real<lower=0> alpha0;
  real<lower=0> period_over_2pi;
}

parameters {
  real mu[n_groups];
  real<lower = 0> sigma[n_groups];

  real nu[n_groups];
  real<lower = 0> kappa[n_groups];

  real <lower=0,upper=1> v[n_groups]; // stick breaking v
}

transformed parameters{
  simplex [n_groups] pmix;
  pmix[1] = v[1];
  // stick-break process based on The BUGS book Chapter 11 (p.294)
  for(j in 2:(n_groups-1)){
      pmix[j] = v[j]*(1-v[j-1]) * pmix[j-1]/v[j-1];
  }
  pmix[n_groups]=1-sum(pmix[1:(n_groups-1)]); // to make a simplex.
}

model {
  real contributions[n_groups];

  // priors
  mu ~ cauchy(0, 10);
  sigma ~ cauchy(0, 10);

  nu ~ cauchy(0, 10);
  kappa ~ cauchy(0, 10);

  v ~ beta(1, alpha0);


  // likelihood
  for(i in 1:n_data) {
    for(k in 1:n_groups) {
      if (kappa[k] < 100) {
        contributions[k] = log(pmix[k]) +
          normal_lpdf(y[i, 1] | mu[k], sigma[k]) +
          von_mises_lpdf(y[i, 2]/period_over_2pi | nu[k], kappa[k]);
      } else {
        // for stability, when kappa is large,
        // use normal distribution
        contributions[k] = log(pmix[k]) +
          normal_lpdf(y[i, 1] | mu[k], sigma[k]) +
          normal_lpdf(y[i, 2]/period_over_2pi | nu[k], 1/sqrt(kappa[k]));
      }
    }
    target += log_sum_exp(contributions);
    // log jacobian unnecessary since just a scale transform
  }
}


generated quantities {
  real log_lik[n_data];
  real contributions[n_groups];

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
  }
}

