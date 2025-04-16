data {
  int<lower=1> N; // Number of data points
  vector[N] x1; // Location of observations of input variable 1
  vector[N] x2; //Location of observations of input variable 2
  vector[N] x3; //Time
  vector[N] y; // Observation values
}

parameters {
  real<lower=0> sigma; // Standard deviation of the Gaussian process
  real<lower=0> length_scale1; // Length scale of K for input variable 1
  real<lower=0> length_scale2;
  real<lower=0> length_scale3;
  real<lower=0> p;
}

model {
  vector[N] mu = rep_vector(0,N); //Prior mean vector

  matrix[N, N] K; // Covariance matrix 

  // Compute the covariance matrix using a squared exponential kernel
  for (i in 1:N) {
    for (j in 1:N) {
      K[i, j] = ((exp(-square(x1[i] - x1[j]) / (2 * square(length_scale1)))) * (exp(-square(x2[i] - x2[j]) / (2 * square(length_scale2)))) * (exp((-2/square(length_scale3)) * square(sin(pi() * abs(x3[i]-x3[j])/p)))));
      if (i == j){
      K[i,j] = K[i,j] + 1e-7;
      }
    }
  }

  K = square(sigma) * K; //Add variance

  // Priors
  sigma ~ normal(0,1) T[0,];
  length_scale1 ~ gamma(2,1);
  length_scale2 ~ gamma(2,1);
  length_scale3 ~ gamma(2,1);
  p ~ normal(0.75,0.1) T[0,];

  // Likelihood
  y ~ multi_normal(mu, K);
}
