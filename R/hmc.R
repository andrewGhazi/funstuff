#### Let's just try writing an HMC from scratch for a 2D Gaussian ----

library(tidyverse)
library(mvnfast)

# T_optimal = H(q,p)^((2-beta) / (2beta)) # see pg 34 of Betancourt, Conceptual intro to HMC
# for a Gaussian target, beta = 2
# so for my test T_opt is directly proportional to the Hamiltonian

# K = -log(pi(p|q))
# V = -log(pi(q))
# K = 1/2 t(p) %*% solve(M) %*% p + log(det(M)) + cons # Euclidean-Gaussian energy

# dq/dt = +dK/dp (partials)
# dp/dt = -dK/dq - dV/dq

# Going to assume the mass matrix M = diag(2) for now
# so pi(p|q) = mvn(p | 0, Sigma = M)



H = function(pos, mom){
  -dmvn(mom, mu = c(0,0), sigma = diag(2), log = TRUE) + #K(p|q)
    -dmvn(pos, mu = c(0,0), sigma = diag(2), log = TRUE) # V(q)
}
dVdq = function(pos){
  c(2*pos[1], 2*pos[2]) # I think this is right for a paraboloid
}

n_iter = 10000

pos = rmvn(1, mu = c(0,0), sigma = diag(2)); print(pos)
mom = rmvn(1, mu = c(0,0), sigma = diag(2))
eps = .1 # just a guess

hmc_mat = matrix(ncol = 4,
                 nrow = n_iter)
hmc_mat[1,1:2] = pos
hmc_mat[1,3:4] = mom

for (i in 2:n_iter){
  new_mom = rmvn(1, mu = c(0,0), sigma = diag(2)) 
  
  iter_H = H(hmc_mat[i-1,1:2], new_mom )
  T_opt = 1/2 * iter_H
  n_lf = floor(T_opt / eps) + 1
  lf_mat = matrix(ncol = 4,
                  nrow = n_lf)
  
  lf_mat[1,1:2] = hmc_mat[i-1,1:2]
  lf_mat[1,3:4] = new_mom
  for (n in 2:n_lf){ # leapfrog integrator loop
    mom_half = lf_mat[n-1,3:4] - eps/2 * dVdq(lf_mat[n-1,1:2])
    lf_mat[n,1:2] = lf_mat[n-1,1:2] + eps * mom_half
    lf_mat[n,3:4] = mom_half - eps/2*dVdq(lf_mat[n,1:2])

  }
  # lf_mat[,1:2] %>%
  #   as_tibble() %>%
  #   mutate(lf_step = 1:n()) %>%
  #   ggplot(aes(V1, V2)) +
  #   geom_point(aes(color = lf_step)) +
  #   coord_equal() +
  #   scale_color_viridis_c() +
  #   theme_light()
  lf_draw = sample(1:n_lf, size = 1)
  pos_prop = lf_mat[lf_draw,1:2]
  mom_prop = -lf_mat[lf_draw,3:4] # negation necessary for reversibility
  
  Q_iter = exp(-iter_H)
  prop_H = H(pos_prop, mom_prop)
  Q_prop = exp(-prop_H)
  
  if (Q_prop > Q_iter){
    mh_check = TRUE
  } else {
    mh_prob = exp(prop_H - iter_H)
    mh_check = runif(1) < mh_prob
  }
  # mh_prob = min(1, Q_prop / Q_iter)
  
  
  if (mh_check){
    hmc_mat[i,1:2] = pos_prop
    hmc_mat[i,3:4] = mom_prop
  } else{
    hmc_mat[i,1:2] = hmc_mat[i-1,1:2]
    hmc_mat[i,3:4] = hmc_mat[i-1,3:4]
  }
}

hmc_mat %>%
  as_tibble() %>% 
  ggplot(aes(V1, V2)) +
  geom_point(alpha = .1) + 
  theme_light() + 
  coord_equal()



#### Now let's adapt the HMC from BDA3 for LR ----

log_post_theta = function(theta, y, x){
  alpha = theta[1]
  beta = theta[2]
  lsigma = theta[3]
  
  log_prior = dnorm(alpha, sd = 10, log = TRUE) + # N(0, 10)
    dnorm(beta, sd = 1, log = TRUE) + # N(0,1)
    dexp(exp(lsigma), rate = 1, log = TRUE) # sigma ~ Exp(1)
  
  ll = sum(dnorm(y, mean = alpha + beta*x, sd = exp(lsigma), log = TRUE))
  
  return(log_prior + ll)
}


gradient_th_numerical = function(theta, y, x){
  alpha = theta[1]
  beta = theta[2]
  lsigma = theta[3]
  sigma = exp(lsigma)
  
  d = length(theta)
  e = .0001
  diff = rep(NA, d)
  for (k in 1:d){
    th_hi = theta
    th_lo = theta
    th_hi[k] = theta[k] + e
    th_lo[k] = theta[k] + e
    diff[k] = (log_post_theta(th_hi, y, x) - log_post_theta(th_lo, y, x)) / (2*e)
  }
  
  return(diff)
}