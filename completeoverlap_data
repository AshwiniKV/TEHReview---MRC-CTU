#Simulated data

p = 3; n = 100; n_burn = 2000; n_sim = 120

FOR i = 1,...,n_sim
  Set SEED as i
  x1 = x2 = x3 = N(0, 1)
  mu = 0.9x1 + 0.8x2 + 0.7x3
  pi = pnorm(mu)
  tau = (0.5(x1 >-3/4) + 0.25(x2 > 0.2) + 0.25(x3 > 3/4))
  ate_tau = mean(tau)
  y_noiseless = mu + tau z
  y = y_noiseless + rnorm(z


