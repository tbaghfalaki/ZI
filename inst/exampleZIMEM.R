data(long.data)

ZIMEM(
  FixedY = Y1 ~ x1 + x2, RandomY = ~obstime, GroupY = ~id,
  FixedZ = ~ x1 + x2, RandomZ = ~obstime, GroupZ = ~id,
  data = long.data, n.chains = 1, n.iter = 200, n.burnin = 100,
  n.thin = 1, family = "Poisson"
)

