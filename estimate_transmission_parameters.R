# Chain binomial parameter estimation with age-dependent risks

log1mexp <- function(x) { ifelse(x < log(0.5), log1p(-exp(x)), log(-expm1(x))) }

# Negative log-likelihood: community and household transmission depend on age
negll <- function(par, dat, eps = 1e-10) {
  delta0 <- par[1]; delta1 <- par[2]
  alpha0 <- par[3]
  gamma  <- c(0, par[4:6])
  beta   <- c(0, par[7:9])

  age <- dat$agegrp
  p_comm <- exp(delta0 + gamma[age] + delta1 * dat$cases)
  p_hh   <- exp(alpha0 + beta[age])

  p_comm <- pmin(pmax(p_comm, eps), 1 - eps)
  p_hh   <- pmin(pmax(p_hh,   eps), 1 - eps)

  p_tot  <- 1 - (1 - p_comm) * (1 - p_hh) ^ dat$n_inf
  p_tot  <- pmin(pmax(p_tot,  eps), 1 - eps)

  -sum(dat$event * log(p_tot) + (1 - dat$event) * log1mexp(log(p_tot)))
}

set.seed(123)
n <- 1000
true_params <- c(
  delta0 = -5.5,
  delta1 = 0.1,
  alpha0 = -2,
  gamma_sibling = 0.4,
  gamma_parent = 0.2,
  gamma_elder = 0.6,
  beta_sibling = 0.3,
  beta_parent = 0.1,
  beta_elder = 0.5
)

dat <- data.frame(
  agegrp = sample(1:4, n, replace = TRUE),
  cases  = runif(n),
  n_inf  = sample(0:3, n, replace = TRUE)
)

gamma <- c(0, true_params["gamma_sibling"], true_params["gamma_parent"], true_params["gamma_elder"])
beta  <- c(0, true_params["beta_sibling"], true_params["beta_parent"], true_params["beta_elder"])
p_comm <- exp(true_params["delta0"] + gamma[dat$agegrp] + true_params["delta1"] * dat$cases)
p_hh   <- exp(true_params["alpha0"] + beta[dat$agegrp])
eps <- 1e-10
p_comm <- pmin(pmax(p_comm, eps), 1 - eps)
p_hh   <- pmin(pmax(p_hh,   eps), 1 - eps)
p_tot  <- 1 - (1 - p_comm) * (1 - p_hh) ^ dat$n_inf
p_tot  <- pmin(pmax(p_tot,  eps), 1 - eps)
dat$event <- rbinom(n, 1, p_tot)

init <- rep(0, 9)
fit <- optim(init, negll, dat = dat, method = "L-BFGS-B")

cat("True parameters:\n")
print(true_params)
cat("\nEstimated parameters:\n")
print(fit$par)

