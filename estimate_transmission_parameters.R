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

# Simulate dataset from the household transmission model ---------------------
source("simulation_analysis_gamma_test.R")

# True parameter values used to generate the data
true_params <- c(
  p.comm.base.infant.fix    = 0.001,
  p.comm.multiplier.sibling = 1.5,
  p.comm.multiplier.parent  = 0.8,
  p.comm.multiplier.elder   = 0.5,
  p.hh.base.infant          = 0.1,
  p.hh.multiplier.sibling   = 0.4,
  p.hh.multiplier.parent    = 0.2,
  p.hh.multiplier.elder     = 0.3
)

# Simulate multiple households and combine into a single dataset
sim_list <- lapply(1:200, function(hh) {
  do.call(sim.hh.func.fixed, c(list(N = hh), as.list(true_params)))
})
sim_dat <- do.call(rbind, sim_list)

# Prepare data for likelihood calculation
role_map <- c(infant = 1, sibling = 2, adult = 3, elder = 4)
sim_dat$agegrp <- role_map[sim_dat$role]

# Number of other infected household members
hh_infections <- ave(sim_dat$infected, sim_dat$HH, FUN = function(x) sum(x))
sim_dat$n_inf <- pmax(hh_infections - sim_dat$infected, 0)

sim_dat$cases <- 0  # community cases placeholder
dat <- sim_dat[, c("agegrp", "cases", "n_inf", "infected")]
names(dat)[4] <- "event"

# Fit parameters --------------------------------------------------------------
init <- rep(0, 9)
fit <- optim(init, negll, dat = dat, method = "L-BFGS-B")

cat("True parameters:\n")
print(true_params)
cat("\nEstimated parameters:\n")
print(fit$par)

