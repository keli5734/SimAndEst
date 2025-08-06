# Chain binomial parameter estimation with age-dependent risks

log1mexp <- function(x) { ifelse(x < log(0.5), log1p(-exp(x)), log(-expm1(x))) }

# Negative log-likelihood: estimate baseline probabilities and age multipliers
negll <- function(par, dat, eps = 1e-10) {
  p_comm_base <- par[1]
  p_comm_mult <- c(1, par[2:4])
  p_hh_base   <- par[5]
  p_hh_mult   <- c(1, par[6:8])

  age <- dat$agegrp
  p_comm <- p_comm_base * p_comm_mult[age]
  p_hh   <- p_hh_base   * p_hh_mult[age]

  p_comm <- pmin(pmax(p_comm, eps), 1 - eps)
  p_hh   <- pmin(pmax(p_hh,   eps), 1 - eps)

  p_tot  <- 1 - (1 - p_comm) * (1 - p_hh) ^ dat$n_inf
  p_tot  <- pmin(pmax(p_tot, eps), 1 - eps)

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

dat <- sim_dat[, c("agegrp", "n_inf", "infected")]
names(dat)[3] <- "event"

# Fit parameters --------------------------------------------------------------
init <- c(0.01, 1, 1, 1, 0.1, 1, 1, 1)
lower <- c(0, rep(0,3), 0, rep(0,3))
upper <- c(1, rep(5,3), 1, rep(5,3))
fit <- optim(init, negll, dat = dat, method = "L-BFGS-B", lower = lower, upper = upper)

cat("True parameters:\n")
print(true_params)
cat("\nEstimated parameters:\n")
names(fit$par) <- names(true_params)
print(fit$par)

