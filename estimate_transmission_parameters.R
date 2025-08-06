source("simulation_analysis_gamma_test.R")

# Function to simulate a set of households
simulate_households <- function(n_households, ...) {
  do.call(rbind, lapply(1:n_households, function(i) sim.hh.func.fixed(N = i, ...)))
}

# Summary statistics: infection prevalence by role
summary_stats <- function(df) {
  tapply(df$infected, df$role, mean)
}

# Generate synthetic data with known parameters
set.seed(123)
true_params <- c(p.comm.base.infant.fix = 0.001, p.hh.base.infant = 0.1)
obs <- simulate_households(
  n_households = 50,
  p.comm.base.infant.fix = true_params[1],
  p.hh.base.infant = true_params[2]
)
obs_stats <- summary_stats(obs)

# Objective function to minimize
objective <- function(par) {
  sim <- simulate_households(
    n_households = 50,
    p.comm.base.infant.fix = par[1],
    p.hh.base.infant = par[2]
  )
  sim_stats <- summary_stats(sim)
  sum((sim_stats - obs_stats)^2)
}

# Initial guess and parameter bounds
init <- c(0.005, 0.2)
fit <- optim(
  par = init,
  fn = objective,
  method = "L-BFGS-B",
  lower = c(0, 0),
  upper = c(1, 1)
)

cat("True parameters:\n")
print(true_params)
cat("\nEstimated parameters:\n")
print(fit$par)
