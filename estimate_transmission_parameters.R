source("simulation_analysis_gamma_test.R")

# Function to simulate a set of households
simulate_households <- function(n_households, ...) {
  do.call(rbind, lapply(1:n_households, function(i) sim.hh.func.fixed(N = i, ...)))
}

# Summary statistics: infection and detection prevalence by role
# Ensures all roles are represented; missing roles return 0
summary_stats <- function(df) {
  roles <- c("infant", "sibling", "adult", "elder")
  stats <- sapply(roles, function(r) {
    subset <- df[df$role == r, ]
    infected <- if (nrow(subset) == 0) 0 else mean(subset$infected)
    detected <- if (nrow(subset) == 0) 0 else mean(subset$detected.infected)
    c(infected = infected, detected = detected)
  })
  c(stats["infected", ], stats["detected", ])
}

# Generate synthetic data with known parameters
set.seed(123)
true_params <- c(
  p.comm.base.infant.fix = 0.001,
  p.comm.multiplier.sibling = 1.2,
  p.comm.multiplier.parent = 0.8,
  p.comm.multiplier.elder = 1.5,
  p.hh.base.infant = 0.1,
  p.hh.multiplier.sibling = 1.3,
  p.hh.multiplier.parent = 0.7,
  p.hh.multiplier.elder = 1.1
)
obs <- simulate_households(
  n_households = 50,
  p.comm.base.infant.fix = true_params["p.comm.base.infant.fix"],
  p.comm.multiplier.sibling = true_params["p.comm.multiplier.sibling"],
  p.comm.multiplier.parent = true_params["p.comm.multiplier.parent"],
  p.comm.multiplier.elder = true_params["p.comm.multiplier.elder"],
  p.hh.base.infant = true_params["p.hh.base.infant"],
  p.hh.multiplier.sibling = true_params["p.hh.multiplier.sibling"],
  p.hh.multiplier.parent = true_params["p.hh.multiplier.parent"],
  p.hh.multiplier.elder = true_params["p.hh.multiplier.elder"]
)
obs_stats <- summary_stats(obs)

# Objective function to minimize
objective <- function(par) {
  sim <- simulate_households(
    n_households = 50,
    p.comm.base.infant.fix = par[1],
    p.comm.multiplier.sibling = par[2],
    p.comm.multiplier.parent = par[3],
    p.comm.multiplier.elder = par[4],
    p.hh.base.infant = par[5],
    p.hh.multiplier.sibling = par[6],
    p.hh.multiplier.parent = par[7],
    p.hh.multiplier.elder = par[8]
  )
  sim_stats <- summary_stats(sim)
  if (any(is.na(sim_stats))) {
    return(Inf)
  }
  sum((sim_stats - obs_stats)^2)
}

# Initial guess and parameter bounds
init <- c(0.005, 1, 1, 1, 0.2, 1, 1, 1)
fit <- optim(
  par = init,
  fn = objective,
  method = "L-BFGS-B",
  lower = c(0, 0, 0, 0, 0, 0, 0, 0),
  upper = c(1, 5, 5, 5, 1, 5, 5, 5)
)

cat("True parameters:\n")
print(true_params)
cat("\nEstimated parameters:\n")
print(fit$par)
