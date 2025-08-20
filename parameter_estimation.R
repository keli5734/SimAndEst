###############################################################################
## 1.  Load raw synthetic data and summarise per individual                  ##
###############################################################################
library(data.table)
library(dplyr)


source("generate_synthetic_data_standarlised.R")
source("summary_synthetic_data.R")

set.seed(124412135)
n_households <- 50  # increased sample size to reduce estimator bias

raw_dt <- pbapply::pblapply(
  1:n_households,
  sim.hh.func.fixed,
  tests.per.week = 1,  # once per week baseline testing
  
  p.comm.base.infant.fix = 7.148217e-05,
  p.comm.multiplier.sibling = 4.331956e+00,
  p.comm.multiplier.parent = 1.835466e+00,
  p.comm.multiplier.elder = 2,
  
  p.hh.base.infant = 2.888953e-01,
  p.hh.multiplier.sibling = 5.267686e-01,
  p.hh.multiplier.parent = 8.008933e-01,
  p.hh.multiplier.elder = 6.008933e-01,
  
  p.imm.base.sibling = 1e-10,
  p.imm.base.parent = 1e-10,
  p.imm.base.elder = 1e-10,
  
  partial.immunity.infant = 1e-10,
  partial.immunity.sibling = 1e-10,
  partial.immunity.parent = 1e-10,
  partial.immunity.elder = 1e-10,
  
  delay = TRUE,
  duration.latent = 2,
  duration.infect.inf = 6,
  multiplier.dur.sibpar = 0.5,
  p.detect = 0.999,
  amplitude = 2.6581*0,
  phase = -0.408,
  start_date = as.Date("2024-09-21"),         
  end_date   = as.Date("2025-04-17")
)

raw_dt <- rbindlist(raw_dt) 

dt <- summarize_synthetic_data(raw_dt)  %>% arrange(HH)

dt <- as.data.table(dt)

setnames(dt, "individual_ID", "indiv.index")

###############################################################################
## 2.  Define a common study calendar                                        ##
###############################################################################
study_start <- as.Date("2024-09-21")            # same as before
study_end   <- as.Date("2025-04-17")
 
#   If your relative days are counted from *day 0 = study_start*
#   you can convert them directly:
conv <- function(x) fifelse(is.na(x), as.Date(NA), study_start + as.integer(x))

dt[, T_FP_date    := conv( infection.detected.start)]   # first pos PCR
dt[, T_LP_date    := conv( infection.detected.end)]     # last  pos PCR
dt[, last_neg_date:= conv(last_negative)]               # last  neg PCR
dt[, inf_date     := as.Date(NA)]                       # to be imputed
dt[, inf_start_date := T_FP_date]                       # will refine
dt[, inf_end_date   := T_LP_date]

## – If your relative days use *another* origin, adjust `study_start`
##   or subtract the correct offset before calling `conv()`.

## ---------------------------------------------------------------------------
## 3.  Recover the infectious window from the string list                     ##
## ---------------------------------------------------------------------------
parse_vec <- function(s){
  if (is.na(s)) return(NA_integer_)
  as.integer(strsplit(s, ",")[[1]])
}

dt[, `:=`( inf_win_start = sapply( infection.infectious.day,
                                  \(s) min(parse_vec(s), na.rm = TRUE)),
           inf_win_end   = sapply( infection.infectious.day,
                                  \(s) max(parse_vec(s), na.rm = TRUE)) )]

dt[is.finite(inf_win_start), 
   inf_start_date := study_start + inf_win_start]
dt[is.finite(inf_win_end),   
   inf_end_date   := study_start + inf_win_end]

## ---------------------------------------------------------------------------
## 4.  Household‑level variables                                              ##
## ---------------------------------------------------------------------------
dt[, infected := !is.na(T_FP_date)]
dt[, is_index := FALSE]

## earliest detected infection in each HH = index
dt[infected == TRUE,
   is_index := T_FP_date == min(T_FP_date, na.rm = TRUE),
   by = HH]

## replace role with static age bands (infant, toddler, …) if needed
dt[, age_cat :=
     fifelse(role == "infant",  1L,
             fifelse(role == "sibling", 2L,       # all “siblings” treated alike
                     fifelse(role == "adult",   3L,
                             fifelse(role == "elder",   4L, NA_integer_))))]

if (any(is.na(dt$age_cat)))
  stop("Unknown role labels remain – please extend the mapping.")

## ---------------------------------------------------------------------------
## 5.  Relative‑day indices required by the likelihood                        ##
## ---------------------------------------------------------------------------
dt[, `:=`( inf_day_rl            = as.integer(inf_date       - study_start),
           infectious_day_rl     = as.integer(inf_start_date - study_start),
           infectious_end_day_rl = as.integer(inf_end_date   - study_start))]

## non‑infected individuals: leave as NA
dt[!infected == T, c("inf_day_rl","infectious_day_rl",
                     "infectious_end_day_rl") := NA]

## ---------------------------------------------------------------------------
## 6.  Prepare columns matching your earlier pipeline                         ##
## ---------------------------------------------------------------------------
setnames(dt, "HH", "ID_hh")                         # same identifier names
dt[, ID_indiv := sprintf("HH%03d_%02d", ID_hh, indiv.index)]



## 3·3  Infection flags & index
dt[, infected := !is.na(T_FP_date)]
dt[, is_index := FALSE]
dt[infected == TRUE,
   is_index := T_FP_date == min(T_FP_date, na.rm = TRUE),
   by = ID_hh]

## 3·4  Observation window (full study for everybody)
dt[, obs_start_date := study_start]
dt[, obs_end_date   := study_end]



###############################################################################
## 4.  Relative‑day indices
###############################################################################
dt[, `:=`( inf_day_rl            = as.integer(inf_date       - study_start),
           infectious_day_rl     = as.integer(inf_start_date - study_start),
           infectious_end_day_rl = as.integer(inf_end_date   - study_start))]




###############################################################################
## 5.  Delay distribution parameters  (same for imputation & simulation)
###############################################################################
latent_par  <- list(shape = 2, scale = 1)
report_par  <- list(shape = 1, scale = 1.5)
infect_par  <- list(shape = 3, scale = 2)

rtrunc_gamma <- function(n, shape, scale, upper){
  upper[upper<=0] <- 1e-8
  qgamma(runif(n)*pgamma(upper, shape, rate=1/scale),
         shape = shape, rate = 1/scale)
}

###############################################################################
## 6.  Likelihood with 4 age strata         (δ0 α0 γ2 γ3 γ4 β2 β3 β4)
###############################################################################
log1mexp <- function(x) ifelse(x>-0.693, log(-expm1(x)), log1p(-exp(x)))

negll <- function(par, dat, lambda = 0.01, eps = 1e-10){


  delta0 <- par[1]

  alpha0 <- par[2]
  gamma  <- c(0, par[3:5])      # γ2 γ3 γ4     length = 4
  beta   <- c(0, par[6:8])      # β2 β3 β4

  age <- dat$agegrp
  p_comm <- exp(delta0 + gamma[age])
  p_comm <- pmin(pmax(p_comm, eps), 1 - eps)

  p_hh <- exp(alpha0 + beta)    # per-contact infection prob by source age
  p_hh <- pmin(pmax(p_hh, eps), 1 - eps)

  q <- (1 - p_comm) *
       (1 - p_hh[1])^dat$n_inf_infant *
       (1 - p_hh[2])^dat$n_inf_sibling *
       (1 - p_hh[3])^dat$n_inf_parent *
       (1 - p_hh[4])^dat$n_inf_elder
  p_tot <- pmin(pmax(1 - q, eps), 1 - eps)

  -sum(dat$event * log(p_tot) +
         (1 - dat$event) * log1mexp(log(p_tot))) +
    lambda * (sum(gamma[-1]^2) + sum(beta[-1]^2))
}

###############################################################################
## 7.  Repeated imputation + ML
###############################################################################
n_runs <- 1                               # number of repetitions
theta_mat <- matrix(NA_real_, n_runs, 8)
vcov_list <- vector("list", n_runs)
tmax = -as.integer(as.Date("2024-09-21") - as.Date("2025-04-17"))

start_par <- c(0, -2, rep(0, 6))       # length 8

multi_start_optim <- function(start_par, fn, dat, n_start = 5, ...) {
  best_fit <- NULL
  best_val <- Inf
  for (i in seq_len(n_start)) {
    trial <- start_par  + rnorm(length(start_par), 0, 1)
    fit <- optim(trial, fn = fn, dat = dat, ...)
    if (fit$value < best_val) {
      best_val <- fit$value
      best_fit <- fit
    }
  }
  best_fit
}

dt[, latent_delay := NA_real_]   # or NA if you want a logical column
dt[, report_delay := NA_real_]   # or NA if you want a logical column
dt[, infect_period := NA_real_]   # or NA if you want a logical column



for (m in 1:n_runs){
  imp <- copy(dt)
  
  ## --- draw delays for infected individuals ---------------------------
  idx <- which(imp$infected)
  if (length(idx)){
    max_back <- ifelse(is.na(imp$last_neg_date[idx]),
                       14,
                       as.integer(imp$T_FP_date[idx] - imp$last_neg_date[idx]))
    lat  <- rtrunc_gamma(length(idx), latent_par$shape, latent_par$scale, max_back)
    repd <- rtrunc_gamma(length(idx), report_par$shape, report_par$scale,
                         pmax(max_back - lat, 1e-8))
    infd <- rgamma(length(idx), infect_par$shape, rate = 1/infect_par$scale)
    
    imp$latent_delay[idx] <- lat
    imp$report_delay[idx]   <- repd
    imp$infect_period[idx]  <- infd
    
    imp$inf_date[idx] <- imp$T_FP_date[idx] - (lat + repd)
    has_neg <- !is.na(imp$last_neg_date[idx])
    imp$inf_date[idx][has_neg] <- pmax(imp$inf_date[idx][has_neg],
                                       imp$last_neg_date[idx][has_neg] + 1)
    
    imp$inf_start_date[idx] <- imp$inf_date[idx] + lat
    imp$inf_end_date[idx]   <- pmin(imp$inf_start_date[idx] + infd,
                                    imp$T_LP_date[idx])
    
    imp$inf_day_rl[idx]            <- as.integer(imp$inf_date[idx]       - study_start)
    imp$infectious_day_rl[idx]     <- as.integer(imp$inf_start_date[idx] - study_start)
    imp$infectious_end_day_rl[idx] <- as.integer(imp$inf_end_date[idx]   - study_start)
  }
  
  ## --- build person‑day table -----------------------------------------
  rows <- list()
  for (hh in unique(imp$ID_hh)){
    hhdat <- imp[ID_hh == hh]
    n_inf <- matrix(0L, nrow = tmax + 1, ncol = 4)
    for (j in hhdat[infected==TRUE]$ID_indiv){
      r <- hhdat[ID_indiv == j]
      a <- max(r$infectious_day_rl, 0, na.rm = TRUE)
      b <- min(r$infectious_end_day_rl, tmax, na.rm = TRUE)
      age_idx <- r$age_cat[1]
      if (!is.na(a) && a <= b)
        n_inf[a:b+1L, age_idx] <- n_inf[a:b+1L, age_idx] + 1L
    }
     ## include all individuals (including the index case) so that
    ## community infections contribute information to the likelihood
    for (i in hhdat$ID_indiv){
 
    for (i in hhdat[is_index == FALSE]$ID_indiv){
       rec <- hhdat[ID_indiv == i]; inf_d <- rec$inf_day_rl
      for (d in 0:tmax){
        if (!is.na(inf_d) && d > inf_d) break
        rows[[length(rows)+1]] <- list(
          agegrp = rec$age_cat,
          n_inf_infant = n_inf[d+1L, 1],
          n_inf_sibling = n_inf[d+1L, 2],
          n_inf_parent = n_inf[d+1L, 3],
          n_inf_elder  = n_inf[d+1L, 4],
          event  = as.integer(!is.na(inf_d) && d == inf_d),
          ID_indiv = i)
      }
    }
  }
  long <- rbindlist(rows)
  
  ## --- optimise --------------------------------------------------------
  fit <- multi_start_optim(start_par, fn = negll, dat = long,
                           n_start = 5, method = "BFGS",
                           hessian = FALSE, control = list(maxit = 2e4))
  
  
  theta_mat[m,] <- fit$par
  #vcov_list[[m]] <- vc
  
}

###########################################################################
## 8.  Simple comparison: mean estimate vs. truth  (no Rubin pooling)    ##
###########################################################################
## (Keep rows whose optimization succeeded)
keep <- complete.cases(theta_mat[,1])
theta_mat <- theta_mat[keep,,drop = FALSE]

mean_est <- colMeans(theta_mat)           # average across runs

par_names <- c("delta0","alpha0",
               "gamma2","gamma3","gamma4",
               "beta2","beta3","beta4")

## ---- TRUE VALUES for the 4‑category model ----
true_vec <- c(
  # Community infection
  delta0 = log(7.148217e-05),                        # ≈ -9.546
  gamma2 = log(7.148217e-05 * 4.331956e+00) - log(7.148217e-05),                # sibling multiplier
  gamma3 = log(7.148217e-05 * 1.835466e+00) - log(7.148217e-05),                # parent multiplier
  gamma4 = log(7.148217e-05 * 2) - log(7.148217e-05),                           # elder multiplier

  # Household transmission
  alpha0 = log(0.2888953),                           # ≈ -1.242
  beta2  = log(0.2888953 * 0.5267686) - log(0.2888953), # sibling / infant
  beta3  = log(0.2888953 * 0.8008933) - log(0.2888953), # parent  / infant
  beta4  = log(0.2888953 * 0.6008933) - log(0.2888953)  # elder   / infant
)
true_vec <- true_vec[par_names]           # align order

## ---- summary table ----
result <- data.table(
  Parameter = par_names,
  True      = round(true_vec, 3),
  Estimate  = round(mean_est, 3),
  Bias      = round(mean_est - true_vec, 3))

cat("\n*** Mean of", nrow(theta_mat), "runs ***\n")
print(result)


#saveRDS(result, "result.rds")
