
setwd(dir =  "/Users/ke/Library/CloudStorage/OneDrive-YaleUniversity/Postdoc projects/CN_household_transmisison_study/datasets/")

#################################################################
##  100‑draw multiple‑imputation ML for chain‑binomial model   ##
#################################################################
library(data.table)
library(readxl)
library(tidyverse)
library(lubridate)

# set.seed(123456)

## ------------------------------------------------------------
## -1.  Prepare data                               ##
## ------------------------------------------------------------

# test_data <- readRDS("sample_T_clean_before20250208.rds") # raw data
# sampled_data <- readRDS("dataset_sample.rds")
#
# # corrections due to data entry errors
# test_data[2619,13] <- "Positive"
# test_data[2619,14] <- "33"
#
#
# test_data <- test_data %>%
#   filter(ID_hh %in% unique(sampled_data$ID_hh)) %>%
#   mutate(RSV_positive = ifelse(RSVA_QL == "Negative" & RSVB_QL == "Negative", "Negative", "Positive"))
#
# start_end_dates <- test_data %>%
#   group_by(ID_indiv) %>%
#   summarise(starting_data = min(sample_T_date),
#             end_data = max(sample_T_date))
#
#
# subset_data <- sampled_data %>%
#   select("ID_hh", "ID_indiv", "date_birth", "comorbidity","class_size",
#          "workplace_size", "episode", "pathogen", "ifsecond","T_RN_date",
#          "T_FP_date", "T_LP_date", "ifsym", "onset_date", "resolution_date")
#
# subset_data$starting_date <- start_end_dates$starting_data
# subset_data$end_date <- start_end_dates$end_data
#
# # T_RN_date: sampling date of the most recent negative sample before the first positive sample in this episode
# # T_FP_date: sampling date of first positive sample in this episode
# # T_LP_date: sampling date of last sample in this episode
#
#
#
# # visualize data (can be figure 1)
#
# subset_data %>%
#   ggplot() +
#   geom_segment(aes(x = ID_indiv, xend = ID_indiv,
#                    y = pmin(T_RN_date,T_FP_date, T_LP_date, onset_date, resolution_date, na.rm = T),
#                    yend = pmax(T_RN_date, T_FP_date, T_LP_date, onset_date, resolution_date, na.rm = T)),
#                color = "black", linewidth = 0.8) +
#   geom_point(aes(x = ID_indiv, y = T_RN_date, color = "most recent negative"), size = 2.5, alpha = .7) +
#   geom_point(aes(x = ID_indiv, y = T_FP_date, color = "first positive"), size = 2.5, alpha = .7) +
#   geom_point(aes(x = ID_indiv, y = T_LP_date, color = "last positive"),  size = 2.5, alpha = .7) +
#   geom_point(aes(x = ID_indiv, y = onset_date, color = "symptom onset"), size = 2.5, alpha = .7)+
#   geom_point(aes(x = ID_indiv, y = resolution_date, color = "symptom resolved"), size = 2.5, alpha = 1) +
#   geom_point(data = test_data, aes(x = ID_indiv, y = sample_T_date, group = ID_indiv, color = RSV_positive), size = 1) +
#   scale_color_manual(values = c("most recent negative" = "grey",
#                                 "first positive" = "#cb181d",
#                                 "last positive" = "#fcbba1",
#                                 "symptom onset" = "#006d2c",
#                                 "symptom resolved" = "#c7e9c0",
#                                 "Negative" = "grey",
#                                 "Positive" = "#cb181d")) +
#   theme_bw()+
#   labs(y = "Date", x = "Individual ID", color = "Event Type") + # Label legend
#   coord_flip()

#  saveRDS(subset_data, "data_subset.rds")


## ------------------------------------------------------------
## 0.  Load baseline data                                     ##
## ------------------------------------------------------------
df <- readRDS("data_subset.rds")   %>% filter(ID_hh == "HH019")   # this is clean a clean dataframe
df <- as.data.table(df)

study_start <- as.Date("2024-09-21")
study_end   <- as.Date("2025-04-17")
tmax        <- as.integer(study_end - study_start) + 1

df[, obs_start_date := as.Date(starting_date)]
df[, obs_end_date   := as.Date(end_date)]

df[is.na(obs_start_date), obs_start_date := starting_date]
df[is.na(obs_end_date  ), obs_end_date   := study_end  ]

## dummy surveillance series – replace with real data
case_df <- read_excel("case_data.xlsx") %>%
  mutate(age_clean = str_replace_all(age, "\\s+", ""),  # Remove all white space incl. \n
         age_clean = case_when(
           str_detect(age_clean, "天") ~ as.numeric(str_remove(age_clean, "天")) / 365.25,
           str_detect(age_clean, "月") ~ as.numeric(str_remove(age_clean, "月")) / 12,
           str_detect(age_clean, "岁") ~ as.numeric(str_remove(age_clean, "岁")),
           str_detect(age_clean, "^\\d+/\\d+$") ~ as.numeric(sapply(str_split(age_clean, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
           TRUE ~ as.numeric(age_clean)
         ))
# Suppose your column is named 'my_date'
case_df$Date <- parse_date_time(case_df$Date, 
                              orders = c("Y/m/d", "m/d/y", "m/d/Y"))

# Convert to Date (drop time part)
case_df$Date <- as.Date(case_df$Date)
case_df$test_numeric <- ifelse(case_df$tests == "阳性", 1,
                          ifelse(case_df$tests == "阴性", 0, NA))

case_t <- case_df %>% 
  group_by(Date) %>% 
  summarise(cases = sum(test_numeric, na.rm = T)) %>% 
  ungroup() %>% 
  complete(Date = seq(min(Date), max(Date), by = "day"), fill = list(cases = 0))

cases_t <- case_t$cases  

## Basic covariates -------------------------------------------
df[, age_years := as.numeric(difftime(study_start, date_birth, units = "days"))/365.25]
df[, age_cat  := cut(age_years,
                     breaks = c(-Inf,1,5,18,65,Inf),
                     labels = 1:5, right = FALSE)]
df[, infected := !is.na(T_FP_date)]
 

## ---- REPLACE your current is_index line with this block ----------------
df[, is_index := FALSE]                               # reset

## (1) flag rows explicitly marked as primary
df[!is.na(ifsecond) & ifsecond == FALSE, is_index := TRUE]

## (2) fallback: earliest T_FP_date per household if none flagged
df[, any_primary := any(is_index), by = ID_hh]
df[any_primary == FALSE & infected == TRUE,
   is_earliest := T_FP_date == min(T_FP_date, na.rm = TRUE),
   by = ID_hh]
df[is_earliest == TRUE, is_index := TRUE]
df[, c("any_primary","is_earliest") := NULL]

## numeric age already exists as age_years
## guarantee character IDs for easy reuse
df[, ID_hh    := as.character(ID_hh   )]
df[, ID_indiv := as.character(ID_indiv)]

## ------------------------------------------------------------
## 1.  Delay distributions                                   ##
latent_par  <- list(shape = 2, scale = 1)
report_par  <- list(shape = 1, scale = 1.5)
infect_par  <- list(shape = 3, scale = 2)

rtrunc_gamma <- function(n, shape, scale, upper) {
  # if upper == 0, return 0
  upper[upper <= 0] <- 1e-10
  u <- runif(n)
  qgamma(u * pgamma(upper, shape, rate = 1/scale),
         shape = shape,
         rate  = 1/scale)
}

## ------------------------------------------------------------
## 2.  Storage objects Storage for Monte‑Carlo output.    ##
M <- 100                  # number of Monte‑Carlo imputations
theta_mat <- matrix(NA_real_, M, 11)   # 11 parameters

## ------------------------------------------------------------
## 3.  Helper functions                                       ##
log1mexp <- function(x) ifelse(x > -0.693,
                               log(-expm1(x)), log1p(-exp(x)))

# Directly doing log(1-exp(x)) fails when x is close to 0 because
# exp(x) is almost 1, so the subtraction 1 - exp(x) loses many digits of
# precision (“catastrophic cancellation”).

negll <- function(par, dat, lambda = 0.01, eps = 1e-10) {

  delta0 <- par[1]; delta1 <- par[2]
  alpha0 <- par[3]
  gamma  <- c(0, par[4:7])           # γ2..γ5
  beta   <- c(0, par[8:11])          # β2..β5

  age <- dat$agegrp

  p_comm <- exp(delta0 + gamma[age] + delta1 * dat$cases)
  p_hh   <- exp(alpha0 + beta[age])

  p_comm <- pmin(pmax(p_comm, eps), 1 - eps)
  p_hh   <- pmin(pmax(p_hh,   eps), 1 - eps)

  p_tot  <- 1 - (1 - p_comm) * (1 - p_hh) ^ dat$n_inf
  p_tot  <- pmin(pmax(p_tot,  eps), 1 - eps) # p_tot between 0 and 1

  LL =  -sum(dat$event * log(p_tot) + (1 - dat$event) * log1mexp(log(p_tot)))

  LL + lambda * (sum(gamma[-1]^2) + sum(beta[-1]^2))
}

multi_start_optim <- function(start_par, fn, dat, n_start = 5, ...) {
  best_fit <- NULL
  best_val <- Inf
  for (i in seq_len(n_start)) {
    trial <- start_par + rnorm(length(start_par), 0, 1)
    fit <- optim(trial, fn = fn, dat = dat, ...)
    if (fit$value < best_val) {
      best_val <- fit$value
      best_fit <- fit
    }
  }
  best_fit
}

start_par <- c(-6, 0.02, -2, rep(0, 8))

## ------------------------------------------------------------
## 4.  Monte‑Carlo loop                                       ##
## ------------------------------------------------------------
for (m in 1:M) {
  
  imp <- copy(df)
  
  imp$inf_date = rep(as.Date(NA), nrow(imp))
  imp$inf_start_date = rep(as.Date(NA), nrow(imp))
  imp$inf_end_date = rep(as.Date(NA), nrow(imp))
  
  
  imp$inf_day_rl = rep(NA_integer_, nrow(imp))
  imp$infectious_day_rl = rep(NA_integer_, nrow(imp))
  imp$infectious_end_day_rl = rep(NA_integer_, nrow(imp))
  
  ## ------------------------------------------------------------
  ## 4‑A  Draw delays for every infected person  (with truncation)
  ## ------------------------------------------------------------
  idx <- which(imp$infected)
  if (length(idx)) {
    
    ## --- 1.  compute per‑person back‑bound ---------------------------
    max_back <- ifelse(is.na(imp$T_RN_date[idx]),
                       14,
                       as.integer(imp$T_FP_date[idx] - imp$T_RN_date[idx]))

    ## --- 2.  latent delay truncated at max_back ----------------------
    lat  <- rtrunc_gamma(length(idx), latent_par$shape, latent_par$scale, max_back)
    ## --- 3.  report delay truncated at max_back - latent -------------
    repd <- rtrunc_gamma(length(idx), report_par$shape, report_par$scale,
                         pmax(max_back - lat, 1e-8))
    ## --- 4.  infectious duration (right truncation later) ------------
    infd <- rgamma(length(idx), infect_par$shape, rate = 1 / infect_par$scale)

    ## --- 5.  construct timeline -------------------------------------
    imp$inf_date[idx] <- imp$T_FP_date[idx] - (lat + repd)
    has_neg <- !is.na(imp$T_RN_date[idx])
    imp$inf_date[idx][has_neg] <- pmax(imp$inf_date[idx][has_neg],
                                       imp$T_RN_date[idx][has_neg] + 1)

    imp$inf_start_date[idx] <- imp$inf_date[idx] + lat
    imp$inf_end_date[idx]   <- pmin(imp$inf_start_date[idx] + infd,
                                    imp$T_LP_date[idx])

    ## --- 6.  relative‑day indices -----------------------------------
    imp$inf_day_rl[idx]            <- as.integer(imp$inf_date[idx]       - study_start)
    imp$infectious_day_rl[idx]     <- as.integer(imp$inf_start_date[idx] - study_start)
    imp$infectious_end_day_rl[idx] <- as.integer(imp$inf_end_date[idx]   - study_start)
  }
  
  ## 4‑B  build person‑day table for ALL households ----------
  rows <- list()
  
  for (hh in unique(imp$ID_hh)) {
    
    hhdat <- imp[ID_hh == hh]

    ## --- A. daily count of infectious household members -------
    n_inf <- integer(tmax + 1)
    for (j in hhdat[infected==TRUE]$ID_indiv){
      r <- hhdat[ID_indiv == j]
      a <- max(r$infectious_day_rl, 0, na.rm = TRUE)
      b <- min(r$infectious_end_day_rl, tmax, na.rm = TRUE)
      if (!is.na(a) && a <= b) n_inf[a:b+1L] <- n_inf[a:b+1L] + 1L
    }

    ## --- B. susceptible rows --------------------------------
    for (i in hhdat[is_index == FALSE]$ID_indiv){
      rec <- hhdat[ID_indiv == i]; inf_d <- rec$inf_day_rl
      for (d in 0:tmax){
        if (!is.na(inf_d) && d > inf_d) break
        rows[[length(rows)+1]] <- list(
          agegrp = as.integer(rec$age_cat),
          n_inf  = n_inf[d+1L],
          cases  = cases_t[d+1L],
          event  = as.integer(!is.na(inf_d) && d == inf_d))
      }
    }
  }

  long <- rbindlist(rows)
  long$cases[is.na(long$cases)] <- 0

  ## 4‑C  fit ML --------------------------------------------
  fit <- multi_start_optim(start_par, fn = negll, dat = long,
                           n_start = 5, method = "BFGS",
                           hessian = FALSE, control = list(maxit = 2e4))
  theta_mat[m, ] <- fit$par
}




## ------------------------------------------------------------
## 5.  Simple average across imputations                      ##
## ------------------------------------------------------------

keep <- complete.cases(theta_mat[,1])
theta_mat <- theta_mat[keep,,drop = FALSE]

mean_est <- colMeans(theta_mat)

par_names <- c("delta0","delta1","alpha0",
               paste0("gamma",2:5), paste0("beta",2:5))

result <- data.table(
  Parameter = par_names,
  Estimate  = round(mean_est, 3)
)

cat("\nMean of", nrow(theta_mat), "runs\n")
print(result)
