#' Summarize infection episodes in synthetic data
#'
#' Given the data frame produced by `sim.hh.func.fixed`, compute per-individual
#' infection summaries using only the available six columns: HH, individual_ID,
#' role, test_date, infection_status and community_risk.
#'
#' The function returns one row per individual with
#'   * HH — household identifier
#'   * individual_ID — individual identifier within the household
#'   * n.true.infection — number of infection episodes in the record
#'   * n.detected.infection — number of detected infection episodes (assumed
#'     equal to true episodes as detection is encoded in `infection_status`)
#'   * infection.detected.start — first test day with a positive infection
#'   * infection.detected.end — last test day of the first infection episode
#'   * infection.true.duration — duration (in days) of the first infection
#'     episode based on test days
#'   * infection.infectious.day — comma separated list of test days in which the
#'     individual was infected for the first episode
#'
#' @param df Data frame returned by `sim.hh.func.fixed`.
#' @return Data frame summarising infection episodes for each individual.
#' @examples

summarize_synthetic_data <- function(df) {
  df <- df[order(df$HH, df$individual_ID, df$test_date), ]
  split_df <- split(df, list(df$HH, df$individual_ID), drop = TRUE)
  
  res_list <- lapply(split_df, function(ind_df) {
    status <- ind_df$infection_status
    dates <- ind_df$test_date
    role <- unique(ind_df$role)
    community_risk <- mean(ind_df$community_risk, na.rm = TRUE)
    r <- rle(status)
    idx_end <- cumsum(r$lengths)
    idx_start <- idx_end - r$lengths + 1
    
    starts <- dates[idx_start[r$values == 1]]
    ends   <- dates[idx_end[r$values == 1]]
    n_inf  <- sum(r$values == 1)
    dur    <- ends - starts + 1
    last_negative <- if (any(ind_df$infection_status == 0 &
                             ind_df$test_date < min(ind_df$test_date[ind_df$infection_status == 1])))
      max(ind_df$test_date[ind_df$infection_status == 0 &
                             ind_df$test_date < min(ind_df$test_date[ind_df$infection_status == 1])]) else NA_integer_
    
    infectious_day <- if (length(starts) > 0)
      paste(starts[1]:ends[1], collapse = ",")
    else
      NA_character_
    
    data.frame(
      HH = ind_df$HH[1],
      individual_ID = ind_df$individual_ID[1],
      n.true.infection = n_inf,
      n.detected.infection = n_inf,
      infection.detected.start = if (length(starts) > 0) starts[1] else NA_integer_,
      infection.detected.end = if (length(ends) > 0) ends[1] else NA_integer_,
      infection.true.duration = if (length(dur) > 0) dur[1] else NA_integer_,
      last_negative = last_negative,
      infection.infectious.day = infectious_day,
      community.risk = community_risk,
      row.names = NULL,
      role = role
    )
  })
  
  do.call(rbind, res_list)
}

# if (identical(environment(), globalenv()) && !interactive()) {
#   source("generate_synthetic_data_standarlised.R")
#   example_df <- sim.hh.func.fixed(1)
#   print(summarize_synthetic_data(example_df))
# }


### example
# source("generate_synthetic_data_standarlised.R")
# N_famlies = 2
# df <- sim.hh.func.fixed(1:N_famlies)
# df <- bind_rows(df)
# summary_df <- summarize_synthetic_data(df)
# print(summary_df)





