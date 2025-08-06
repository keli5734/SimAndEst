# Example script to generate 100 households using sim.hh.func.fixed
# Requires packages: pbapply, dplyr

# number of households to simulate
N.families <- 100

# run the household simulations
household_list <- pbapply::pblapply(
  1:N.families,
  sim.hh.func.fixed,
  tests.per.week = 2,  # twice per week baseline testing
  p.comm.base.infant.fix = 7.148217e-05,
  p.comm.multiplier.sibling = 4.331956e+00,
  p.comm.multiplier.parent = 1.835466e+00,
  p.daycare.infant = 0,
  dayc.age.mean = 3.9,
  dayc.age.sd = 0.8,
  p.dayc.multiplier.infant = 1.5,
  p.hh.base.infant = 2.888953e-01,
  p.hh.multiplier.sibling = 5.267686e-01,
  p.hh.multiplier.parent = 8.008933e-01,
  p.imm.base.sibling = 0.5,
  p.imm.base.parent = 0.8,
  partial.immunity.infant = 0.8,
  partial.immunity.sibling = 0.6,
  partial.immunity.parent = 0.4,
  delay = TRUE,
  duration.latent = 4.5,
  duration.infect.inf = 9,
  multiplier.dur.sibpar = 0.5,
  p.detect = 0.999,
  amplitude = 2.6581,
  phase = -0.408
)

# combine the list of household data frames into one
household_df <- dplyr::bind_rows(household_list)

# preview the first few rows
print(head(household_df))
