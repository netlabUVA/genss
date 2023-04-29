library(pomp)
library(genss)

# this will probably become a vignette at some point

my_data <- diag(1, 10, 3)

# discrete time grm
loadings <- matrix(data = c(c(1, 0), c(1, 0), c(1, 0)), nrow = 3, ncol = 2, byrow = TRUE)
alphas <- c(1, 1, 1)
betas <- list(
  c(-2, -1, 0, 1),
  c(-1.5, -.5, .5, 1.5),
  c(-2.5, -1.5, -.5, .5)
)
process_mat <- matrix(data = c(c(.5, .6), c(.7, .8)), nrow = 2, ncol = 2)
process_mat
c(process_mat)
measure_args <- list("loadings" = loadings, "alphas" = alphas, "betas" = betas)
process_args <- list("process_mat" = process_mat)

grm <- statespace_model(
  measure_type = "grm",
  n_states = NCOL(loadings),
  n_measures = NROW(loadings),
  measure_args = measure_args,
  process_args = process_args,
  times = c(1:10, 20:30),
  sim = TRUE
)
spy(grm)

# plot(simulate(grm))

# continous time grm
process_mat <- matrix(data = c(c(-.5, .5), c(.5, -.5)), nrow = 2, ncol = 2)
process_args <- list("process_mat" = process_mat)

grm_con <- statespace_model(
  measure_type = "grm",
  n_states = NCOL(loadings),
  n_measures = NROW(loadings),
  measure_args = measure_args,
  process_args = process_args,
  times = c(1:10, 20:30),
  dt = .001,
  sim = TRUE
)
# plot(simulate(grm_con))

# discrete time lin
loadings_mat <- matrix(data = c(c(1, 0), c(0, 1), c(1, 0)), nrow = 3, ncol = 2, byrow = TRUE)
process_mat <- matrix(data = c(c(.5, .5), c(.5, .5)), nrow = 2, ncol = 2)
sigmas <- c(5, 1, 1)


test <- lin_two_state_statespacemodel(
                              n_measures=NROW(loadings_mat),
                              process_mat=process_mat,
                              loadings_mat = loadings_mat,
                              sigmas=sigmas,
                              times=1:100,
                              sim=FALSE)