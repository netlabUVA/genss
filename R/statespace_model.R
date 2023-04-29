#' An S4 class for fitting and validating state-space models.
#'
#' @slot pomp_obj pomp.
#' @slot data data.frame.
#' @slot measure_type character. List of named parameters for the r and d measure functions.
#' @slot process_type character. List of named parameters for the r and d process functions.
#' @slot measure_args list.
#' @slot process_args list.
#' @slot time_args list.
#' @slot mif2_args list.
#' @slot param_traces list.
#' @slot filtered_states list.
#' @slot mif2_diagnostics data.frame.
#'
#' @export
#'
setClass("statespace_mod", slots = c(
  pomp_obj = "pomp",
  data = "data.frame",
  measure_type = "character",
  process_type = "character",
  measure_args = "list",
  process_args = "list",
  time_args = "list",
  mif2_args = "list",
  param_traces = "list",
  filtered_states = "list",
  mif2_diagnostics = "data.frame"
))


check_general_args <- function(n_measures,
                               process_mat,
                               loadings_mat,
                               times,
                               sim) {
  if (any(dim(process_mat) != c(2, 2))) {
    stop("Incorrect dimension: process matrix should be 2x2.")
  }

  if (any(dim(loadings_mat) != c(n_measures, 2))) {
    stop("Incorrect dimension: loadings matrix should be n_measures x 2.")
  }

  if (sim && is.null(times)) {
    stop("Times vector must be specified for simulation mode.")
  }

  if (!sim && is.null(data)) {
    stop("Data must be specified for estimate mode.")
  }
}


#' Generate a two-state linear measure state-space model object.
#'
#' @param data The provided dataset. Defaults to NULL.
#' @param n_measures The number of measures.
#' @param process_mat A 2 x 2 transition matrix for each state.
#' @param loadings_mat A (n_measure x 2) matrix of state-measure loadings. Cross-loadings allowed.
#' @param sigmas A vector defining standard deviations for each measure.
#' @param times A vector sequence of times.
#' @param n_particles Number of particles used to estimate parameters. Defaults to 500.
#' @param n_mif2_iter Number of iterations on a mif2 run. Defaults to 1000.
#' @param sim Whether or not a model should be fit or data should be simulated. If FALSE, returns a pomp object. Defaults to FALSE.
#'
#' @return An S4 statespace_mod class or a pomp object.
#' @export
#'
#' @examples
lin_two_state_statespacemodel <- function(data = NULL,
                                          n_measures,
                                          process_mat,
                                          loadings_mat,
                                          sigmas,
                                          times = NULL,
                                          n_particles = 500,
                                          n_mif2_iter = 1000,
                                          sim = FALSE) {

  check_general_args(n_measures, process_mat, loadings_mat, times, sim)
  if (length(sigmas) != n_measures) {
    stop("Incorrect dimension: sigmas should have length n_measures.")
  }

  sigma_names <- sprintf("sigma%d", 1:n_measures)
  obs_names <- sprintf("measure%d", 1:n_measures)
  state_names <- c("state1", "state2")

  loading_names <- c()
  for (measure_i in seq(n_measures)) {
    loading_names <- c(loading_names, sprintf("l%d_%d", measure_i, 1:2))
  }
  param_names <- c("a1_1", "a1_2", "a2_1", "a2_2", loading_names, sigma_names)
  params <- c(c(process_mat), c(loadings_mat), sigmas)
  
  names(params) <- param_names

  if (sim) {
    if (is.null(data)) {
      data <- data.frame(matrix(0, nrow = length(times), ncol = n_measures))
    }

    pomp::pomp(
      data = data,
      times = times,
      t0 = 0,
      rmeasure = gen_rmeasure_lin(loadings_mat),
      dmeasure = gen_dmeasure_lin(loadings_mat),
      rprocess = pomp::discrete_time(step.fun = gen_rprocess_lin()),
      rinit = gen_rinit_lin(),
      obsnames = obs_names,
      statenames = state_names,
      paramnames = param_names,
      params = params,
      partrans = parameter_trans(log = sigma_names)
    ) -> pomp_obj

    return(pomp_obj)
  } else {

    pomp::pomp(
      data = data,
      times = times,
      t0 = 0,
      rmeasure = gen_rmeasure_lin(loadings_mat),
      dmeasure = gen_dmeasure_lin(loadings_mat),
      rprocess = pomp::discrete_time(step.fun = gen_rprocess_lin()),
      rinit = gen_rinit_lin(),
      obsnames = obs_names,
      statenames = state_names,
      paramnames = param_names,
      params = params, 
      partrans = parameter_trans(log = sigma_names)
    ) -> pomp_obj

    return(methods::new("statespace_mod",
      pomp_obj = pomp_obj,
      data = data,
      measure_type = "lin",
      process_type = "lin",
      measure_args = list("loadings_mat" = loadings_mat, "sigmas" = sigmas),
      process_args = list(),
      time_args = list("times" = times),
      mif2_args = list("n_particles" = n_particles, "n_mif2_iter" = n_mif2_iter),
      param_traces = list(),
      filtered_states = list(),
      mif2_diagnostics = data.frame()
    ))
  }
}

#' Generate a two-state graded-response measure state-space model object.
#'
#' @param data The provided dataset. Defaults to NULL.
#' @param n_measures The number of measures.
#' @param process_mat A 2 x 2 transition matrix for each state.
#' @param loadings_mat A (n_measure x 2) matrix of state-measure loadings. 1 if on, 0 if not. No cross-loadings.
#' @param alphas A vector defining discrimination values for each measure.
#' @param betas A list of vectors defining category boundaries for each measure.
#' @param times A vector sequence of times.
#' @param n_particles Number of particles used to estimate parameters. Defaults to 500.
#' @param n_mif2_iter Number of iterations on a mif2 run. Defaults to 1000.
#' @param sim Whether or not a model should be fit or data should be simulated. If FALSE, returns a pomp object. Defaults to FALSE.
#'
#' @return An S4 statespace_mod class or a pomp object.
#' @export
#'
#' @examples
grm_two_state_statespacemodel <- function(data = NULL,
                                          n_measures,
                                          process_mat,
                                          loadings_mat,
                                          alphas,
                                          betas,
                                          times = NULL,
                                          n_particles = 500,
                                          n_mif2_iter = 1000,
                                          sim = FALSE) {
  check_general_args(n_measures, process_mat, loadings_mat, times, sim)
  if (length(alphas) != n_measures) {
    stop("Incorrect dimension: alphas should have length n_measures.")
  }
  if (length(betas) != n_measures) {
    stop("Incorrect dimension: betas should have length n_measures.")
  }
  if ((TRUE %in% (rowSums(loadings_mat) > 1)) || TRUE %in% (rowSums(loadings_mat) == 0)) {
    stop("A measure must depend on exactly one state.")
  }

  alpha_names <- sprintf("log_alpha%d", 1:n_measures)
  loading_names <- c()
  for (measure_i in seq(n_measures)) {
    loading_names <- c(loading_names, sprintf("l%d_%d", measure_i, 1:2))
  }

  beta_names <- c()
  log_betas <- betas
  for (measure_i in seq_along(betas)){
    for (beta_i in seq_along(betas[[measure_i]])){
      beta_names <- c(beta_names, sprintf("log_beta%d_%d", measure_i, beta_i))
      
      if (beta_i > 1) {
        log_betas[[measure_i]][[beta_i]] <- log(betas[[measure_i]][[beta_i]] - betas[[measure_i]][[beta_i - 1]])
      }
    }
  }

  obs_names <- sprintf("measure%d", 1:n_measures)
  state_names <- c("state1", "state2")

  loading_names <- c()
  for (measure_i in seq(n_measures)) {
    loading_names <- c(loading_names, sprintf("l%d_%d", measure_i, 1:2))
  }
  param_names <- c("a1_1", "a1_2", "a2_1", "a2_2", loading_names, alpha_names, beta_names)
  params <- c(c(process_mat), c(loadings_mat), log(alphas), unlist(log_betas))
  
  names(params) <- param_names

  if (sim) {
    if (is.null(data)) {
      data <- data.frame(matrix(0, nrow = length(times), ncol = n_measures))
    }
    pomp::pomp(
      data = data,
      times = times,
      t0 = 0,
      rmeasure = gen_rmeasure_grm(loadings_mat, log(alphas), log_betas),
      dmeasure = gen_dmeasure_grm(loadings_mat, log(alphas), log_betas),
      rprocess = pomp::discrete_time(step.fun = gen_rprocess_lin()),
      rinit = gen_rinit_lin(),
      obsnames = obs_names,
      statenames = state_names,
      paramnames = param_names,
      params = params
    ) -> pomp_obj

    return(pomp_obj)
  } else {

    pomp::pomp(
      data = data,
      times = times,
      t0 = 0,
      rmeasure = gen_rmeasure_grm(loadings_mat, log(alphas), log_betas),
      dmeasure = gen_dmeasure_grm(loadings_mat, log(alphas), log_betas),
      rprocess = pomp::discrete_time(step.fun = gen_rprocess_lin()),
      rinit = gen_rinit_lin(),
      obsnames = obs_names,
      statenames = state_names,
      paramnames = param_names,
      params = params
    ) -> pomp_obj

    return(methods::new("statespace_mod",
      pomp_obj = pomp_obj,
      data = data,
      measure_type = "grm",
      process_type = "lin",
      measure_args = list("loadings_mat" = loadings_mat, "alphas" = alphas, "log_alphas" = log(alphas), "betas" = betas, "log_betas" = log_betas),
      process_args = list(),
      time_args = list("times" = times),
      mif2_args = list("n_particles" = n_particles, "n_mif2_iter" = n_mif2_iter),
      param_traces = list(),
      filtered_states = list(),
      mif2_diagnostics = data.frame()
    ))
  }
}