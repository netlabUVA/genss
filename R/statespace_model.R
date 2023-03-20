#' An S4 class for fitting and validating state-space models
#'
#' @slot pomp_obj pomp.
#' @slot data data.frame.
#' @slot measure_type character. List of named parameters for the r and d measure functions.
#' @slot process_type character. List of named parameters for the r and d process functions.
#' @slot measure_args list .
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

#' Generate a state-space model object
#'
#' possible measure args:
#' [grm]
#' - loadings: a (n_measure x n_state) matrix of state-measure loadings. 1 if on, 0 if not. No cross-loadings.
#' - alphas: a vector defining discrimination values for each measure.
#' - betas: a list of vectors defining category boundaries for each measure.
#'
#' [lin]
#' - loadings: a (n_measure x n_state) matrix of state-measure loadings. Cross-loadings allowed.
#' - sigmas: a vector defining standard deviations for each measure.
#'
#' possible process args:
#' [discrete]
#' - process_mat: a (n_state x n_state) transition matrix for each state.
#' - <optional> freed_params: a (n_state x n_state) matrix containing NA where transition values are not free to be estimated.
#' [con]
#' - process_mat: a (n_state x n_state) transition matrix for each state.
#' - <optional> freed_params: a (n_state x n_state) matrix containing NA where transition values are not free to be estimated.
#'
#' @param data The provided dataset. Defaults to NULL.
#' @param n_states The number of states.
#' @param n_measures The number of measures.
#' @param measure_type Measure type selected ("grm" or "lin").
#' @param process_type Process type selected. Defaults to "lin".
#' @param measure_args List of named parameters for the r and d measure functions.
#' @param process_args List of named parameters for the r and d process functions.
#' @param dt Step size for Euler's method. If NULL, discrete time. Defaults to NULL.
#' @param times Vector sequence of times.
#' @param n_particles Number of particles used to estimate parameters. Defaults to 500.
#' @param n_mif2_iter Number of iterations on a mif2 run. Defaults to 1000.
#' @param sim Whether or not a model should be fit or data should be simulated. If FALSE, returns a pomp object. Defaults to FALSE.
#'
#' @return An S4 statespace_mod class or a pomp object.
#' @export
#'
#' @examples
statespace_model <- function(data = NULL,
                             n_states,
                             n_measures,
                             measure_type,
                             process_type = "lin",
                             measure_args,
                             process_args,
                             dt = NULL,
                             times = NULL,
                             n_particles = 500,
                             n_mif2_iter = 1000,
                             sim = FALSE) {
  if (sim) {
    if (is.null(times)) {
      stop("Times vector must be specified for simulation mode.")
    }
  } else {
    if (is.null(data)) {
      stop("Data must be specified for estimate mode.")
    }

    if (!is.null(times)) {
      if (!is.null(dt)) {
        warning("Continuous time observation vector (times) not given. Inferring discrete observations from data.")
      }
      times <- seq_len(NROW(data))
    }
  }

  process_args[["dt"]] <- dt
  process_args[["n_states"]] <- n_states
  process_args[["n_measures"]] <- n_measures

  measure_args[["n_states"]] <- n_states
  measure_args[["n_measures"]] <- n_measures

  do.call(paste("validate_measure_args_", measure_type, sep = ""), list(measure_args))
  do.call(paste("validate_process_args_", process_type, sep = ""), list(process_args))

  obsnames <- sprintf("measure%d", 1:n_measures)
  statenames <- sprintf("state%d", 1:n_states)

  params <- do.call(paste("get_transmeasure_", measure_type, sep = ""), list(measure_args))
  params <- c(params, do.call(paste("get_transprocess_", process_type, sep = ""), list(process_args)))
  paramnames <- names(params)

  rmeas <- do.call(paste("gen_rmeasure_", measure_type, sep = ""), list(measure_args))
  dmeas <- do.call(paste("gen_dmeasure_", measure_type, sep = ""), list(measure_args))
  rproc <- do.call(paste("gen_rprocess_", process_type, sep = ""), list(process_args))
  rinit <- do.call(paste("gen_rinit_", process_type, sep = ""), list(process_args))

  if (is.null(dt)) {
    rproc <- pomp::discrete_time(step.fun = rproc)
  } else {
    rproc <- pomp::euler(step.fun = rproc, delta.t = dt)
  }

  if (sim) {
    if (is.null(data)) {
      data <- data.frame(matrix(0, nrow = length(times), ncol = n_measures))
    }
    pomp::pomp(
      data = data, times = times, t0 = 0, rmeasure = rmeas, dmeasure = dmeas,
      rprocess = rproc, rinit = rinit, obsnames = obsnames, statenames = statenames,
      paramnames = paramnames, params = params
    ) -> pomp_obj

    return(pomp_obj)
  } else {

    # transform alpha within model if grm
    if (measure_type == "grm") {
      n_betas <- sum(lengths(measure_args[["betas"]]))
      alphas <- paramnames[(n_betas + 1):(n_betas + length(measure_args[["alphas"]]))]
      # untransform params
      params[(n_betas + 1):(n_betas + length(measure_args[["alphas"]]))] <- exp(params[(n_betas + 1):(n_betas + length(measure_args[["alphas"]]))])

      pomp::pomp(
        data = data, times = times, t0 = 0, rmeasure = rmeas, dmeasure = dmeas,
        rprocess = rproc, rinit = rinit, obsnames = obsnames, statenames = statenames,
        paramnames = paramnames, params = params, partrans = parameter_trans(log = alphas)
      ) -> pomp_obj
    } else {
      pomp::pomp(
        data = data, times = times, t0 = 0, rmeasure = rmeas, dmeasure = dmeas,
        rprocess = rproc, rinit = rinit, obsnames = obsnames, statenames = statenames,
        paramnames = paramnames, params = params
      ) -> pomp_obj
    }



    return(methods::new("statespace_mod",
      pomp_obj = pomp_obj,
      data = data,
      measure_type = measure_type,
      process_type = process_type,
      measure_args = measure_args,
      process_args = process_args,
      time_args = list("dt" = dt, "times" = times),
      mif2_args = list("n_particles" = n_particles, "n_mif2_iter" = n_mif2_iter),
      param_traces = list(),
      filtered_states = list(),
      mif2_diagnostics = data.frame()
    ))
  }
}
