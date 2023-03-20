#' Generate linear rmeasure C code
#'
#' @param measure_args Named list of measure arguments
#'
#' @return A Csnippet
#' @export
#'
#' @examples
gen_rmeasure_lin <- function(measure_args) {
  loadings <- measure_args[["loadings"]]
  sigmas <- measure_args[["sigmas"]]
  rmeasure_text <- c()

  for (measure_i in seq(NROW(loadings))) {
    line <- sprintf("measure%d = rnorm(", measure_i)

    for (state_i in seq(NCOL(loadings))) {
      if (state_i < NCOL(loadings)) {
        line <- c(line, sprintf("l%d_%d*state%d + ", measure_i, state_i, state_i))
      } else {
        line <- c(line, sprintf("l%d_%d*state%d, exp(sigma%d));", measure_i, state_i, state_i, measure_i))
        rmeasure_text <- c(rmeasure_text, line)
      }
    }
  }
  return(pomp::Csnippet(rmeasure_text))
}

#' Generate linear dmeasure C code
#'
#' @param measure_args Named list of measure arguments
#'
#' @return A Csnippet
#' @export
#'
#' @examples
gen_dmeasure_lin <- function(measure_args) {
  loadings <- measure_args[["loadings"]]
  dmeasure_text <- c("lik = 0;")

  for (measure_i in seq(NROW(loadings))) {
    line <- sprintf("lik += dnorm(measure%d, ", measure_i)

    for (state_i in seq_along(loadings[measure_i, ])) {
      if (state_i < NCOL(loadings)) {
        line <- c(line, sprintf("l%d_%d*state%d + ", measure_i, state_i, state_i))
      } else {
        line <- c(line, sprintf("l%d_%d*state%d, exp(sigma%d), 1);", measure_i, state_i, state_i, measure_i))
        dmeasure_text <- c(dmeasure_text, line)
      }
    }
  }

  dmeasure_text <- c(dmeasure_text, "lik = (give_log) ? lik : exp(lik);")

  return(pomp::Csnippet(dmeasure_text))
}

#' Transform linear parameters
#'
#' @param measure_args Named list of measure arguments
#'
#' @return A named list containing the parameters
#' @export
#'
#' @examples
get_transmeasure_lin <- function(measure_args) {
  sigmas <- measure_args[["sigmas"]]
  loadings <- measure_args[["loadings"]]

  paramnames <- c()
  params <- c()

  for (sigma_i in seq(length(sigmas))) {
    curr_sigma <- sprintf("sigma%d", sigma_i)
    paramnames <- c(paramnames, curr_sigma)

    params <- c(params, log(sigmas[[sigma_i]]))
  }

  for (measure_i in seq(NROW(loadings))) {
    for (state_i in seq(NCOL(loadings))) {
        curr_loading <- sprintf("l%d_%d", measure_i, state_i)
        paramnames <- c(paramnames, curr_loading)

        params <- c(params, loadings[measure_i, state_i])
    }
  }

  return(stats::setNames(params, paramnames))
}

#' Validate linear measure arguments
#'
#' @param measure_args Named list of measure arguments
#'
#' @return NULL
#' @export
#'
#' @examples
validate_measure_args_lin <- function(measure_args) {
  loadings <- measure_args[["loadings"]]
  n_states <- measure_args[["n_states"]]
  n_measures <- measure_args[["n_measures"]]
  sigmas <- measure_args[["sigmas"]]

  if (NCOL(loadings) != n_states || NROW(loadings) != n_measures) {
    stop("State-measure loadings must have dimensions n_measures x n_states.")
  }

  if (length(sigmas) != n_measures) {
    stop("Measure standard deviation must be defined for each measure.")
  }
}
