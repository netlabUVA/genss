#' Generate linear rinit C code
#'
#' @param process_args Named list of process arguments
#'
#' @return A Csnippet
#' @export
#'
#' @examples
gen_rinit_lin <- function(process_args) {
  n_states <- process_args[["n_states"]]

  rinit_text <- c()

  for (state_i in seq(n_states)) {
    rinit_text <- c(rinit_text, sprintf("state%d = rnorm(0,1);", state_i))
  }

  return(pomp::Csnippet(rinit_text))
}

#' Generate linear rprocess C code
#'
#' @param process_args Named list of process arguments
#'
#' @return A Csnippet
#' @export
#'
#' @examples
gen_rprocess_lin <- function(process_args) {
  n_states <- process_args[["n_states"]]
  dt <- process_args[["dt"]]

  text <- tempfile(pattern = "rprocess", fileext = ".txt")
  con <- file(text, "w")

  if (is.null(dt)) {
    for (row_i in seq(n_states)) {
      line <- sprintf("double tstate%d = rnorm(", row_i)

      for (col_i in seq(n_states)) {
        if (col_i < n_states) {
          line <- c(line, sprintf("a%d_%d*state%d + ", row_i, col_i, col_i))
        } else {
          line <- c(line, sprintf("a%d_%d*state%d, 1.0);", row_i, col_i, col_i))
        }
      }
      writeLines(line, con)
    }
    writeLines("\n", con)

    for (state_i in seq(n_states)) {
      writeLines(sprintf("state%d = tstate%d;", state_i, state_i), con)
    }
  } else {
    for (row_i in seq(n_states)) {
      line <- sprintf("double dstate%d = rnorm(", row_i)

      for (col_i in seq(n_states)) {
        if (col_i < n_states) {
          line <- c(line, sprintf("a%d_%d*state%d + ", row_i, col_i, col_i))
        } else {
          line <- c(line, sprintf("a%d_%d*state%d, sqrt(dt));", row_i, col_i, col_i))
        }
      }
      writeLines(line, con)
    }
    writeLines("\n", con)

    for (state_i in seq(n_states)) {
      writeLines(sprintf("state%d = dstate%d*dt;", state_i, state_i), con)
    }
  }

  close(con)
  text_str <- paste(readLines(text), collapse = "\n")
  return(pomp::Csnippet(text_str))
}

#' Transform linear process parameters
#'
#' @param process_args Named list of process arguments
#'
#' @return A named list containing the parameters
#' @export
#'
#' @examples
get_transprocess_lin <- function(process_args) {
  process_mat <- process_args[["process_mat"]]
  freed_params <- process_args[["freed_params"]]

  paramnames <- c()
  params <- c()

  # a
  for (row_i in seq(NROW(process_mat))) {
    for (col_i in seq(NCOL(process_mat))) {
      if (is.null(freed_params) || is.na(freed_params[row_i, col_i])) {
        curr_a <- sprintf("a%d_%d", row_i, col_i)
        paramnames <- c(paramnames, curr_a)

        params <- c(params, process_mat[row_i, col_i])
      }
    }
  }

  return(stats::setNames(params, paramnames))
}

#' Validate linear process arguments
#'
#' @param process_args Named list of process arguments
#'
#' @return NULL
#' @export
#'
#' @examples
validate_process_args_lin <- function(process_args) {
  process_mat <- process_args[["process_mat"]]
  freed_params <- process_args[["freed_params"]]
  n_states <- process_args[["n_states"]]
  dt <- process_args[["dt"]]

  if (NROW(process_mat) != NCOL(process_mat)) {
    stop("Process matrix must be square.")
  }

  # if freed params is not null, dim must match process_mat
  if (!is.null(freed_params) && dim(freed_params) != dim(process_mat)) {
    stop("Freed parameters matrix must match the dimension of the process matrix.")
  }
}
