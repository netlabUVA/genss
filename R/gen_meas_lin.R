gen_rmeasure_lin <- function(loadings) {
  rmeasure_text <- c()

  for (measure_i in seq_len(NROW(loadings))) {
    line <- sprintf("measure%d = rnorm(", measure_i)

    for (state_i in seq_len(NCOL(loadings))) {
      if (state_i < NCOL(loadings)) {
        line <- c(line, sprintf("l%d_%d*state%d + ", measure_i, state_i, state_i))
      } else {
        line <- c(line, sprintf("l%d_%d*state%d, sigma%d);", measure_i, state_i, state_i, measure_i))
        rmeasure_text <- c(rmeasure_text, line)
      }
    }
  }
  return(pomp::Csnippet(rmeasure_text))
}

gen_dmeasure_lin <- function(loadings) {
  dmeasure_text <- c("lik = 0;")

  for (measure_i in seq_len(NROW(loadings))) {
    line <- sprintf("lik += dnorm(measure%d, ", measure_i)

    for (state_i in seq_along(loadings[measure_i, ])) {
      if (state_i < NCOL(loadings)) {
        line <- c(line, sprintf("l%d_%d*state%d + ", measure_i, state_i, state_i))
      } else {
        line <- c(line, sprintf("l%d_%d*state%d, sigma%d, 1);", measure_i, state_i, state_i, measure_i))
        dmeasure_text <- c(dmeasure_text, line)
      }
    }
  }

  dmeasure_text <- c(dmeasure_text, "lik = (give_log) ? lik : exp(lik);")

  return(pomp::Csnippet(dmeasure_text))
}
