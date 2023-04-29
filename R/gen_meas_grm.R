gen_rmeasure_grm <- function(loadings, alphas, betas) {
  text <- tempfile(pattern = "rmeasure", fileext = ".txt")
  con <- file(text, "w")

  # pull out arrays of beta values to iterate
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("const double *log_beta%dt = &log_beta%d_1;", measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # transform alphas
  for (alpha_i in seq_along(alphas)) {
    writeLines(sprintf("double alpha%dt = exp(log_alpha%d);", alpha_i, alpha_i), con)
  }
  writeLines("\n", con)

  # create beta arrays
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("double beta%d[%d] = {log_beta%dt[0]};", measure_i, length(betas[[measure_i]]), measure_i), con)
  }
  writeLines("\n", con)

  # compute beta values
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("for(int i = 1; i < %d; i++){", length(betas[[measure_i]])), con)
    writeLines(sprintf("beta%d[i] = beta%d[i-1] + exp(log_beta%dt[i]);}", measure_i, measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # compute GRM cumulative densities
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("double p%d[%d] = {1.0};", measure_i, length(betas[[measure_i]]) + 2), con)
  }
  writeLines("\n", con)

  for (measure_i in seq_along(betas)) {
    # find the loaded state at the current measure
    for (state_i in seq_along(loadings[measure_i, ])) {
      if (loadings[measure_i, state_i] == 1) {
        loaded_i <- state_i
        break
      }
    }
    writeLines(sprintf("for(int i = 1; i < %d; i++){", length(betas[[measure_i]]) + 1), con)
    writeLines(sprintf("p%d[i] = 1/(1+exp(-alpha%dt * (state%d - beta%d[i-1])));}", measure_i, measure_i, loaded_i, measure_i), con)
  }
  writeLines("\n", con)

  # compute GRM probabilities
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("double p%ddiff[%d] = {0.0};", measure_i, length(betas[[measure_i]]) + 1), con)
  }
  writeLines("\n", con)

  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("for(int i = 0; i < %d; i++){", length(betas[[measure_i]]) + 1), con)
    writeLines(sprintf("p%ddiff[i] = p%d[i] - p%d[i+1];}", measure_i, measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # simulate measures from probabilites
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("int measure%d_temp = 0;", measure_i), con)
    writeLines(sprintf("int perm%d[%d] = {0};", measure_i, length(betas[[measure_i]]) + 1), con)
    writeLines(sprintf("rmultinom(1, &p%ddiff[0], %d, &perm%d[0]);", measure_i, length(betas[[measure_i]]) + 1, measure_i), con)
  }
  writeLines("\n", con)

  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("for(int i = 0; i < %d; i++){", length(betas[[measure_i]]) + 1), con)
    writeLines(sprintf("if(perm%d[i]==1){measure%d_temp = i+1;}}", measure_i, measure_i), con)
  }
  writeLines("\n", con)

  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("measure%d = (double) measure%d_temp;", measure_i, measure_i), con)
  }

  close(con)
  text_str <- paste(readLines(text), collapse = "\n")
  return(pomp::Csnippet(text_str))
}


gen_dmeasure_grm <- function(loadings, alphas, betas) {
  text <- tempfile(pattern = "dmeasure", fileext = ".txt")
  con <- file(text, "w")

  # pull out measures and convert to integers
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("int measure%di = (int) measure%d;", measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # pull out arrays of beta values to iterate
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("const double *log_beta%dt = &log_beta%d_1;", measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # transform alphas
  for (alpha_i in seq_along(alphas)) {
    writeLines(sprintf("double alpha%dt = exp(log_alpha%d);", alpha_i, alpha_i), con)
  }
  writeLines("\n", con)

  # create beta arrays
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("double beta%d[%d] = {log_beta%dt[0]};", measure_i, length(betas[[measure_i]]), measure_i), con)
  }
  writeLines("\n", con)

  # compute beta values
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("for(int i = 1; i < %d; i++){", length(betas[[measure_i]])), con)
    writeLines(sprintf("beta%d[i] = beta%d[i-1] + exp(log_beta%dt[i]);}", measure_i, measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # compute GRM cumulative densities
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("double p%d[%d] = {1.0};", measure_i, length(betas[[measure_i]]) + 2), con)
  }
  writeLines("\n", con)

  for (measure_i in seq_along(betas)) {
    # find the loaded state at the current measure
    for (state_i in seq_along(loadings[measure_i, ])) {
      if (loadings[measure_i, state_i] == 1) {
        loaded_i <- state_i
        break
      }
    }
    writeLines(sprintf("for(int i = 1; i < %d; i++){", length(betas[[measure_i]]) + 1), con)
    writeLines(sprintf("p%d[i] = 1/(1+exp(-alpha%dt * (state%d - beta%d[i-1])));}", measure_i, measure_i, loaded_i, measure_i), con)
  }
  writeLines("\n", con)

  # compute GRM probabilities
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("double p%ddiff[%d] = {0.0};", measure_i, length(betas[[measure_i]]) + 1), con)
  }
  writeLines("\n", con)

  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("for(int i = 0; i < %d; i++){", length(betas[[measure_i]]) + 1), con)
    writeLines(sprintf("p%ddiff[i] = p%d[i] - p%d[i+1];}", measure_i, measure_i, measure_i), con)
  }
  writeLines("\n", con)

  # get probabilities and calculate log likelihood
  writeLines("lik = 0;", con)
  for (measure_i in seq_along(betas)) {
    writeLines(sprintf("lik += log(p%ddiff[measure%di-1]);", measure_i, measure_i), con)
  }
  writeLines("lik = (give_log) ? lik : exp(lik);", con)

  close(con)
  text_str <- paste(readLines(text), collapse = "\n")
  return(pomp::Csnippet(text_str))
}