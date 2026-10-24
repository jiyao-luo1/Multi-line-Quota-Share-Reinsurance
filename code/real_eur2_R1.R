## European data analysis: Case 2

case2_script_directory <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L && nzchar(sub("^--file=", "", file_arg))) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg), winslash = "/", mustWork = FALSE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

CASE2_SCRIPT_DIR <- case2_script_directory()
previous_skip_main <- Sys.getenv(
  "REINSURANCE_SKIP_MAIN", unset = NA_character_
)
Sys.setenv(REINSURANCE_SKIP_MAIN = "1")
source(file.path(CASE2_SCRIPT_DIR, "real_eur1_R1.R"))

if (is.na(previous_skip_main)) {
  Sys.unsetenv("REINSURANCE_SKIP_MAIN")
} else {
  Sys.setenv(REINSURANCE_SKIP_MAIN = previous_skip_main)
}

case2_main <- function() {
  B <- as.integer(Sys.getenv("REINSURANCE_B", unset = "1000"))
  ncore <- as.integer(Sys.getenv("REINSURANCE_NCORE", unset = "14"))
  N_sim <- as.integer(Sys.getenv("REINSURANCE_NSIM", unset = "51"))

  detected_cores <- parallel::detectCores(logical = TRUE)
  if (is.finite(detected_cores)) ncore <- min(ncore, detected_cores)

  settings <- list(
    alpha = 0.99,
    beta_vec = c(0.05, 0.03, 0.02),
    tilde_beta = 0.05863,
    tilde_beta_vec = c(0.10, 0.06, 0.04),
    mu0 = -1
  )

  X <- load_loss_matrix(
    empirical_data_file(), seed = 9L
  )
  start_time <- proc.time()[["elapsed"]]
  sim.list.2 <- run_case_analysis(
    X, settings, B = B, ncore = ncore, N_sim = N_sim, seed = 92000L
  )
  elapsed <- proc.time()[["elapsed"]] - start_time
  attr(sim.list.2, "elapsed_seconds") <- elapsed

  save(
    sim.list.2,
    file = file.path(empirical_output_directory(), "real_eur2_R1.Rda")
  )

  diagnostics <- attr(sim.list.2, "diagnostics")
  cat("Case 2 completed in", round(elapsed, 1), "seconds.\n")
  cat(
    "Baseline theta estimates:",
    paste(round(diagnostics$baseline_theta, 6), collapse = ", "), "\n"
  )
  cat(
    "Point-estimate c range:",
    paste(round(diagnostics$point_c_range, 6), collapse = " to "), "\n"
  )
  cat(
    "Bootstrap c range:",
    paste(round(diagnostics$bootstrap_c_range, 6), collapse = " to "), "\n"
  )
  cat(
    "Basic 90% interval range:",
    paste(
      round(diagnostics$basic_90_interval_range, 6), collapse = " to "
    ), "\n"
  )
}

if (is.na(previous_skip_main) || previous_skip_main != "1") {
  case2_main()
}
