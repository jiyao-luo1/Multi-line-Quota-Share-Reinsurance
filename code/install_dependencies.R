# Run once from the repository root: Rscript code/install_dependencies.R
packages <- c("snowfall", "rlecuyer", "copula", "nleqslv", "mvtnorm", "QRM")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
print(data.frame(package = packages,
                 version = vapply(packages, function(p) as.character(packageVersion(p)),
                                  character(1))))
