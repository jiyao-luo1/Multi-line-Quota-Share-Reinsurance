# Supplementary Code S1: cross-coverage dependence only.
# Base R only. No network use, row-level output, or identifier export.
# Usage: Rscript Code/Data_Dependence_Rev8.R private-data.Rda output-directory
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[1L] else "Code/data_eur_new.Rda"
output <- if (length(args) >= 2L) args[2L] else "Empirical_Rev8"
if (!file.exists(input)) stop("Supply the locally held private Rda file.")
env <- new.env(parent = baseenv())
load(input, envir = env)
if (!exists("data_eur", envir = env, inherits = FALSE)) stop("Expected object data_eur.")
dat <- env$data_eur
cols <- c("CLAIM_TPL", "CLAIM_CAR", "CLAIM_UMC")
if (!is.data.frame(dat) || !all(cols %in% names(dat))) stop("Missing loss columns.")
if (!all(vapply(dat[cols], is.numeric, logical(1)))) stop("Loss columns must be numeric.")
X <- as.matrix(dat[cols])
if (nrow(X) < 2L || any(!is.finite(X)) || any(X < 0)) stop("Expected finite nonnegative losses.")
n <- nrow(X)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
pearson <- cor(X)
# cor(method='spearman') is Pearson correlation of average ranks, including zero ties.
spearman <- cor(X, method = "spearman")
pair <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
dep <- data.frame(line_i = pair[, 1L], line_j = pair[, 2L],
                  pearson = pearson[pair], spearman = spearman[pair])
tail_list <- list()
for (q in c(.95, .99)) {
  threshold <- apply(X, 2L, quantile, probs = q, type = 7, names = FALSE)
  exceed <- sweep(X, 2L, threshold, ">")
  counts <- colSums(exceed)
  joint <- crossprod(1L * exceed)
  joint_pair <- joint[pair]
  tail_list[[as.character(q)]] <- data.frame(q = q, line_i = pair[, 1L], line_j = pair[, 2L],
    threshold_i = threshold[pair[, 1L]], threshold_j = threshold[pair[, 2L]],
    exceed_i = counts[pair[, 1L]], exceed_j = counts[pair[, 2L]], joint_count = joint_pair,
    joint_probability = joint_pair / n,
    conditional_j_given_i = joint_pair / counts[pair[, 1L]],
    conditional_i_given_j = joint_pair / counts[pair[, 2L]],
    independence_conditional_j = counts[pair[, 2L]] / n)
  dep[[paste0("chi_", q)]] <- joint_pair / counts[pair[, 1L]]
}
tail <- do.call(rbind, tail_list)
rownames(tail) <- NULL
write.csv(dep, file.path(output, "dependence_summary.csv"), row.names = FALSE)
write.csv(tail, file.path(output, "tail_exceedance_counts.csv"), row.names = FALSE)
fmt <- function(z, digits = 4L) formatC(z, format = "f", digits = digits)
rows <- vapply(seq_len(nrow(dep)), function(k) paste0(
  "$(", dep$line_i[k], ",", dep$line_j[k], ")$ & ",
  paste(fmt(as.numeric(dep[k, 3:6])), collapse = " & "), " \\\\"), character(1L))
writeLines(c("\\begin{table}[!ht]", "\\centering", "\\small",
  "\\caption{\\cnew{Empirical dependence between coverage types, using all observations including zero losses. The finite-threshold conditional exceedance probabilities $\\widehat\\chi_{ij}(q)$ are defined in the text.}}",
  "\\label{tab:dependence}", "{\\color{red}", "\\begin{tabular}{lrrrr}", "\\toprule",
  "Pair $(i,j)$ & Pearson & Spearman & $\\widehat\\chi_{ij}(.95)$ & $\\widehat\\chi_{ij}(.99)$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}}", "\\end{table}"),
  file.path(output, "dependence_table.tex"))
writeLines(c("Supplementary Code S1: dependence summaries only.",
  "Data are private and are not distributed with this code.",
  "All rows, including zeros, enter the pairwise dependence summaries.",
  "Spearman uses average ranks for ties; tail quantiles use R type 7 and strict >.",
  "Chi(q) is a finite-threshold conditional probability, not an asymptotic tail coefficient.",
  "No p-values or independence-based uncertainty intervals are computed for these descriptive summaries.",
  "Only the three claim columns enter the calculations; no individual records are exported.",
  capture.output(sessionInfo())), file.path(output, "computation_notes.txt"))
cat("Saved aggregate summaries and manuscript tables for", n, "records.\n")
