# Run DGCA (McKenzie 2016, BMC Syst Biol) on one dataset and return the feature
# pairs it calls as differentially correlated.
#
# Usage: Rscript dgca_adapter.R <in.json> <out.json> <alpha> [nperm]
#
# in.json carries a and b as features x samples. out.json receives a list of
# zero-based [i, j] pairs. An empty list means DGCA ran and called nothing; a
# non-zero exit means it could not run.
#
# DGCA is the modern successor of DiffCorr: it computes Spearman correlations
# per condition, converts the difference to a z score, and, with nPerm > 0,
# calibrates that z against permutations of the condition labels and reports an
# empirical FDR. We use the empirical FDR when permutations are requested and
# the Benjamini-Hochberg adjusted p otherwise, so DGCA is scored at the same
# false-discovery level as PhenoRewire and DiffCorr.

args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]; outp <- args[2]; alpha <- as.numeric(args[3])
nperm <- if (length(args) >= 4) as.integer(args[4]) else 0L

suppressMessages(library(jsonlite))
suppressMessages(library(DGCA))

d <- fromJSON(inp)
A <- d$a; B <- d$b
if (is.list(A)) A <- do.call(rbind, A)
if (is.list(B)) B <- do.call(rbind, B)
A <- as.matrix(A); B <- as.matrix(B)

mat <- cbind(A, B)
rownames(mat) <- paste0("f", seq_len(nrow(mat)))
colnames(mat) <- paste0("s", seq_len(ncol(mat)))

design <- cbind(state_a = c(rep(1L, ncol(A)), rep(0L, ncol(B))),
                state_b = c(rep(0L, ncol(A)), rep(1L, ncol(B))))
rownames(design) <- colnames(mat)

res <- suppressWarnings(suppressMessages(
  ddcorAll(inputMat = mat, design = design, compare = c("state_a", "state_b"),
           corrType = "spearman", adjust = "BH", nPerm = nperm,
           sigThresh = 1, corSigThresh = 1, heatmapPlot = FALSE)
))

# ddcorAll returns Gene1, Gene2, correlations and either "pValDiff_adj"
# (BH on the asymptotic p) or "empPVals"/"pValDiff_adj" when permuting.
col <- if (nperm > 0 && "empPVals" %in% colnames(res)) "empPVals" else "pValDiff_adj"
if (!col %in% colnames(res)) stop("unexpected DGCA output columns: ", paste(colnames(res), collapse = ","))
p <- suppressWarnings(as.numeric(res[[col]]))
if (nperm > 0 && col == "empPVals") p <- p.adjust(p, method = "BH")
keep <- res[!is.na(p) & p <= alpha, c("Gene1", "Gene2"), drop = FALSE]

if (nrow(keep) == 0) {
  write_json(list(), outp)
} else {
  i <- as.integer(sub("^f", "", keep[[1]])) - 1L
  j <- as.integer(sub("^f", "", keep[[2]])) - 1L
  write_json(cbind(pmin(i, j), pmax(i, j)), outp)
}
