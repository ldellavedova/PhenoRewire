# Run DiffCorr (Fukushima 2013, Gene) on one simulated dataset and return the
# feature pairs it calls as differentially correlated.
#
# Usage: Rscript diffcorr_adapter.R <in.json> <out.json> <alpha>
#
# in.json carries a and b as features x samples. out.json receives a list of
# zero-based [i, j] pairs. An empty list means DiffCorr ran and called nothing;
# a non-zero exit means it could not run.
#
# comp.2.cc.fdr writes a headerless 11-column table. Positions: 1-2 molecule
# names, 3 r in state 1, 4 p1, 5 r in state 2, 6 p2, 7 p for the difference,
# 8 r1 - r2, 9-10 adjusted p1/p2, 11 adjusted p for the difference. Column 11
# is thresholded, so DiffCorr is compared with PhenoRewire at a matched
# false-discovery level. Spearman correlation, BH adjustment.

args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]; outp <- args[2]; alpha <- as.numeric(args[3])

suppressMessages(library(jsonlite))
suppressMessages(library(DiffCorr))

d <- fromJSON(inp)
A <- d$a; B <- d$b
if (is.list(A)) A <- do.call(rbind, A)
if (is.list(B)) B <- do.call(rbind, B)
rownames(A) <- rownames(B) <- paste0("f", seq_len(nrow(A)))

tmp <- tempfile(fileext = ".txt")
invisible(capture.output(
  comp.2.cc.fdr(output.file = tmp, data1 = A, data2 = B,
                method = "spearman", p.adjust.methods = "BH",
                threshold = 1.0, save = TRUE)))
tb <- read.delim(tmp, header = FALSE)
if (ncol(tb) < 11) stop("unexpected DiffCorr output: ", ncol(tb), " columns")
keep <- tb[!is.na(tb[[11]]) & tb[[11]] <= alpha, c(1, 2), drop = FALSE]

if (nrow(keep) == 0) {
  write_json(list(), outp)
} else {
  i <- as.integer(sub("^f", "", keep[[1]])) - 1L
  j <- as.integer(sub("^f", "", keep[[2]])) - 1L
  write_json(cbind(pmin(i, j), pmax(i, j)), outp)
}
