source("cell_test.R")
flattop <- function(n) { u <- runif(n); ifelse(u < .6, runif(n, 120000, 160000), ifelse(u < .8, 120000 - abs(rnorm(n, 0, 15000)), 160000 + abs(rnorm(n, 0, 25000)))) }
args <- commandArgs(TRUE); n <- as.integer(args[1]); R <- as.integer(args[2]); s0 <- as.integer(args[3])
set.seed(s0); pv <- replicate(R, cell_test_base(flattop(n), rep(1, n), 5000, B = 199, seed = sample.int(1e6, 1), qlo = .05)["p"])
write(pv, file = sprintf("output/pv_flattop_%d_%d.txt", n, s0), ncolumns = 1)
