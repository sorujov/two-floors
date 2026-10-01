source("cell_test.R")
gens <- list(
  normal = function(n) rnorm(n, 157000, 16000),
  lognormal = function(n) exp(rnorm(n, log(180000), .3)),
  flattop = function(n) { u <- runif(n); ifelse(u < .6, runif(n, 120000, 160000), ifelse(u < .8, 120000 - abs(rnorm(n, 0, 15000)), 160000 + abs(rnorm(n, 0, 25000)))) },
  skewpeak = function(n) { u <- runif(n); ifelse(u < .5, 130000 + abs(rnorm(n, 0, 4000)) , 130000 + rexp(n, 1/40000)) })
args <- commandArgs(TRUE); g <- args[1]; n <- as.integer(args[2]); R <- as.integer(args[3])
set.seed(n + nchar(g)); pv <- replicate(R, cell_test_base(gens[[g]](n), rep(1, n), 5000, B = 199, seed = sample.int(1e6, 1), qlo = .05)["p"])
cat("SIZE", g, "n", n, "R", R, "rej05", mean(pv <= .05), "rej10", mean(pv <= .10), "\n")
