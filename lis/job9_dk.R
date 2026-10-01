## Glass-Barriers job 9 (Denmark): two modes among transfer-reliant older singles among one-person households. Output: test statistics and scalar summaries only.
suppressPackageStartupMessages({ library(lissyrtools); library(purrr) })
## Base-R cell test for LISSY (no compiled code). Delta = sup_lambda {E2 - E1} on a lambda grid.
dlt <- function(p, lam) { K <- length(p); g0 <- c(0, cumsum(p) - lam * seq_len(K))
  e1 <- max(g0 - cummin(g0)); left <- cummax(g0 - cummin(g0)); rg <- rev(g0); right <- rev(cummax(cummax(rg) - rg))
  max(left + right) - e1 }
Delta <- function(p, ngrid = 1000) { lv <- sort(unique(c(p, seq(min(p), max(p), length.out = ngrid))))
  max(vapply(lv, function(l) dlt(p, l), 0)) }
argvalley <- function(p, ngrid = 1000) { lv <- sort(unique(c(p, seq(min(p), max(p), length.out = ngrid))))
  d <- vapply(lv, function(l) dlt(p, l), 0); lam <- lv[which.max(d)]
  K <- length(p); g <- cumsum(p) - lam * seq_len(K); g0 <- c(0, g)
  left <- cummax(g0 - cummin(g0))[-1]; rg <- rev(g0); right <- rev(cummax(cummax(rg) - rg))[-1]
  t <- which.max(left[-K] + right[-1]); gl <- g[1:t]; j1 <- which.max(gl - cummin(gl))
  gr <- g[(t + 1):K]; i2 <- t + which.max(rev(cummax(rev(gr))) - gr); c(j1, i2) }
pava_inc <- function(y) { n <- length(y); if (!n) return(y); v <- numeric(n); ww <- numeric(n); m <- 0
  for (i in 1:n) { m <- m + 1; v[m] <- y[i]; ww[m] <- 1
    while (m > 1 && v[m - 1] > v[m]) { v[m - 1] <- (ww[m - 1] * v[m - 1] + ww[m] * v[m]) / (ww[m - 1] + ww[m]); ww[m - 1] <- ww[m - 1] + ww[m]; m <- m - 1 } }
  rep(v[1:m], ww[1:m]) }
unimodal_proj <- function(p) { K <- length(p); best <- NULL
  for (m in 1:K) { f <- numeric(K)
    if (m > 1) f[1:(m - 1)] <- pava_inc(p[1:(m - 1)])
    if (m < K) f[(m + 1):K] <- rev(pava_inc(rev(p[(m + 1):K])))
    f[m] <- max(p[m], if (m > 1) f[m - 1] else 0, if (m < K) f[m + 1] else 0)
    sse <- sum((p - f)^2); if (is.null(best) || sse < best$sse) best <- list(f = f, sse = sse) }
  f <- pmax(best$f, 0); f / sum(f) }
pool_adjacent <- function(p, n, kappa) { K <- length(p); blk <- integer(K); b <- 1; blk[1] <- 1
  for (k in 2:K) { if (abs(p[k] - p[k - 1]) > kappa * sqrt((p[k] + p[k - 1]) / n)) b <- b + 1; blk[k] <- b }
  ave(p, blk) }
## cells of width w centred at multiples of w over [qlo, 95th pct]; returns test, valley and the two mode centres
cell_test_base <- function(y, wt, w, B = 199, seed = 1, qlo = .01) {
  set.seed(seed); ne <- sum(wt)^2 / sum(wt^2)
  lo <- as.numeric(quantile(y, qlo)); hi <- as.numeric(quantile(y, .95))
  k <- floor(y / w + 0.5); keep <- y >= lo & y <= hi; kk <- k[keep]; ww <- wt[keep]
  ks <- min(kk):max(kk); p <- as.numeric(tapply(ww, factor(kk, levels = ks), sum)); p[is.na(p)] <- 0; p <- p / sum(p)
  K <- length(p); obs <- Delta(p); f <- unimodal_proj(pool_adjacent(p, ne, sqrt(log(ne))))
  st <- vapply(seq_len(B), function(b) Delta(rmultinom(1, round(ne), f)[, 1] / round(ne)), 0)
  pv <- (1 + sum(st >= obs - 1e-14)) / (B + 1)
  out <- c(K = K, p = pv, delta = obs, lo = NA, hi = NA, pslo = NA, pshi = NA, m1 = NA, m2 = NA)
  if (obs > 0) { v <- argvalley(p); c0 <- ks * w
    below <- mean(y < lo); out["lo"] <- c0[v[1]] - w / 2; out["hi"] <- c0[v[2]] + w / 2
    out["pslo"] <- below + (mean(keep)) * (if (v[1] > 1) sum(p[1:(v[1] - 1)]) else 0)
    out["pshi"] <- below + (mean(keep)) * sum(p[1:v[2]])
    out["m1"] <- c0[which.max(p[1:v[1]])]; out["m2"] <- c0[v[2] - 1 + which.max(p[v[2]:K])] }
  out }
FP <- c(dk04 = 111924, dk07 = 119244, dk10 = 133272, dk13 = 141984, dk16 = 149544, dk19 = 159000)
hv <- c("dhi", "hitotal", "hipubsoc", "hi31", "hi311", "hi32", "hi33", "hi46", "hilabour", "hicapital", "nhhmem", "hwgt")
hh <- lissyuse(data = c("dk"), vars = hv, from = 2004, to = 2019)
pp <- lissyuse(data = c("dk"), vars = c("age", "sex"), from = 2004, to = 2019)
pc <- function(tag, r) tryCatch(cat(tag, paste(names(r), signif(r, 5), sep = "=", collapse = " "), "\n"), error = function(e) cat(tag, "ERROR", conditionMessage(e), "\n"))
z <- function(x) ifelse(is.na(x), 0, x)
qq <- function(x) paste(round(quantile(x, c(.1, .25, .5, .75, .9), na.rm = TRUE)), collapse = "/")
walk(intersect(names(FP), substr(names(hh), 1, 4)), function(nm) {
  h <- hh[[paste0(nm, "h")]]; for (v in setdiff(hv, names(h))) h[[v]] <- NA
  h <- h[!is.na(h$nhhmem) & h$nhhmem == 1 & !is.na(h$dhi) & h$dhi > 0, ]
  p <- pp[[paste0(nm, "p")]]; m <- match(h$hid, p$hid); h$age <- p$age[m]; h$sex <- p$sex[m]
  pub <- z(h$hi31) + z(h$hi32) + z(h$hipubsoc); oth <- z(h$hi33) + z(h$hilabour) + z(h$hicapital)
  wt <- if (all(is.na(h$hwgt))) rep(1, nrow(h)) else h$hwgt; F <- FP[[nm]]
  ro <- h$age >= 65 & !is.na(h$age) & (pub + oth) > 0 & oth <= 0.1 * (pub + oth)
  hb <- ro & z(h$hi46) > 0; nb <- ro & z(h$hi46) <= 0
  cat("HB", nm, "reliant n", sum(ro), "share_with_housing_benefit", round(sum(wt[hb]) / sum(wt[ro]), 3),
      "hi46_q_among_recipients", qq(h$hi46[hb]), "\n")
  cat("HI311", nm, "reliant hi311 q", qq(h$hi311[ro]), "share_hi311_within2pct_of_FP", round(sum(wt[ro & abs(z(h$hi311) - F) <= .02 * F]) / sum(wt[ro]), 3),
      "hi32 q", qq(h$hi32[ro]), "hipubsoc q", qq(h$hipubsoc[ro]), "\n")
  cat("DHI", nm, "withHB", qq(h$dhi[hb]), "noHB", qq(h$dhi[nb]), "\n")
  cat("SEX", nm, "share_women_reliant", round(sum(wt[ro & h$sex %in% 2]) / sum(wt[ro]), 3), "\n")
  if (sum(hb) > 300) pc(paste("CELL", nm, "reliant_withHB trim.05 w5000"), cell_test_base(h$dhi[hb], wt[hb], 5000, qlo = .05))
  if (sum(nb) > 300) pc(paste("CELL", nm, "reliant_noHB trim.05 w5000"), cell_test_base(h$dhi[nb], wt[nb], 5000, qlo = .05))
  g <- h$hitotal; pc(paste("CELL", nm, "reliant_gross trim.05 w5000"), cell_test_base(g[ro & !is.na(g)], wt[ro & !is.na(g)], 5000, qlo = .05))
  x <- z(h$hi311); pc(paste("CELL", nm, "reliant_hi311 trim.05 w5000"), cell_test_base(x[ro & x > 0], wt[ro & x > 0], 5000, qlo = .05))
})
