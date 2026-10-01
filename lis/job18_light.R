## Glass-Barriers job 18 (Denmark, Norway): light power, grid shift, reweighting, low-income lines among one-person households. Output: test statistics and scalar summaries only.
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
cell_test_base <- function(y, wt, w, B = 199, seed = 1, qlo = .01, off = 0) {
  set.seed(seed); ne <- sum(wt)^2 / sum(wt^2)
  lo <- as.numeric(quantile(y, qlo)); hi <- as.numeric(quantile(y, .95))
  k <- floor((y - off) / w + 0.5); keep <- y >= lo & y <= hi; kk <- k[keep]; ww <- wt[keep]
  ks <- min(kk):max(kk); p <- as.numeric(tapply(ww, factor(kk, levels = ks), sum)); p[is.na(p)] <- 0; p <- p / sum(p)
  K <- length(p); obs <- Delta(p); f <- unimodal_proj(pool_adjacent(p, ne, sqrt(log(ne))))
  st <- vapply(seq_len(B), function(b) Delta(rmultinom(1, round(ne), f)[, 1] / round(ne)), 0)
  pv <- (1 + sum(st >= obs - 1e-14)) / (B + 1)
  out <- c(K = K, p = pv, delta = obs, lo = NA, hi = NA, pslo = NA, pshi = NA, m1 = NA, m2 = NA)
  if (obs > 0) { v <- argvalley(p); c0 <- ks * w + off
    below <- mean(y < lo); out["lo"] <- c0[v[1]] + w / 2; out["hi"] <- c0[v[2]] + w / 2
    out["pslo"] <- below + (mean(keep)) * sum(p[1:v[1]])
    out["pshi"] <- below + (mean(keep)) * sum(p[1:v[2]])
    out["m1"] <- c0[which.max(p[1:v[1]])]; out["m2"] <- c0[v[2] - 1 + which.max(p[v[2]:K])] }
  out }
hv <- c("dhi", "hipubsoc", "hi31", "hi32", "hi33", "hi46", "hilabour", "hicapital", "nhhmem", "hwgt")
pc <- function(tag, r) tryCatch(cat(tag, paste(names(r), signif(r, 5), sep = "=", collapse = " "), "\n"), error = function(e) cat(tag, "ERROR", conditionMessage(e), "\n"))
z <- function(x) ifelse(is.na(x), 0, x)
wq <- function(x, w, q) { o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w) / sum(w); x[which(cw >= q)[1]] }
for (nm in c("dk13", "dk19", "no19")) { cc <- substr(nm, 1, 2); yr <- 2000 + as.integer(substr(nm, 3, 4)); a0 <- if (cc == "dk") 65 else 67
  hh <- lissyuse(data = cc, vars = hv, from = yr, to = yr); pp <- lissyuse(data = cc, vars = c("age"), from = yr, to = yr)
  H <- hh[[paste0(nm, "h")]]; p <- pp[[paste0(nm, "p")]]; rm(hh, pp)
  H <- H[!is.na(H$dhi) & !is.na(H$nhhmem) & H$nhhmem > 0, ]
  Wp <- (if (all(is.na(H$hwgt))) 1 else H$hwgt) * H$nhhmem; med <- wq(H$dhi / sqrt(H$nhhmem), Wp, .5); L50 <- .5 * med; L60 <- .6 * med
  h <- H[H$nhhmem == 1 & H$dhi > 0, ]; h$age <- p$age[match(h$hid, p$hid)]; rm(p, H)
  pub <- z(h$hi31) + z(h$hi32) + z(h$hipubsoc); oth <- z(h$hi33) + z(h$hilabour) + z(h$hicapital); tot <- pub + oth
  wt <- if (all(is.na(h$hwgt))) rep(1, nrow(h)) else h$hwgt
  old <- h$age >= a0 & !is.na(h$age); hb <- z(h$hi46) > 0; ro <- old & tot > 0 & oth <= .1 * tot; r2 <- ro & z(h$hi32) <= .1 * tot
  fr <- function(i, j) round(sum(wt[i & j]) / sum(wt[i]), 3)
  cat("LINES", nm, "median_equiv", round(med), "L50", round(L50), "L60", round(L60), "old<L50", fr(old, h$dhi < L50), "old<L60", fr(old, h$dhi < L60),
      "r2HB<L50", fr(r2 & hb, h$dhi < L50), "r2noHB<L50", fr(r2 & !hb, h$dhi < L50), "r2HB<L60", fr(r2 & hb, h$dhi < L60), "r2noHB<L60", fr(r2 & !hb, h$dhi < L60), "\n")
  cat("HB", nm, "reliant_share", fr(old, ro), "hb_among_reliant", fr(ro, hb), "med_hi46_recip", round(median(h$hi46[ro & hb])), "\n")
  if (cc == "dk") {
    pc(paste("SHIFT", nm, "r2 off2500"), cell_test_base(h$dhi[r2], wt[r2], 5000, qlo = .05, off = 2500))
    pc(paste("SHIFT", nm, "old off2500"), cell_test_base(h$dhi[old], wt[old], 5000, qlo = .05, off = 2500))
    set.seed(7); io <- which(old & !ro); ir <- which(ro); k <- round(.06 / .94 * length(io))
    pv <- vapply(1:10, function(b) { s <- c(io, sample(ir, k)); cell_test_base(h$dhi[s], wt[s], 5000, seed = b, qlo = .05)[["p"]] }, 0)
    cat("REWEIGHT", nm, "reliant share .06, 10 reps: rej05", mean(pv <= .05), "median p", median(pv), "\n")
    i2 <- which(r2); set.seed(11)
    for (n in c(150, 300, 700)) { pv <- vapply(1:20, function(b) { s <- sample(i2, n, prob = wt[i2]); cell_test_base(h$dhi[s], rep(1, n), 5000, seed = b, qlo = .05)[["p"]] }, 0)
      cat("POWER", nm, "r2 n", n, "20 reps rej05", mean(pv <= .05), "rej10", mean(pv <= .10), "\n") }
  }
}
