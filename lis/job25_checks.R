## Glass-Barriers job 25: per-wave checks, alternative scale, France early waves, Danish refined parts. Output: test statistics and shares only.
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
sh <- function(x, wt, a, b) round(sum(wt[x >= a & x < b]) / sum(wt), 4)
for (cc in c("nl", "au", "uk", "fr", "dk")) { X <- NULL
  yrs <- if (cc == "dk") c(2016, 2019) else 2010:2021
  for (yr in yrs) {
    hh <- try(suppressMessages(lissyuse(data = cc, vars = hv, from = yr, to = yr)), silent = TRUE); nm <- paste0(cc, substr(yr, 3, 4))
    if (inherits(hh, "try-error") || is.null(hh[[paste0(nm, "h")]])) next
    pp <- try(suppressMessages(lissyuse(data = cc, vars = c("age"), from = yr, to = yr)), silent = TRUE)
    if (inherits(pp, "try-error") || is.null(pp[[paste0(nm, "p")]])) next
    h <- hh[[paste0(nm, "h")]]; p <- pp[[paste0(nm, "p")]]; rm(hh, pp)
    h <- h[!is.na(h$nhhmem) & h$nhhmem == 1 & !is.na(h$dhi) & h$dhi > 0, ]; h$age <- p$age[match(h$hid, p$hid)]; rm(p)
    h <- h[h$age >= 65 & !is.na(h$age), ]
    pub <- z(h$hi31) + z(h$hi32) + z(h$hipubsoc); oth <- z(h$hi33) + z(h$hilabour) + z(h$hicapital); tot <- pub + oth
    ro <- tot > 0 & oth <= .1 * tot; hb <- z(h$hi46) > 0; r2 <- ro & z(h$hi32) <= .1 * tot
    s0 <- median(h$dhi[ro & !hb]); s1 <- median(h$dhi[ro]); wt <- if (all(is.na(h$hwgt))) rep(1, nrow(h)) else h$hwgt; wt <- wt / mean(wt)
    d <- data.frame(yr = yr, x = 100 * h$dhi / s0, xa = 100 * h$dhi / s1, x1 = 100 * (h$dhi - z(h$hi46)) / s0, hb = hb, r2 = r2, wt = wt)[ro, ]
    if (cc %in% c("uk", "au", "dk")) pc(paste("WAVE", nm, "n", nrow(d), "w4"), cell_test_base(d$x, d$wt, 4, qlo = .05))
    if (cc == "nl") cat("NLWAVE", nm, "n", nrow(d), "valley_103_115", sh(d$x, d$wt, 103, 115), "left_91_103", sh(d$x, d$wt, 91, 103), "right_115_127", sh(d$x, d$wt, 115, 127), "\n")
    X <- rbind(X, d); rm(h) }
  if (cc %in% c("nl", "au", "uk", "dk")) pc(paste("ALTSCALE", cc, "w2"), cell_test_base(X$xa, X$wt, 2, qlo = .05))
  if (cc == "fr") { i <- X$yr <= 2014; cat("FR early n", sum(i), "\n")
    pc("CELL fr1014 w2", cell_test_base(X$x[i], X$wt[i], 2, qlo = .05)); pc("CELL fr1014 w4", cell_test_base(X$x[i], X$wt[i], 4, qlo = .05))
    j <- i & X$x1 > 0; pc("CELL fr1014 minusHB w2", cell_test_base(X$x1[j], X$wt[j], 2, qlo = .05)) }
  if (cc == "dk") { i <- X$r2 & X$x1 > 0; pc("CELL dk refined minusHB w2", cell_test_base(X$x1[i], X$wt[i], 2, qlo = .05))
    i <- X$r2 & !X$hb; pc("CELL dk refined noHB w2", cell_test_base(X$x[i], X$wt[i], 2, qlo = .05))
    i <- X$r2 & X$hb; pc("CELL dk refined HB w2", cell_test_base(X$x[i], X$wt[i], 2, qlo = .05)) }
}
