## Glass-Barriers job 17 (Denmark): binned distribution for a figure (5,000 DKK cells, cells with fewer than 30 observations suppressed). Aggregate shares only.
suppressPackageStartupMessages({ library(lissyrtools) })
hv <- c("dhi", "hipubsoc", "hi31", "hi32", "hi33", "hi46", "hilabour", "hicapital", "nhhmem", "hwgt")
z <- function(x) ifelse(is.na(x), 0, x)
for (nm in c("dk10", "dk19")) { yr <- 2000 + as.integer(substr(nm, 3, 4))
  hh <- lissyuse(data = "dk", vars = hv, from = yr, to = yr); pp <- lissyuse(data = "dk", vars = c("age"), from = yr, to = yr)
  h <- hh[[paste0(nm, "h")]]; p <- pp[[paste0(nm, "p")]]; rm(hh, pp)
  h <- h[!is.na(h$nhhmem) & h$nhhmem == 1 & !is.na(h$dhi) & h$dhi > 0, ]; h$age <- p$age[match(h$hid, p$hid)]; rm(p)
  pub <- z(h$hi31) + z(h$hi32) + z(h$hipubsoc); oth <- z(h$hi33) + z(h$hilabour) + z(h$hicapital); tot <- pub + oth
  wt <- if (all(is.na(h$hwgt))) rep(1, nrow(h)) else h$hwgt
  old <- h$age >= 65 & !is.na(h$age); hb <- z(h$hi46) > 0; r2 <- old & tot > 0 & oth <= .1 * tot & z(h$hi32) <= .1 * tot
  br <- seq(47500, 302500, by = 5000); ctr <- (br[-1] + br[-length(br)]) / 2
  for (g in c("old", "r2_hb", "r2_nohb", "old_other")) {
    i <- switch(g, old = old, r2_hb = r2 & hb, r2_nohb = r2 & !hb, old_other = old & !r2)
    x <- h$dhi[i]; v <- wt[i]; k <- cut(x, br, labels = FALSE)
    cnt <- tabulate(k, nbins = length(ctr)); sh <- tapply(v, factor(k, levels = seq_along(ctr)), sum); sh[is.na(sh)] <- 0; sh <- sh / sum(v)
    sh[cnt < 30] <- NA
    cat("FIG", nm, g, "n", sum(i), "share_of_old", round(sum(v) / sum(wt[old]), 3), "\n")
    cat("FIGSHARES", nm, g, paste(ctr / 1000, round(sh, 4), sep = ":", collapse = " "), "\n") }
}
