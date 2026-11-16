# Reference results from vegan::adonis2 for test/vegan/compare.jl.
# Usage: Rscript reference.R OUTDIR
suppressMessages(library(vegan))
out <- commandArgs(TRUE)[1]
set.seed(1)
n <- 80
dat <- data.frame(a = sample(c("x", "y", "z"), n, TRUE), b = sample(c("u", "v"), n, TRUE), c = rnorm(n))
Y <- matrix(rpois(n * 10, 5), n, 10)
Y[dat$a == "x", 1:3] <- Y[dat$a == "x", 1:3] + 4
D <- vegdist(Y, "bray")
perm <- shuffleSet(n, 999)
w <- function(x, f) write.csv(x, file.path(out, f), row.names = FALSE)
w(dat, "data.csv"); w(Y, "Y.csv"); w(unclass(perm), "perm.csv")
tab <- function(r) data.frame(term = rownames(r), Df = r$Df, SumOfSqs = r$SumOfSqs, R2 = r$R2, F = r$F, P = r[["Pr(>F)"]])
w(tab(adonis2(D ~ a + b + c + a:b, dat, permutations = perm, by = "terms")), "terms.csv")
w(tab(adonis2(D ~ a + b + c, dat, permutations = perm, by = "margin")), "margin.csv")
w(tab(adonis2(D ~ a + b + c, dat, permutations = perm)), "model.csv")
cat("vegan", as.character(packageVersion("vegan")), "\n")
