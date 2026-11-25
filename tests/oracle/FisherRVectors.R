#!/usr/bin/env Rscript
# Writes R fisher.test vectors for src/stats/fisher (FisherFlavor::r): the 2x2 tables of
# tests/golden/fisher_vectors.tsv (the jar's vectors), then 20,000 random tables shaped like
# variant strand tables (fixed seed). The matrix is built as the scripts in tests/scripts/ build
# it: matrix(c(a, b, c, d), nrow = 2). The random set exists because uniroot's stopping test is a
# threshold: function values that differ in the last bits can stop it one step apart, which
# moves the 5th decimal of the estimate (found on a real panel; fixed by summing in long double
# as R's sum() does). A third block (10,000 tables, deep coverage against 0-2 reads in the other
# cells) holds one-sided p-values that sit on a 5-dp tie such as 41/64, where R's phyper()
# series and a plain sum of densities round differently (found on the HG008 tumour/normal panel).
#
# Columns: a b c d, then, each rounded to 5 dp and printed as write.table prints it under
# options(scipen = 999) (15 significant digits):
#   two-sided p, conditional-MLE estimate   (teststrandbias.R, testsomatic.R strand tests)
#   p and estimate for alternative="greater", p and estimate for alternative="less"
#                                            (testsomatic.R tumour-vs-normal test)
#
# Usage: Rscript tests/oracle/FisherRVectors.R > tests/golden/fisher_r_vectors.tsv
tables <- read.table("tests/golden/fisher_vectors.tsv", sep = "\t", colClasses = "character")[, 1:4]
options(scipen = 999)
fmt <- function(v) format(round(v, 5), digits = 15)
set.seed(20261101)
n <- 20000
alt_max <- sample(c(3, 10, 40, 200), n, replace = TRUE)
ref_max <- sample(c(20, 200, 1000, 3000), n, replace = TRUE)
extra <- cbind(sapply(ref_max, function(r) sample.int(r + 1, 1) - 1), sapply(ref_max, function(r) sample.int(r + 1, 1) - 1),
               sapply(alt_max, function(a) sample.int(a + 1, 1) - 1), sapply(alt_max, function(a) sample.int(a + 1, 1) - 1))
n3 <- 10000
deep <- cbind(sample.int(401, n3, replace = TRUE) - 1, sample.int(3, n3, replace = TRUE) - 1,
              sample.int(401, n3, replace = TRUE) - 1, sample.int(3, n3, replace = TRUE) - 1)
all_tables <- rbind(as.matrix(sapply(tables, as.numeric)), extra, deep)
for (i in seq_len(nrow(all_tables))) {
    t <- as.numeric(all_tables[i, ])
    m <- matrix(t, nrow = 2)
    h <- fisher.test(m)
    g <- fisher.test(m, alternative = "greater")
    l <- fisher.test(m, alternative = "less")
    cat(paste(c(t, fmt(h$p.value), fmt(h$estimate), fmt(g$p.value), fmt(g$estimate),
                fmt(l$p.value), fmt(l$estimate)), collapse = "\t"), "\n", sep = "")
}
