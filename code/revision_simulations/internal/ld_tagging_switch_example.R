# A 16-day example in which LD tagging creates an apparent sign switch.

set.seed(20260731)
n <- 10000L
days <- 0:15

# Generate two haplotypes per individual. A and B are approximately
# independent, while the non-causal tag SNP C is in positive LD with A and
# negative LD with B.
h_a <- rbinom(2L * n, size = 1L, prob = 0.5)
h_b <- rbinom(2L * n, size = 1L, prob = 0.5)
p_c <- 0.5 + 0.45 * (h_a - h_b)
h_c <- rbinom(2L * n, size = 1L, prob = p_c)
combine_haplotypes <- function(x) x[seq_len(n)] + x[n + seq_len(n)]

genotype <- data.frame(
  A = combine_haplotypes(h_a),
  B = combine_haplotypes(h_b),
  C = combine_haplotypes(h_c)
)

# A has a positive early effect and B has a positive late effect. C has no
# direct causal effect at any time.
true_effect_a <- exp(-0.5 * ((days - 2) / 2.5)^2)
true_effect_b <- exp(-0.5 * ((days - 13) / 2.5)^2)
expression <- outer(genotype$A, true_effect_a) +
  outer(genotype$B, true_effect_b) +
  matrix(rnorm(n * length(days)), nrow = n)

fit_c_effect <- function(day_index, conditional = FALSE) {
  response <- expression[, day_index]
  fit <- if (conditional) {
    lm(response ~ C + A + B, data = genotype)
  } else {
    lm(response ~ C, data = genotype)
  }
  c(
    estimate = coef(fit)[["C"]],
    standard_error = sqrt(vcov(fit)["C", "C"])
  )
}

marginal_fit_c <- t(vapply(seq_along(days), fit_c_effect, numeric(2)))
conditional_fit_c <- t(vapply(
  seq_along(days), fit_c_effect, numeric(2), conditional = TRUE
))

marginal_effect_c <- marginal_fit_c[, "estimate"]
marginal_se_c <- marginal_fit_c[, "standard_error"]
conditional_effect_c <- conditional_fit_c[, "estimate"]
conditional_se_c <- conditional_fit_c[, "standard_error"]

result <- data.frame(
  day = days,
  true_effect_a = true_effect_a,
  true_effect_b = true_effect_b,
  marginal_effect_c = marginal_effect_c,
  marginal_se_c = marginal_se_c,
  conditional_effect_c = conditional_effect_c,
  conditional_se_c = conditional_se_c
)

stopifnot(
  marginal_effect_c[which.max(true_effect_a)] > 0,
  marginal_effect_c[which.max(true_effect_b)] < 0,
  max(abs(conditional_effect_c)) < 0.1,
  all(marginal_se_c > 0),
  all(conditional_se_c > 0)
)

print(round(cor(genotype), 3))
print(result, digits = 3, row.names = FALSE)

variant_colors <- c(C = "#2a78d6", A = "#eb6834", B = "#1baf7a")
output_file <- file.path(
  "output", "revision_simulations", "internal", "ld_tagging_switch_example.png"
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

png(output_file, width = 1400, height = 600, res = 150)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.5, 1.0))

plot(days, true_effect_a,
  type = "l", lwd = 3, col = variant_colors[["A"]],
  ylim = c(0, 1.05), xlab = "Day", ylab = "True causal effect",
  main = "Distinct early and late causal variants"
)
lines(days, true_effect_b, lwd = 3, lty = 2, col = variant_colors[["B"]])
lines(days, rep(0, length(days)), lwd = 2, lty = 3, col = "#777777")
legend("topright",
  legend = c("A: early causal", "B: late causal", "C: no direct effect"),
  col = c(variant_colors[["A"]], variant_colors[["B"]], "#777777"),
  lty = c(1, 2, 3), lwd = c(3, 3, 2), bty = "n"
)

y_limit <- 1.15 * max(abs(c(
  marginal_effect_c - marginal_se_c,
  marginal_effect_c + marginal_se_c,
  conditional_effect_c - conditional_se_c,
  conditional_effect_c + conditional_se_c
)))
plot(days, marginal_effect_c,
  type = "n",
  ylim = c(-y_limit, y_limit), xlab = "Day", ylab = "Effect estimate for C",
  main = "Tag SNP C appears to switch"
)
abline(h = 0, col = "#c3c2b7")
arrows(
  days, marginal_effect_c - marginal_se_c,
  days, marginal_effect_c + marginal_se_c,
  angle = 90, code = 3, length = 0.035, lwd = 1.5,
  col = adjustcolor(variant_colors[["C"]], alpha.f = 0.75)
)
arrows(
  days, conditional_effect_c - conditional_se_c,
  days, conditional_effect_c + conditional_se_c,
  angle = 90, code = 3, length = 0.035, lwd = 1.25,
  col = adjustcolor("#777777", alpha.f = 0.75)
)
lines(days, marginal_effect_c, lwd = 3, col = variant_colors[["C"]])
points(days, marginal_effect_c, pch = 16, col = variant_colors[["C"]])
lines(days, conditional_effect_c, lwd = 2, lty = 2, col = "#777777")
points(days, conditional_effect_c, pch = 1, col = "#777777")
legend("topright",
  legend = c("Marginal C", "C conditional on A and B"),
  col = c(variant_colors[["C"]], "#777777"),
  pch = c(16, 1), lty = c(1, 2), lwd = c(3, 2), bty = "n"
)

dev.off()
cat("\nSaved figure to", output_file, "\n")
