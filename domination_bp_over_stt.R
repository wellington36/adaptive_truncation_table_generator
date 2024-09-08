library(rstan)

expose_stan_functions("stan/domination_bp_over_stt.stan")

a <- c(2, 1.1, 1.01, 1.001, 1.0001)

normalising_btt_terms <- numeric(length(a))
normalising_bp_terms <- numeric(length(a))
percentage <- numeric(length(a))
for (i in seq_along(a)) {
  normalising_btt_terms[i] <- normalising_f_btt(a[i], 2^-52);
  normalising_bp_terms[i] <- normalising_f_bp(a[i], 2^-52);
  percentage[i] = normalising_bp_terms[i]/normalising_btt_terms[i] - 1;
}

############
table <- data.frame(
  Column1 = normalising_btt_terms,
  Column2 = normalising_bp_terms,
  Column3 = percentage
)

# Set the column names
colnames(table) <- c("Sum-to-threshold", "Bounding-pairs", "Percentage")

# Set the row names
rownames(table) <- c("a=2", "a=1.1", "a=1.01", "a=1.001", "a=1.0001")

# Print the table
print(table)