# Load packages
library(epitools)

# Read file
file1 <- read.csv("Tcon vs Treg pvalue adn cohensd statistical test.csv")

# Rename columns to cleaner names
colnames(file1) <- c("subsystem", "Pos_d_lt0", "Pos_d_gt0", "Neg_d_lt0", "Neg_d_gt0")

# Replace NA with 0
file1[is.na(file1)] <- 0

# Initialize results list
results_list <- list()

# Loop over rows
for (i in 1:nrow(file1)) {
  
  # Extract counts
  a <- as.numeric(file1$Pos_d_gt0[i])  # p<0.05_Pos_d>0
  b <- as.numeric(file1$Neg_d_gt0[i])  # p<0.05_Neg_d>0
  c <- as.numeric(file1$Pos_d_lt0[i])  # p<0.05_Pos_d<0
  d <- as.numeric(file1$Neg_d_lt0[i])  # p<0.05_Neg_d<0
  
  # Build raw matrix
  mat <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE,
                dimnames = list(c("Pos", "Neg"), c("d>0", "d<0")))
  
  # Fisher Test (raw table)
  fisher_res <- fisher.test(mat)
  
  # Haldane–Anscombe correction
  mat_corr <- mat + 0.5
  chisq_res <- chisq.test(mat_corr)
  or_res <- oddsratio(mat_corr, method="wald")
  
  # Store results
  results_list[[i]] <- data.frame(
    Subsystem = file1$subsystem[i],
    Count_a = a, Count_b = b, Count_c = c, Count_d = d,
    Fisher_p = fisher_res$p.value,
    ChiSquare_statistic = unname(chisq_res$statistic),
    ChiSquare_p = chisq_res$p.value,
    Odds_Ratio = or_res$measure[2,1],
    CI_Low = or_res$measure[2,2],
    CI_High = or_res$measure[2,3]
  )
}

# Combine all into one dataframe
results <- do.call(rbind, results_list)

# Save to CSV
write.csv(results, "Statistical_test Haldane applied to all.csv", row.names = FALSE)

cat("✅ Results saved to 'Statistical_test_results_per_subsystem.csv'\n")
