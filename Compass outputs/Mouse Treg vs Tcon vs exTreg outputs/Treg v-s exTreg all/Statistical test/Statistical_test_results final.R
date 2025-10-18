# Load packages
library(epitools)

# Read file
file1 <- read.csv("Tregs vs exTregs pvalue adn cohensd statistical test.csv")

# Fix column names
colnames(file1) <- c("subsystem", "Pos_d_lt0", "Pos_d_gt0", "Neg_d_lt0", "Neg_d_gt0")

# Ensure numeric
file1$Pos_d_lt0 <- as.numeric(file1$Pos_d_lt0)
file1$Pos_d_gt0 <- as.numeric(file1$Pos_d_gt0)
file1$Neg_d_lt0 <- as.numeric(file1$Neg_d_lt0)
file1$Neg_d_gt0 <- as.numeric(file1$Neg_d_gt0)

results_list <- vector("list", nrow(file1))

for (i in seq_len(nrow(file1))) {
  a <- file1$Pos_d_gt0[i]  # p<0.05_Pos_d>0
  b <- file1$Neg_d_gt0[i]  # p<0.05_Neg_d>0
  c <- file1$Pos_d_lt0[i]  # p<0.05_Pos_d<0
  d <- file1$Neg_d_lt0[i]  # p<0.05_Neg_d<0
  
  # Build 2x2 matrix
  mat <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE,
                dimnames = list(c("Pos","Neg"), c("d>0","d<0")))
  
  # Fisher’s test always on raw counts
  fisher_res <- fisher.test(mat)
  
  # Decide whether to apply Haldane correction
  if (any(mat == 0)) {
    mat_corr <- mat + 0.5
    used_correction <- TRUE
  } else {
    mat_corr <- mat
    used_correction <- FALSE
  }
  
  # Chi-square test
  chisq_res <- suppressWarnings(chisq.test(mat_corr))
  
  # Odds ratio + CI
  or_res <- oddsratio(mat_corr, method="wald")
  
  results_list[[i]] <- data.frame(
    Subsystem = file1$subsystem[i],
    Count_a = a, Count_b = b, Count_c = c, Count_d = d,
    Fisher_p = fisher_res$p.value,
    ChiSquare_statistic = unname(chisq_res$statistic),
    ChiSquare_p = chisq_res$p.value,
    Odds_Ratio = or_res$measure[2,1],
    CI_Low = or_res$measure[2,2],
    CI_High = or_res$measure[2,3],
    HaldaneApplied = used_correction,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, results_list)

# Save output
write.csv(results, "Statistical_test_results.csv", row.names = FALSE)
cat("✅ Results saved to 'Statistical_test_results_per_subsystem.csv'\n")
