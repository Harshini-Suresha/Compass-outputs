library(tidyverse)

# Load data
df <- read.csv("Tcon vs exTreg pvalue adn cohensd statistical test.csv", check.names = FALSE)

# Ensure numeric conversion and NA handling
df <- df %>%
  mutate(
    `p<0.05_Pos_d<0` = as.numeric(`p<0.05_Pos_d<0`),
    `p<0.05_Neg_d<0` = as.numeric(`p<0.05_Neg_d<0`),
    `p<0.05_Pos_d>0` = as.numeric(`p<0.05_Pos_d>0`),
    `p<0.05_Neg_d>0` = as.numeric(`p<0.05_Neg_d>0`)
  ) %>%
  replace_na(list(
    `p<0.05_Pos_d<0` = 0,
    `p<0.05_Neg_d<0` = 0,
    `p<0.05_Pos_d>0` = 0,
    `p<0.05_Neg_d>0` = 0
  ))

# Statistical tests and OR
results <- df %>%
  rowwise() %>%
  mutate(
    a = ifelse(`p<0.05_Pos_d<0` == 0, 0.0001, `p<0.05_Pos_d<0`),
    b = ifelse(`p<0.05_Neg_d<0` == 0, 0.0001, `p<0.05_Neg_d<0`),
    c = ifelse(`p<0.05_Pos_d>0` == 0, 0.0001, `p<0.05_Pos_d>0`),
    d = ifelse(`p<0.05_Neg_d>0` == 0, 0.0001, `p<0.05_Neg_d>0`),
    
    # Chi-square test
    chi_p = tryCatch({
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      suppressWarnings(chisq.test(tbl, correct = FALSE)$p.value)
    }, error = function(e) NA_real_()),
    
    # Fisher's exact test
    fisher_p = tryCatch({
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      fisher.test(tbl)$p.value
    }, error = function(e) NA_real_()),
    
    # Odds Ratio and CI
    odds_ratio = tryCatch({
      (a * d) / (b * c)
    }, error = function(e) NA_real_()),
    
    se_log_or = tryCatch({
      sqrt(1/a + 1/b + 1/c + 1/d)
    }, error = function(e) NA_real_()),
    
    log_or = log(odds_ratio),
    ci_low = exp(log_or - 1.96 * se_log_or),
    ci_high = exp(log_or + 1.96 * se_log_or)
  ) %>%
  ungroup() %>%
  select(-se_log_or, -log_or)

# Save to file
write.csv(results, "Normalised 0.0001.csv", row.names = FALSE)

# Preview result
head(results)
